import tqdm
import numpy as np
import scipy

def data_matrix(p3_basis, uvtn, idptr, v_num, F, edge_global_info):
    def assemble_matrix_single(i, idptr, basis, uv, v_num, F, edge_global_info):
        global_id = global_id_for_face(
            i, v_num, F, edge_global_info)
        x, y = uv[:, 0], uv[:, 1]  # uvtn[idptr:idptr_1].T
        cols = global_id.repeat(len(x))
        data = np.concatenate([f(x, y) for f in basis])
        rows = np.tile(np.arange(idptr[0], idptr[1]), 10)
        return rows, cols, data
    a_rows, a_cols, a_data, rhs_t = [], [], [], []
    for i in range(F.shape[0]):
        r, c, d = assemble_matrix_single(i, (idptr[i:i+2]), basis=p3_basis, uv=uvtn[idptr[i]:idptr[i+1], :2],
                                         v_num=v_num, F=F, edge_global_info=edge_global_info)
        a_rows.append(r)
        a_cols.append(c)
        a_data.append(d)
    f_num = F.shape[0]
    e_num = edge_global_info[0].max() + 1
    A = scipy.sparse.coo_matrix((np.concatenate(a_data),
                                 (np.concatenate(a_rows).astype(np.int), np.concatenate(a_cols).astype(np.int))),
                                shape=(len(uvtn), v_num+f_num + 2*e_num))
    return A.tocsr()

def offset_from_ctrlpt_operator(value_at_points, all_global_ids, SVI, SVJ, all_dir, numcp, f_num):
    '''build Sparse matrix A such that 
        A*sol = h*D, on upsampled points, and height*Dir is the offset.
        Basically, this is differentiation against {newbase + ((sol[all_global_ids]@values_at_points).reshape(-1,1)[SVI]*new_dir}'''
    dh_s_r, dh_s_c, dh_s_d = [], [], []
    newdir = all_dir.reshape(-1, 3)[SVI]
    for i in tqdm.trange(numcp):
        nzfaces, nzids = np.where(all_global_ids == i)
        h_nz = value_at_points[nzids]
        srows = SVJ[np.ravel_multi_index([np.repeat(nzfaces, value_at_points.shape[1]), np.tile(np.arange(value_at_points.shape[1]), len(nzfaces))],
                                         (f_num, value_at_points.shape[1]))]
        h_nz_data = h_nz.flatten()
        dh_s_id = np.ravel_multi_index([np.tile(srows, (3, 1)).T.flatten(),  # [a,a,a,b,b,b ....; 0,1,2,0,1,2]
                                        np.tile(np.arange(3), (srows.shape[0], 1)).flatten()], (SVI.shape[0], 3))

        dh_s_id, uid = np.unique(dh_s_id, return_index=True)
        cols_i = np.ones(len(dh_s_id))*i
        data_i = ((newdir[srows])*h_nz_data[:, None]).flatten()[uid]
        dh_s_r.append(dh_s_id)
        dh_s_c.append(cols_i)
        dh_s_d.append(data_i)
    sp_dh_D_s = scipy.sparse.coo_matrix((np.concatenate(dh_s_d),
                                         (np.concatenate(dh_s_r), np.concatenate(dh_s_c))),
                                        shape=(newdir.shape[0]*3, numcp)).tocsr()
    return sp_dh_D_s


def spline_operator(value_at_points, all_global_ids, SVI, SVJ, all_dir, numcp, f_num):
    '''build Sparse matrix A such that 
        A*ctrl = V
        Basically, this is differentiation against ctrl[all_global_ids]@value_at_points '''
    assert value_at_points.shape[0]==all_global_ids.shape[1]
    dh_s_r, dh_s_c, dh_s_d = [], [], []
    for i in tqdm.trange(numcp):
        nzfaces, nzids = np.where(all_global_ids == i)
        h_nz = value_at_points[nzids]
        srows = SVJ[np.ravel_multi_index([np.repeat(nzfaces, value_at_points.shape[1]), np.tile(np.arange(value_at_points.shape[1]), len(nzfaces))],
                                         (f_num, value_at_points.shape[1]))]
        h_nz_data = h_nz.flatten()
        dh_s_id = srows

        dh_s_id, uid = np.unique(dh_s_id, return_index=True)
        cols_i = np.ones(len(dh_s_id))*i
        data_i = (h_nz_data[:, None]).flatten()[uid]
        dh_s_r.append(dh_s_id)
        dh_s_c.append(cols_i)
        dh_s_d.append(data_i)
    sp_dh_D_s = scipy.sparse.coo_matrix((np.concatenate(dh_s_d),
                                         (np.concatenate(dh_s_r), np.concatenate(dh_s_c))),
                                        shape=(len(SVI), numcp)).tocsr()
    return sp_dh_D_s


def perpendicular_direction_op(V, F, TT,TTi, edge_global_info, all_global_ids, basis, basis_du, basis_dv):
    '''
    Matrix A such that A*ctrl = Dir
    where Dir is the dr/dv, if the control point is categorized to du.
    The operators are used for approximate continuity optimizations.
    '''
    rows, cols = [], []
    opp_cols = []
    along_d, ortho_d, oppo_d = [], [], []
    pos_d = []
    edge_index, edge_flag = edge_global_info

    # upsample and take values.
    level = 6
    sample_u = [np.arange(level)/(level-1), 1 -
                np.arange(level)/(level-1), np.zeros(level)]
    sample_v = [np.zeros(level), np.arange(
        level)/(level-1), 1-np.arange(level)/(level-1)]
    values_du = [np.vstack(basis_du(sample_u[j], sample_v[j])).T
                 for j in range(3)]
    values_dv = [np.vstack(basis_dv(sample_u[j], sample_v[j])).T
                 for j in range(3)]
    values = [np.vstack(basis(sample_u[j], sample_v[j])).T for j in range(3)]
    along_dir = np.array([[1, 0], [-1/np.sqrt(2), 1/np.sqrt(2)], [0, -1]])
    ortho_dir = np.array([[0, 1], [-1/np.sqrt(2), -1/np.sqrt(2)], [1, 0]])

    cur = 0
    for i, f in enumerate(F):
        cur_ids = all_global_ids[i]#, len(V), F, edge_global_info)
        for j in range(3):
            if edge_flag[i, j] == 0: # only do correct side
                continue
            ii, jj = TT[i,j], TTi[i,j]
            if ii == -1: # skip boundary edges
                continue
            opp_ids = all_global_ids[ii]
            # rows: adding new, cols: to ctrl, and data is the basis value there.
            ortho_data = ortho_dir[j, 0]*values_du[j] + ortho_dir[j, 1]*values_dv[j]
            along_data = along_dir[j, 0]*values_du[j] + along_dir[j, 1]*values_dv[j]
            oppo_data = (ortho_dir[jj, 0]*values_du[jj] + ortho_dir[jj, 1]*values_dv[jj])[::-1]

            rows.append(np.repeat(np.arange(cur, cur+level), len(cur_ids)))
            cur = cur + level
            cols.append(np.tile(cur_ids, level))
            opp_cols.append(np.tile(opp_ids, level))
            ortho_d.append(ortho_data.flatten())
            along_d.append(along_data.flatten())
            oppo_d.append(oppo_data.flatten())
            pos_d.append(values[j].flatten())

    rows, cols = np.concatenate(rows), np.concatenate(cols)
    numcp = len(V)+ len(F) + 2*(edge_index.max()+1)
    ortho_mat = scipy.sparse.coo_matrix((np.concatenate(ortho_d),
                                         (rows, cols)), shape=(cur, numcp))
    oppo_mat = scipy.sparse.coo_matrix((np.concatenate(oppo_d),
                                         (rows, np.concatenate(opp_cols))), shape=(cur, numcp))
    along_mat = scipy.sparse.coo_matrix((np.concatenate(along_d),
                                         (rows, cols)), shape=(cur, numcp))
    pos_mat = scipy.sparse.coo_matrix((np.concatenate(pos_d),
                                         (rows, cols)), shape=(cur, numcp))
    return along_mat, ortho_mat, oppo_mat, pos_mat

