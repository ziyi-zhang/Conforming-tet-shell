import numpy as np
import torch
import igl
import math


def edge_globals(F):
    """
    edge_global_index: #F*3, each corner to a global edge index
    edge_global_flag: #F*3, default 0 opposite, 1 means true (half-)edge.
    """
    FF, FFi = igl.triangle_triangle_adjacency(F)
    edge_global_index, edge_global_flag = -np.ones_like(F), np.zeros_like(F)
    cnt = 0
    for f in range(F.shape[0]):
        for e in range(3):
            if edge_global_index[f][e] >= 0:
                continue  # already assigned, pass
            edge_global_index[f, e] = cnt
            edge_global_flag[f, e] = 1
            if FF[f, e] != -1:  # ignore border
                edge_global_index[FF[f, e], FFi[f, e]] = cnt
            cnt += 1
    return edge_global_index, edge_global_flag


def global_id_for_face_cubic(i, v_num, F, edge_info):
    '''Note that this is a specific ordering for cubic (10 nodes)'''
    edge_global_index, edge_global_flag = edge_info
    global_id = np.zeros(10, dtype=np.int)
    global_id[:3] = F[i]  # local 0,1,2
    for e in range(3):  # local 2*e + (3,4)
        global_id[2*e+3] = v_num+F.shape[0] + 2*edge_global_index[i, e]
        global_id[2*e+4] = v_num+F.shape[0] + 2*edge_global_index[i, e] + 1
        if edge_global_flag[i, e] == 0:  # swap two
            global_id[2*e+3], global_id[2*e +
                                        4] = global_id[2*e+4], global_id[2*e+3]
    global_id[-1] = v_num + i  # local 9
    return global_id


def global_id_for_face_quadratic(i, v_num, F, edge_info):
    '''Note that this is a specific ordering for quadratic (6 nodes)'''
    edge_global_index, edge_global_flag = edge_info
    global_id = np.zeros(6, dtype=np.int)
    global_id[:3] = F[i]  # local 0,1,2
    for e in range(3):  # local e + 3
        global_id[e+3] = v_num + edge_global_index[i, e]
    return global_id


def lagrange_basis():
    p3_basis_dx = [lambda x, y:-13.5*x**2 - 27.0*x*y + 18*x - 13.5*y**2 + 18*y - 5.5,
                   lambda x, y:4.5*x**2 + 0.5*x*(18*x - 9) - 4.5*x + 1.0,
                   lambda x, y: 0*x,
                   lambda x, y:13.5*x *
                   (x + y - 1) + 4.5*x*(3*x + 3*y - 2) +
                   (3*x + 3*y - 2)*(4.5*x + 4.5*y - 4.5),
                   lambda x, y:-13.5*x**2 - 13.5*x*y - 4.5*x *
                   (6*x + 3*y - 4) + 18.0*x + 4.5*y - 4.5,
                   lambda x, y:13.5*x*y + 4.5*y*(3*x - 1),
                   lambda x, y:4.5*y*(3*y - 1),
                   lambda x, y:-4.5*y*(3*y - 1),
                   lambda x, y:13.5*y*(x + y - 1) + 4.5*y*(3*x + 3*y - 2),
                   lambda x, y:-27*x*y - 27*y*(x + y - 1)]
    p3_basis_dy = [lambda x, y:-13.5*x**2 - 27.0*x*y + 18*x - 13.5*y**2 + 18*y - 5.5,
                   lambda x, y:0*x,
                   lambda x, y:4.5*y**2 + 0.5*y*(18*y - 9) - 4.5*y + 1.0,
                   lambda x, y:13.5*x*(x + y - 1) + 4.5*x*(3*x + 3*y - 2),
                   lambda x, y:-4.5*x*(3*x - 1),
                   lambda x, y:4.5*x*(3*x - 1),
                   lambda x, y:13.5*x*y + 4.5*x*(3*y - 1),
                   lambda x, y:-13.5*x*y + 4.5*x - 13.5*y**2 -
                   4.5*y*(3*x + 6*y - 4) + 18.0*y - 4.5,
                   lambda x, y:13.5*y *
                   (x + y - 1) + 4.5*y*(3*x + 3*y - 2) +
                   (3*x + 3*y - 2)*(4.5*x + 4.5*y - 4.5),
                   lambda x, y:-27*x*y - 27*x*(x + y - 1)]
    p3_basis = [
        lambda x, y: -27.0/2.0*x**2*y + 9*x**2 - 27.0/2.0*y**2*x + 9*y**2 -
        9.0/2.0*x**3 + 18*x*y - 11.0/2.0*x - 9.0/2.0*y**3 - 11.0/2.0*y + 1,
        lambda x, y: (1.0/2.0)*x*(9*x**2 - 9*x + 2),
        lambda x, y: (1.0/2.0)*y*(9*y**2 - 9*y + 2),
        lambda x, y: (9.0/2.0)*x*(x + y - 1)*(3*x + 3*y - 2),
        lambda x, y: -9.0/2.0*x*(3*x**2 + 3*x*y - 4*x - y + 1),
        lambda x, y: (9.0/2.0)*x*y*(3*x - 1),
        lambda x, y: (9.0/2.0)*x*y*(3*y - 1),
        lambda x, y: -9.0/2.0*y*(3*x*y - x + 3*y**2 - 4*y + 1),
        lambda x, y: (9.0/2.0)*y*(x + y - 1)*(3*x + 3*y - 2),
        lambda x, y: -27*x*y*(x + y - 1)]
    for x, y in [(0, 0), (1, 0), (0, 1), (1/3, 0), (2/3, 0), (2/3, 1/3), (1/3, 2/3), (0, 2/3), (0, 1/3), (1/3, 1/3)]:
        print(np.argmax([f(x, y) for f in p3_basis]),
              np.max([f(x, y) for f in p3_basis]))
    return p3_basis


def bernstein_basis(d=3):
    # https://www.iue.tuwien.ac.at/phd/heitzinger/node17.html
    def basis3(u, v): return [(1-u-v)**3,  u**3,  v**3,
                              3*((1-u-v)**2)*u,  3*((1-u-v)) *
                              (u**2),  (3*((u**2)*(v))),
                              (3*((u)*(v**2))), 3*((1-u-v))*(v**2), 3*((1-u-v)**2)*v, 6*((1-u-v)*u*v)]

    def basis3_du(u, v): return [-3*(-u - v + 1)**2,
                                 3*u**2,
                                 np.zeros_like(u),
                                 3*u*(2*u + 2*v - 2) + 3*(-u - v + 1)**2,
                                 -3*u**2 + 2*u*(-3*u - 3*v + 3),
                                 6*u*v,
                                 3*v**2,
                                 -3*v**2,
                                 3*v*(2*u + 2*v - 2),
                                 -6*u*v + 6*v*(-u - v + 1)]

    def basis3_dv(u, v): return [-3*(-u - v + 1)**2,
                                 np.zeros_like(v),
                                 3*v**2,
                                 3*u*(2*u + 2*v - 2),
                                 -3*u**2,
                                 3*u**2,
                                 6*u*v,
                                 -3*v**2 + 2*v*(-3*u - 3*v + 3),
                                 3*v*(2*u + 2*v - 2) + 3*(-u - v + 1)**2,
                                 -6*u*v + 6*u*(-u - v + 1)]

    def basis2(x, y): return [(-x - y + 1)**2, x**2,
                              y**2,  2*x*(-x - y + 1), 2*x*y, 2*y*(-x - y + 1)]
    if d == 3:
        return basis3, basis3_du, basis3_dv
    else:
        assert d == 2
        return basis2, None, None


def jacobians_sample(P, basis_d, p3_codec):
    I4 = np.eye(4)
    tetra20_codec = '000,111,222,333,001,011,112,122,022,002,033,003,233,223,133,113,012,013,023,123'
    gmsh_codec = np.vstack([I4[int(i[0])]+I4[int(i[1])]+I4[int(i[2])]
                            for i in tetra20_codec.split(',')]).astype(np.int)

    def invert_permutation(p):
        s = np.empty_like(p)
        s[p] = np.arange(p.size)
        return s

    p3_nodes = (P.reshape(-1, 20, 3))[:, np.lexsort(gmsh_codec.T),
                                      :][:, invert_permutation(np.lexsort(p3_codec.T)), :]

    dxyz = np.array([np.asarray([[basis_d[j](*i) for i in pts]
                                 for j in range(3)]) for pts in p3_nodes])

    jacs = np.linalg.det(np.einsum('Nujk,Nkv->Njuv', dxyz, p3_nodes))
    return jacs


def bunny_dataquery(B, T, F, order=2, vis_lv=4):
    import seism
    import tqdm
    basis, basis_du, basis_dv = bernstein_basis(order)
    edge_info = edge_globals(F)
    newV, newF, newE, SVI, SVJ, values_at_points, all_base, all_dir, all_global_ids = upsampled(
        B, T, F, basis, edge_info, level=vis_lv)
    pc = seism.PrismCage('../buildr/bunny.offfit.h5')
    aabb = seism.AABB(pc.refV, pc.refF)
    queries = -np.ones((len(newV), 3))
    for i, (v, n) in enumerate(zip(tqdm.tqdm(newV), all_dir.reshape(-1, 3)[SVI])):
        side1 = aabb.segment_hit(v, v+n, True)
        queries[i] = side1[:-1]
    qu, qv = queries[:, 1:].T
    query_data = np.einsum(
        'qed,qe->qd', pc.refV[pc.refF[queries[:, 0].astype(np.int)]], np.vstack([1-qu-qv, qu, qv]).T)
    return query_data


def tetmesh_from_shell(base, top, F):
    '''reminder: in general, this require F to be after roll'''
    tetra_splits = (np.array([0, 3, 4, 5, 1, 4, 2, 0, 2, 5, 0, 4]).reshape(-1, 4),
                    np.array([0, 3, 4, 5, 1, 4, 5, 0, 2, 5, 0, 1]).reshape(-1, 4))
    vnum = len(base)
    T = []
    for f in F:
        tet_c = tetra_splits[0] if f[1] > f[2] else tetra_splits[1]
        T.append((tet_c // 3)*vnum + f[tet_c % 3])
    return np.vstack([base, top]), np.vstack(T)

import numba
@numba.jit(nopython=True)
def rowwise_norm(V):
    norms = np.zeros(V.shape[0])
    for i in range(V.shape[0]):
        norms[i] = np.linalg.norm(V[i])
    return norms


@numba.jit(nopython=True)
def mips_energy(V, f):
    e1 = V[f[:, 1]] - V[f[:, 0]]
    e2 = V[f[:, 2]] - V[f[:, 0]]
    e1_len = rowwise_norm(e1)
    e2_x = np.sum(e1*e2, axis=1) / e1_len
    e2_y = rowwise_norm(e2 - (e2_x / e1_len).reshape(-1, 1) * e1)
    a, b, c = e1_len, -e1_len / \
        np.sqrt(3)+e2_x*2/np.sqrt(3), 2*e2_y/np.sqrt(3)  # Jacobian = (a,b;0,c)
    jac_det = a*c
    jac_frob = a**2 + b**2 + c**2
    return jac_frob/jac_det


@torch.jit.script
def mips_energy_pt(V, f):
    e1 = V[f[:, 1]] - V[f[:, 0]]
    e2 = V[f[:, 2]] - V[f[:, 0]]
    e1_len = torch.norm(e1, dim=1, p=2)
    e2_x = torch.sum(e1*e2, dim=1) / e1_len
    e2_y = torch.norm(e2 - (e2_x / e1_len)[:, None] * e1, dim=1, p=2)
    a, b, c = e1_len, -e1_len / \
        math.sqrt(3)+e2_x*2/math.sqrt(3), 2*e2_y / \
        math.sqrt(3)  # Jacobian = (a,b;0,c)
    jac_det = a*c
    jac_frob = a**2 + b**2 + c**2
    return jac_frob/jac_det


def sp2pt_sparse(along_mat):
    return torch.sparse_coo_tensor(torch.LongTensor(np.vstack([along_mat.row, along_mat.col])), torch.DoubleTensor(along_mat.data))


def upsampled(sB, sT, F, basis, edge_global_info, level=5):
    base, direction = sB, sT-sB
    usV, usF = igl.upsample(
        np.array([[0, 0.], [1, 0], [0, 1]]), np.array([[0, 1, 2]]), level)
    bnd0 = igl.boundary_loop(usF)
    usE = np.vstack([bnd0, np.roll(bnd0, -1)]).T
    u, v = usV[:, :1], usV[:, 1:]
    values_at_nodes = np.hstack(basis(u, v)).T
    if values_at_nodes.shape[0] == 6:
        global_id_for_face = global_id_for_face_quadratic
    if values_at_nodes.shape[0] == 10:
        global_id_for_face = global_id_for_face_cubic

    all_global_ids = np.vstack(
        [global_id_for_face(i, len(sB), F, edge_global_info).astype(int) for i in range(len(F))])

    all_base = np.einsum('fjk,bj->fbk', base[F], np.hstack([1-u-v, u, v]))
    all_dir = np.einsum('fjk,bj->fbk', direction[F], np.hstack([1-u-v, u, v]))

    F_combine = np.vstack([usF+len(u)*i for i in range(len(F))])
    E_combine = np.vstack([usE+len(u)*i for i in range(len(F))])

    newV, SVI, SVJ, _ = igl.remove_duplicate_vertices(
        all_base.reshape(-1, 3), F_combine, 1e-10)
    newF, newE = SVJ[F_combine], SVJ[E_combine]
    # assert len(igl.boundary_loop(newF)) == 0
    return newV, newF, newE, SVI, SVJ, values_at_nodes, all_base, all_dir, all_global_ids

def curved_amips(cp, deri_u, deri_v):
    DU, DV = np.einsum('fed,se->fsd', cp, deri_u), np.einsum('fed,se->fsd', cp, deri_v)

    e1_len = np.linalg.norm(DU,axis=2)
    e2_x = (DU*DV).sum(axis=2)/e1_len
    e2_y = np.linalg.norm(DV - DU * (e2_x/e1_len)[:,:,None],axis=2)

    inv_ref = np.array([[1,-1/np.sqrt(3)], [0, 2/np.sqrt(3)]])

    quality = ((e1_len**2 + (e2_x - e2_y/np.sqrt(3))**2 + (e2_y**2)*(4/3))/
               (e1_len * 2*e2_y/np.sqrt(3)))
