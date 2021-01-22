if __package__:
    from . import utils, fem_generator
else:
    import utils
    import fem_generator
import numpy as np
import meshio


def conforming_tetegen(mB, mF):
    import pymesh
    tetgen = pymesh.tetgen()
    tetgen.points = mB
    tetgen.triangles = mF
    tetgen.verbosity = 3
    tetgen.coarsening = False
    tetgen.split_boundary = False
    tetgen.run()
    assert np.linalg.norm(
        tetgen.mesh.vertices[:len(mB)]-mB) < 1e-10, "preserving order"
    return tetgen


def surface_to_curved_layer(B, F, cp, bern2elevlag, tetgen=None):
    '''
    bern2elevlag comes from `cumin_datagen.py`
    '''
    tri_o4 = fem_generator.basis_info(order=4, nsd=2)
    def codec_to_n(co): return [k for i, j in enumerate(co) for k in [i]*j]
    tri_o4_n = np.array([codec_to_n(c) for c in tri_o4['codec']])

    val_at_samples = np.einsum('fed,eE->fEd', cp, bern2elevlag)
    _, Tbt = utils.tetmesh_from_shell(B, B, F)

    tet_o4_n = fem_generator.codecs()['tetra35'][1]  # gmsh codec

    elev_base = B[F[:, tri_o4_n]].mean(axis=2).reshape(-1, 3)
    elev_dir = val_at_samples.reshape(-1, 3) - elev_base

    nodes_height = (Tbt[:, tet_o4_n]//len(B)).mean(axis=2)
    triple2index = {tuple(sorted(f)): i for i, f in enumerate(
        F[:, tri_o4_n].reshape(-1, tri_o4_n.shape[1]))}
    tetnode_index = [triple2index[tuple(sorted(t))] for t in (
        Tbt[:, tet_o4_n] % len(B)).reshape(-1, tet_o4_n.shape[1])]

    tetnodebase = elev_base[tetnode_index]
    tetnodedir = elev_dir[tetnode_index] * (nodes_height.reshape(-1, 1))
    combined_tuple = np.sort(Tbt[:, tet_o4_n].reshape(-1, 4))
    combined_verts = tetnodebase + tetnodedir
    if tetgen is not None:
        (tetV, tetT) = tetgen
        tg_vox = tetT
        tg_nodes = tetV[tg_vox[:, tet_o4_n]].mean(axis=2).reshape(-1, 3)
        tg_vox = np.where(tg_vox < len(B), tg_vox, tg_vox + len(B))
        tg_tuples = np.sort(tg_vox[:, tet_o4_n].reshape(-1, 4))
        combined_tuple = np.vstack([combined_tuple, tg_tuples])
        combined_verts = np.vstack([combined_verts, tg_nodes])
    uniq_tuples, uniq_ind, uniq_inv = np.unique(combined_tuple,
                                                axis=0, return_index=True, return_inverse=True)

    assert np.all(combined_tuple[uniq_ind] == uniq_tuples)
    num_tets = len(combined_tuple)//35
    return combined_verts[uniq_ind], uniq_inv[np.tile(np.arange(35), num_tets).reshape(num_tets, -1) + np.arange(num_tets)[:, None]*35]
    #mesh = meshio.Mesh(combined_verts[uniq_ind],
    #                   [("tetra35", uniq_inv[np.tile(np.arange(35), num_tets).reshape(num_tets, -1) + np.arange(num_tets)[:, None]*35])])
    
    # if tetgen is not None:
        # assert len(np.unique(mesh.cells[0].data[:,:4])) == len(tetgen[0]) + len(B)
    #return mesh

def sol_to_curved_layer(B, T, F, tV, tT, tF, sol):
    _, Tbt = utils.tetmesh_from_shell(B, T, F)

    edge_info = utils.edge_globals(F)
    basis, _, _ = utils.bernstein_basis(2)
    all_global_ids = np.vstack([utils.global_id_for_face_quadratic(
        i, len(B), F, edge_info) for i in range(len(F))])

    codec = utils.codecs()
    _, tri15_codec_n, tri15_op = codec['tri15']
    _, tetra35_codec_n, _ = codec['tetra35']

    val_at_samples = np.vstack(basis(tri15_op[:, 1], tri15_op[:, 2]))
    sampled_heights = ((np.clip(sol, 0, 1)
                        [all_global_ids])@val_at_samples).flatten()

    elev_base = B[F[:, tri15_codec_n]].mean(axis=2).reshape(-1, 3)
    elev_dir = (T-B)[F[:, tri15_codec_n]].mean(axis=2).reshape(-1, 3)
    nodes_height = (Tbt[:, tetra35_codec_n]//len(B)).mean(axis=2)
    triple2index = {tuple(sorted(f)): i for i, f in enumerate(
        F[:, tri15_codec_n].reshape(-1, tri15_codec_n.shape[1]))}
    tetnode_index = [triple2index[tuple(sorted(t))] for t in (
        Tbt[:, tetra35_codec_n] % len(B)).reshape(-1, tetra35_codec_n.shape[1])]

    tetnodebase = elev_base[tetnode_index]
    tetnodedir = elev_dir[tetnode_index] * \
        (sampled_heights[tetnode_index, None]*nodes_height.reshape(-1, 1))
    all_tuples = np.sort(Tbt[:, tetra35_codec_n].reshape(-1, 4))
    # meshio.write('test.msh',
    #              meshio.Mesh((tetnodebase + tetnodedir),
    #                          [("tetra35",np.arange(4719*35, 4720*35)[None,:])]
    #                          ))
    # return
    uniq_tuples, uniq_ind, uniq_inv = np.unique(
        all_tuples, axis=0, return_index=True, return_inverse=True)
    assert np.all(all_tuples[uniq_ind] == uniq_tuples)

    meshio.write('test.msh',
                 meshio.Mesh((tetnodebase + tetnodedir)[uniq_ind],
                             [("tetra35",
                               uniq_inv[np.tile(np.arange(35), len(Tbt)).reshape(
                                   len(Tbt), -1) + np.arange(len(Tbt))[:, None]*35]  # [4719: 4720]
                               )]
                             ))

def test_connected(p4T):
    TT,TTi = igl.tet_tet_adjacency(p4T[:,:4])
    for i,t in enumerate(TT):
        for j in t:
            if j != -1:
                assert(len(np.intersect1d(p4T[i], p4T[j]))>=15)

def reorder_tetra(m):
    gmsh_cod = (fem_generator.codecs()['tetra35'][-1]*4).astype(np.int)

    auto_cod = fem_generator.tuple_gen(order=4, var_n=3)

    reorder = np.lexsort(
            np.array(auto_cod).T)[fem_generator.invert_permutation(np.lexsort(gmsh_cod.T))]

    assert np.all(np.array(auto_cod)== gmsh_cod[fem_generator.invert_permutation(reorder)])

    return m.points, m.cells[0].data[:,fem_generator.invert_permutation(reorder)]


if __name__ == '__main__':
    import sys
    import h5py
    import igl
    import meshio
    import os

    assert(len(sys.argv) > 1)
    input_filename = sys.argv[1]
    curveFolder = './'
    if (len(sys.argv) > 2):
        curveFolder = sys.argv[2]
    interior = True

    str_ = '/home/ziyi/TetShell/data/0112_exterior/' + os.path.basename(input_filename)
    output_filename_msh = os.path.splitext(str_)[0] + 'stitch.msh'
    output_filename_h5  = os.path.splitext(str_)[0] + 'stitch.h5'
    curve_filename = os.path.basename(input_filename)
    curve_filename = os.path.join(curveFolder, os.path.splitext(curve_filename)[0])
    curve_filename = curve_filename[:-1] + '5'  # remove the annoying '_'

    # read params
    with h5py.File("/home/ziyi/TetShell/code/Conforming-tet-shell/python/curve/data/tri_o3_lv3.h5",'r') as fp:
        bern2elevlag = fp['bern2elevlag'][()]
    # read tet mesh
    print('input mesh: {}'.format(input_filename))
    mesh = meshio.read(input_filename)
    V_msh = mesh.points
    T_msh = mesh.cells[0][1]
    label = mesh.cell_data['label'][0]

    # read curved mesh
    print('curved mesh: {}'.format(curve_filename))
    with h5py.File(curve_filename, 'r') as fp:
        cp = fp['complete_cp'][()]
        mT = fp['mtop'][()]
        mB = fp['mbase'][()]
        mV = fp['mV'][()]
        mF = fp['mF'][()]

    #Vb, Tb = surface_to_curved_layer(mB, mF, cp, bern2elevlag, (V_msh, T_msh[label == 1]))
    #Vt, Tt = surface_to_curved_layer(mT, mF, cp, bern2elevlag, (V_msh, T_msh[label == 2]))
    #Tt += len(Vb)
    #m = meshio.Mesh(np.vstack((Vb, Vt)), [("tetra35", np.vstack((Tb, Tt)))])  # with duplicated vertices
    if interior:
        V, T = surface_to_curved_layer(mB, mF, cp, bern2elevlag, (V_msh, T_msh[label == 1]))
        print('interior={}   T_msh size={}'.format(interior, T_msh[label == 1].shape[0]))
    else:
        N = mV.shape[0]
        Vb = V_msh[0:N, :]
        Vt = V_msh[N:2*N, :]
        V_msh[0:N, :] = Vt
        V_msh[N:2*N, :] = Vb
        mask_b = T_msh < N
        mask_t = np.logical_and(T_msh>=N, T_msh<2*N)
        T_msh[mask_b] += N
        T_msh[mask_t] -= N

        #V_msh = V_msh[mV.shape[0]:, :]
        #T_msh -= mV.shape[0]
        V, T = surface_to_curved_layer(mV, mF, cp, bern2elevlag, (V_msh, T_msh[label == 2]))
        print('interior={}   T_msh size={}'.format(interior, T_msh[label == 2].shape[0]))
    m = meshio.Mesh(V, [("tetra35", T)])

    # write as msh to output_filename_msh
    print('result written to: {}'.format(output_filename_msh))
    meshio.write(output_filename_msh, m)

    # write as h5 to output_filename_h5
    print('result written to: {}'.format(output_filename_h5))
    hP, hC = reorder_tetra(m)
    with h5py.File(output_filename_h5, 'w') as fp:
        fp['lagr'] = hP
        fp['cells'] = hC
