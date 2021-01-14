#!/usr/bin/env python
import pdb
import numpy as np
from scipy.optimize.optimize import fmin_powell
import torch as th
import igl
import scipy
import tqdm
import scipy.optimize
from sksparse.cholmod import cholesky
if __package__:
    from .fem_generator import codecs
    from .utils import *
else:
    from utils import *
    from fem_generator import codecs
import quadpy
import meshio
import pathlib
import time
datapath = pathlib.Path(__file__).parent.absolute() / 'data'

device = 'cpu'
def mips_p4(nodes35):
    '''
    nodes35: (num_tets, 35, 3)
    '''
    if mips_p4.first_run:
        with np.load(datapath/f'p4_q{mips_p4.rule}_dxyz.npz') as npl:
            mips_p4.dxyz, mips_p4.weights, pts = map(
                lambda x: th.from_numpy(npl[x]).to(device), ['dxyz', 'weights', 'points'])
        mips_p4.first_run = False
    jacs = (mips_p4.dxyz@(nodes35.unsqueeze(1))
            ).transpose(1, 2).reshape(-1, 3, 3)
    dets = th.sum(jacs[:, :, 0] *
                  th.cross(jacs[:, :, 1], jacs[:, :, 2]), dim=-1)
    if dets.min() <= 0:
        return th.tensor([1e10]).to(device)
    frob2 = th.sum(jacs.reshape(-1, 9)**2, dim=-1)
    mipses = (frob2/dets**(2/3)).reshape(len(nodes35), -1)
    return (mipses*mips_p4.weights).sum(axis=1)


mips_p4.first_run = True
mips_p4.rule = 'wsj'


def mips_p1(nodes4):
    """MIPS energy of P1 elements

    Args:
        nodes4 (num_tets, 4, 3): all positions

    Returns:
        array: num_tets double
    """
    jacs = nodes4[:, 1:] - nodes4[:, :1]
    dets = th.sum(
        jacs[:, :, 0]*th.cross(jacs[:, :, 1], jacs[:, :, 2], dim=1), dim=-1)
    frob2 = th.sum(jacs.reshape(-1, 9)**2, dim=-1)
    mipses = (frob2/dets**(2/3))
    return mipses


def vertex_tetra_adjacency(n: int, pT, dim=3):
    """Adjacency from vertex to (possibly high order) tetras
    Warning: this function assumes vertex ordered first

    Args:
        n (int): number of vertices
        pT (ndarray, T by node_total):
        dim (int, optional): Defaults to 3.

    Returns:
        list of list: tetras around each vertex
    """
    VT = [[] for _ in range(n)]
    for i, t in enumerate(pT):
        for tt in t[:dim+1]:  # this depends on the vertex-first ordering
            VT[tt].append(i)
    return VT


def vertex_node_adjacency(n: int, pT, codec_n):
    """for P4, the affected nodes of a vertex, excluding the link, which is supposed to be fixed during optimization

    Args:
        n (int): number of vertices (not nodes).
        pT (ndarray): Face/Cell Array 
        codec_n ([type]): integral codec for the ordering (contribution based)

    Returns:
        list of list: nodes around each vertex one ring (exclude vertices)
    """
    assert pT.shape[1] == codec_n.shape[0]
    num = codec_n.shape[1]
    local_affected = [np.where([i in j for j in codec_n])[0]
                      for i in range(num)]
    VN = [set([i]) for i in range(n)]
    for t in pT:
        for i, v in enumerate(t[:num]):
            VN[v] |= set(map(int, t[local_affected[i]]))
    assert pT[0, 1] not in VN[pT[0, 0]], 'vertices should not be contained'
    return VN


def energy_for_one_ring(v: int, VT, n_pos, complete_tuples):
    energy = 0

    # p4
    p4_tets = th.tensor(VT[v])
    localnodes = n_pos[complete_tuples[p4_tets]].mean(dim=2)
    energy += mips_p4(localnodes).sum()

    if np.isnan(energy.item()):  # flip detection at these samples.
        return 1e10
    return energy


def test_dxyz():
    # nodes during 3D computation are all Lagrange now.
    rule = 'uniform5'
    with np.load(datapath/f'p4_q{rule}_dxyz.npz') as npl:
        dxyz, weights, points = npl['dxyz'], npl['weights'], npl['points']
    codec = codecs()['tetra35']
    tetra35_op = codec[2]
    x, y, z = points[:, 1], points[:, 2], points[:, 3]

    # with trivial tet, all the detereminants should be 1.
    dets = np.linalg.det(np.einsum('dqL,LD->qDd', dxyz, tetra35_op[:, 1:]))
    assert np.linalg.norm(dets - 1) < 1e-10

    space_warp = (tetra35_op[:, 1:]+1)**2
    dets = np.linalg.det(np.einsum('dqL,LD->qDd', dxyz, space_warp))
    assert np.linalg.norm(dets - 8*(x+1)*(y+1)*(z+1)) < 1e-10

    space_warp = np.vstack(
        [tetra35_op[:, 1]**2, tetra35_op[:, 2]**3, tetra35_op[:, 3]**4]).T
    dets = np.linalg.det(np.einsum('dqL,LD->qDd', dxyz, space_warp))
    assert np.linalg.norm(dets - 2*x*3*y**2*4*z**3) < 1e-10


def test_p4_energy_ad():
    with np.load('mixed_mesh.npz') as npl:
        # lag_pos includes internal linear and dependent nodes, and the dependecy relationship is revealed by dp_trans
        p4T, p1T, lag_pos, dp_trans, sf_nodes, uniq_inv = map(
            lambda x: th.from_numpy(npl[x]), ('p4T', 'p1T', 'pos', 'dp_trans', 'sf_nodes', 'uniq_inv'))
    complete_tuples = dp_trans[th.take(
        uniq_inv, (35*th.arange(len(p4T))[:, None]+th.arange(35)))]

    VT = [[24, 25, 26, 555, 556, 557, 702, 703, 704,
           711, 712, 713, 966, 967, 968, 990, 991, 992]]
    lag_pos.requires_grad_()
    e = energy_for_one_ring(0, VT, lag_pos, complete_tuples)
    e.backward()
    assert abs(e.item() - 110.56881646233114 -
               99.24246388164019) < 1e-10, 'empirical, with respect the right tetrahedron, may require transformation later'


def find_nodes_on_surface(p4T, cod):
    bf, _ = igl.tet_tet_adjacency(p4T[:, :4])
    tt_spec = np.array([[0, 1, 2], [0, 1, 3], [1, 2, 3], [2, 0, 3]])
    tets, faces = np.where(bf == -1)
    marker = np.ones_like(p4T[:, :4])
    for t, f in zip(tets, tt_spec[faces]):
        marker[t, f] = 0
    mask = (marker[:, cod].sum(axis=2) < 4)
    return np.unique(p4T[mask])


def test_optim(out_file='block_opt.msh'):
    import meshio
    mesh = meshio.read('/home/zhongshi/public/curved/block.msh')
    lag_pos = mesh.points
    p4T = mesh.cells[0].data
    cod = codecs()['tetra35'][1]  # gmsh codec
    sf_nodes = find_nodes_on_surface(p4T, cod)

    # vertices excluding sf. i.e. dof to optimize
    vertices = list(set(p4T[:, :4].flatten()) - set(sf_nodes))

    p4T, lag_pos = map(th.from_numpy, (p4T, lag_pos))
    p4T = p4T.to(device)
    lag_pos = lag_pos.to(device)
    mips_p4.rule = 'uniform20'

    print('number of vertices', len(vertices))
    codec_35n = codecs()['tetra35'][1]
    VN = [th.tensor(list(v))
          for v in vertex_node_adjacency(len(lag_pos), p4T, codec_35n)]
    VT = vertex_tetra_adjacency(len(lag_pos), p4T)
    for _ in range(20):
        for v in tqdm.tqdm(vertices):
            # clone is important, otherwise aliasing error
            free_vars = lag_pos[VN[v]].clone().requires_grad_()
            print(f'free {free_vars.shape}')
            optimizer = th.optim.LBFGS(
                [free_vars], lr=1e-2, line_search_fn='strong_wolfe')
            p4indices = p4T[VT[v]].long()

            def closure(free_vars, grad=True):
                lag_pos[VN[v]] = free_vars
                e = mips_p4(lag_pos[p4indices]).sum()
                if np.isnan(e.item()) or e > 1e9:
                    return 1e10  # infinity, only activated during line search
                if grad:
                    e.backward(retain_graph=True)
                return e
            print(f'vertex {v} {closure(free_vars)}')
            for _ in range(10):
                optimizer.zero_grad()
                e = closure(free_vars)
                ini_val = e.item()
                g = free_vars.grad.data
                alpha = 1e-4
                with torch.no_grad():
                    for _ in range(20):
                        cur_val = closure(free_vars - alpha*g, False)
                        if ini_val > cur_val:
                            print(cur_val.item(),end=' ')
                            break
                        alpha = alpha /2
                    else:
                        print('LS fail')
                        break
                free_vars.data -= alpha*g.data
            lag_pos = lag_pos.detach()
            print()
        meshio.write(out_file,meshio.Mesh(points=lag_pos.cpu().numpy(), cells=[('tetra35', p4T.cpu().numpy())]))


def write_lags(file='lag_pos.npy', rule='uniform50'):
    with np.load('mixed_mesh.npz') as npl:
        p4T, p1T, lag0, dp_trans, sf_nodes, uniq_inv = map(
            lambda x: th.from_numpy(npl[x]), ('p4T', 'p1T', 'pos', 'dp_trans', 'sf_nodes', 'uniq_inv'))
    lag1 = th.from_numpy(np.load(file))
    mips_p4.rule = rule
    p4tets = th.arange(len(p4T))
    complete_tuples = dp_trans[th.take(uniq_inv,
                                       (35*p4tets[:, None]+th.arange(35)))]
    print(complete_tuples.shape, complete_tuples.sum())
    _, idp_nodes = np.unique(
        dp_trans[:, 0], return_index=True)  # indeprendent nodes
    dep_nodes = np.asarray(
        list(set(range(len(dp_trans))) - set(idp_nodes)))  # dependent nodes
    lag1[dep_nodes] = lag1[dp_trans[dep_nodes]].mean(axis=-2)
    for n_pos in [lag0, lag1]:
        localnodes = n_pos[complete_tuples].mean(axis=-2)
        print(mips_p4(localnodes).sum()/len(p4tets))
    meshio.write('temp1.msh',
                 meshio.Mesh(lag1.numpy(), [('tetra35', p4T.numpy())]))


if __name__ == '__main__':
    import fire
    fire.Fire()
