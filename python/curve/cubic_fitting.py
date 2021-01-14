import numba
import math
import time
import osqp
import torch
import numpy as np
import igl
import scipy
import tqdm
import scipy.optimize
from sksparse.cholmod import cholesky
import pdb
if __package__:
    from .utils import *
    from .assembler import *
else:
    from utils import *
    from assembler import *


import pathlib
datapath = pathlib.Path(__file__).parent.absolute() / 'data'


def sphere_data_prepare(sB, sT, F):
    def modify_sqrt(V, F):
        V = np.vstack([V, V[F[-1]].mean(axis=0)])
        F = np.vstack([F[:-1], [[4, 0, 8], [0, 2, 8], [2, 4, 8]]])
        return V, F

    def modify_sub(V, F):
        V = np.vstack([V, V[1:3].mean(axis=0)])
        F = np.vstack([[[0, 1, 8], [0, 8, 2], [1, 3, 8], [3, 2, 8]], F[2:]])
        return V, F

    def simple_spherical_intersect(sB, sT):
        sD = sT-sB
        dd = np.sum(sD*sD, axis=1)
        db = 2*np.sum(sD*sB, axis=1)
        bb = np.sum(sB*sB, axis=1) - 1
        ray_steps = (-db + np.sqrt(db**2 - 4*dd*bb))/(2*dd)
        # assert np.linalg.norm(sB + ray_steps[:,None]*(sD),axis=1) all one
        return ray_steps

    # data generation
    triangle_uv = []
    #idptr = [0]
    for u in np.linspace(0, 1, num=20):
        for v in np.linspace(0, 1-u, num=20):
            triangle_uv.append((u, v))
    triangle_uv = np.array(triangle_uv)
    # idptr.append(len(all_uv))

    tri_bc = np.hstack(
        [(1-triangle_uv.sum(axis=1, keepdims=True)), triangle_uv])

    idptr = np.arange(len(F)+1)*len(triangle_uv)
    uvtn = np.zeros(((len(F))*len(tri_bc), 6))

    for i, f in enumerate(F):
        steps = simple_spherical_intersect(tri_bc@sB[f], tri_bc@sT[f])
        uvtn[idptr[i]:idptr[i+1], :2] = triangle_uv
        uvtn[idptr[i]:idptr[i+1], 2] = steps
        uvtn[idptr[i]:idptr[i+1], 3:] = tri_bc@sB[f] * \
            (1-steps[:, None]) + tri_bc@sT[f]*(steps[:, None])
    return uvtn


def approximate_continuity():
    V, F = igl.read_triangle_mesh(
        '/home/zhongshi/Workspace/libigl/tutorial/data/cube.obj')
    V, F = modify_sqrt(V, F)
    level = 3
    V, F = igl.upsample(V, F, level)
    basis, basis_du, basis_dv = bernstein_basis()
    sB = V/np.sqrt(3)
    sT = V*np.sqrt(3)*2
    edge_info = edge_globals(F)
    TT, TTi = igl.triangle_triangle_adjacency(F)
    along_mat, ortho_mat, oppo_mat, pos_mat = perpendicular_direction_op(
        V, F, TT, TTi, edge_info, basis, basis_du, basis_dv)
    newV, newF, newE, SVI, SVJ, values_at_points, all_base, all_dir, all_global_ids = upsampled(
        sB, sT, F, basis, edge_globals(F), level=8-level)

    massmat = igl.massmatrix(newV, newF, igl.MASSMATRIX_TYPE_VORONOI).data
    diagmat = scipy.sparse.diags(massmat)
    sp_dh_D_s = spline_operator(
        values_at_points, all_global_ids, SVI, SVJ, all_dir, all_global_ids.max()+1, len(F))
    all_v = newV/(np.linalg.norm(newV, axis=1)[:, None])
    sol0 = []
    for j in range(3):
        prob = osqp.OSQP()
        prob.setup(sp_dh_D_s.T@diagmat@sp_dh_D_s,
                   -sp_dh_D_s.T@diagmat@all_v[:, j].flatten(),
                   verbose=False)
        prob_res = prob.solve()
        sol0.append(prob_res.x)
    sol0 = np.vstack(sol0).T

    along_pt = sp2pt_sparse(along_mat)
    ortho_pt = sp2pt_sparse(ortho_mat)
    oppo_pt = sp2pt_sparse(oppo_mat)
    sol_tensor = torch.DoubleTensor(sol0)
    sol_tensor.requires_grad_()
    optimizer = torch.optim.LBFGS([sol_tensor])

    def closure():
        optimizer.zero_grad()
        loss = ((torch.sum(torch.mul(along_pt@sol_tensor, torch.cross(ortho_pt @
                                                                      sol_tensor, oppo_pt@sol_tensor)), dim=1)**2).sum())
        loss.backward()
        return loss
    for _ in range(100):
        print(optimizer.step(closure))
    np.savez(f's3c_lv{level}_v8_actn.npz', newF=newF,
             opearator=sp_dh_D_s, sol=sol_tensor.detach().numpy())
    import ipdb
    ipdb.set_trace()


def heightfield_bilaplacian():
    V, F = igl.read_triangle_mesh(
        '/home/zhongshi/Workspace/libigl/tutorial/data/cube.obj')
    V, F = modify_sqrt(V, F)
    level = 3
    V, F = igl.upsample(V, F, level)
    basis, basis_du, basis_dv = bernstein_basis()
    sB = V/np.sqrt(3)
    sT = V*np.sqrt(3)*2

    edge_info = edge_globals(F)
    newV, newF, newE, SVI, SVJ, values_at_points, all_base, all_dir, all_global_ids = upsampled(
        sB, sT, F, basis, edge_info, level=8-level)

    massmat = igl.massmatrix(newV, newF).data
    diagmat = scipy.sparse.diags(massmat)
    mass3 = scipy.sparse.diags(np.repeat(massmat, 3))
    sp_dh_D_s = offset_from_ctrlpt_operator(
        values_at_points, all_global_ids, SVI, SVJ, all_dir, all_global_ids.max()+1, len(F))
    all_v = newV/(np.linalg.norm(newV, axis=1)[:, None])
    newbase = all_base.reshape(-1, 3)[SVI]
    newdir = all_dir.reshape(-1, 3)[SVI]
    def update_V(s): return newbase + (sp_dh_D_s@s).reshape(-1, 3)
    sol0 = scipy.sparse.linalg.spsolve(
        sp_dh_D_s.T@mass3@sp_dh_D_s, sp_dh_D_s.T@mass3@(all_v.flatten() - newbase.flatten()))

    operator = sp_dh_D_s
    sol1 = sol0

    def smooth_iter(sol1):
        V1 = update_V(sol1)
        laplacian = -igl.cotmatrix(V1, newF).tocoo()
        m = igl.massmatrix(V1, newF)
        invmass3 = scipy.sparse.diags(1/np.repeat(m.data, 3))

        lap3 = scipy.sparse.coo_matrix((np.concatenate([laplacian.data, laplacian.data, laplacian.data]),
                                        (np.concatenate([laplacian.row*3, laplacian.row*3+1, laplacian.row*3+2]),
                                         np.concatenate([laplacian.col*3, laplacian.col*3+1, laplacian.col*3+2]))))

        print('True', ((m.data**(-1))[:, None] *
                       (laplacian@V1)**2).sum(), end=' ')
        assert not np.any(np.isnan(laplacian.data))
        w = 1e-3
        L_w = w*lap3@invmass3@lap3

        A = operator.T@(L_w + scipy.sparse.eye(lap3.shape[0]))@operator
        A = A.tocsc()
        print('sol', sol1.max(), sol1.min())
        newsol1 = []

        lower = -10*np.ones(len(sol1))
        upper = 10*np.ones(len(sol1))
        b = operator.T@(L_w@newbase.flatten() - operator@sol1)
        lower[:3] = sol1[:3]
        upper[:3] = sol1[:3]
        prob = osqp.OSQP()
        prob.setup(A, b,
                   A=scipy.sparse.eye(len(sol0)).tocsc(),
                   l=lower,
                   u=upper,
                   verbose=False)
        prob.warm_start(x=sol1)
        prob_res = prob.solve()
        sol1 = prob_res.x
        assert not np.any(np.isnan(sol1))
        return sol1
    for it in range(100):
        sol1 = smooth_iter(sol1)
        if it % 10 == 0:
            np.savez(f's3c_lv{level}_v8_height.npz', newF=newF, newE=newE, V1=update_V(
                sol1), operator=sp_dh_D_s, V0=update_V(sol0), sol=sol1)

def heightfield_bunny(warm_start=True, vis_lv=5, step=1e-2, w=1e-3, num_iter=60, prefix='1e-1', constrain=False, order=2,
                      model='bunny'):
    """
    bunny, fit height field onto the shell.

    Args:
        w (double): willmore regularization weight (larger means smoother)
        step (double): flow step size
        vis_lv (int): level of subdivsion for the evaluation of smoothness
        warm_start (bool): warm start
    """
    print(locals())
    with np.load(f'curve/data/{model}_data.npz') as npl:
        query_data, sB, sT, F, vis_lv_ = npl['query_data'], npl['B'], npl['T'], npl['F'], npl['vis_lv']
        assert vis_lv == vis_lv_, 'setting has to be consistent with ser'
    basis, basis_du, basis_dv = bernstein_basis(order)

    edge_info = edge_globals(F)
    newV, newF, newE, SVI, SVJ, values_at_points, all_base, all_dir, all_global_ids = upsampled(
        sB, sT, F, basis, edge_info, level=vis_lv)

    TT, TTi = igl.triangle_triangle_adjacency(newF)
    FN = igl.per_face_normals(newV, newF, np.ones(3))
    angles = np.einsum('fd,fed->fe', FN, FN[TT])
    feid = np.array(np.where(angles <= 0.5)).T
    vert_on_feature = np.unique(np.concatenate([newF[feid[:, 0], feid[:, 1]-2].flatten(),
                                                newF[feid[:, 0], feid[:, 1]].flatten()]))
    feature_flag = np.ones((3, len(newV)))
    # feature_flag[:, vert_on_feature] = 0
    feature_flag = feature_flag.flatten()

    massmat = igl.massmatrix(newV, newF).data
    mass3 = scipy.sparse.diags(np.repeat(massmat, 3))
    sp_dh_D_s = offset_from_ctrlpt_operator(
        values_at_points, all_global_ids, SVI, SVJ, all_dir, all_global_ids.max()+1, len(F))
    newbase = all_base.reshape(-1, 3)[SVI]
    newdir = all_dir.reshape(-1, 3)[SVI]
    def update_V(s): return newbase + (sp_dh_D_s@s).reshape(-1, 3)
    sol0 = scipy.sparse.linalg.spsolve(sp_dh_D_s.T@mass3@sp_dh_D_s, 
        sp_dh_D_s.T@mass3@(query_data.flatten() - newbase.flatten()))
    sol_list = [sol0]
    if warm_start:
        # ,newF=newF, newE=newE, V1=update_V(sol1), operator=sp_dh_D_s, V0 = update_V(sol0), sol=sol1)
        with np.load(f'{prefix}{model}_lv{vis_lv}_height.npz') as npl:
            sol_list = list(npl['sol_list'])
    sol1 = sol_list[-1]

    D = sp_dh_D_s
    V0 = update_V(sol0)

    def smooth_iter(sol1):
        V1 = update_V(sol1)
        laplacian = -igl.cotmatrix(V1, newF).tocoo()
        m = igl.massmatrix(V1, newF)
        mass3 = scipy.sparse.diags(feature_flag*np.repeat(m.data, 3))
        invmass3 = scipy.sparse.diags(feature_flag*1/np.repeat(m.data, 3))

        lap3 = scipy.sparse.coo_matrix((np.concatenate([laplacian.data, laplacian.data, laplacian.data]),
                                        (np.concatenate([laplacian.row*3, laplacian.row*3+1, laplacian.row*3+2]),
                                         np.concatenate([laplacian.col*3, laplacian.col*3+1, laplacian.col*3+2]))))

        willmore = ((m.data**(-1))[:, None]*(laplacian@V1)**2).sum()/4
        proximity = (((V1-V0))**2).sum()/2
        print(
            f'Willmore {willmore} Proximity {proximity} Total {willmore + proximity/w} ')
        assert not np.any(np.isnan(laplacian.data))
        # mass -> identity
        mass3 = scipy.sparse.eye(mass3.shape[0])
        L_w = lap3@invmass3@lap3/2
        P = step*D.T@(L_w + mass3/w)@D + scipy.sparse.eye(len(sol0))
        q = step*D.T@(mass3@D@sol0/w - L_w@newbase.flatten()) + sol1
        if not constrain:
            return scipy.sparse.linalg.spsolve(P, q)
        prob = osqp.OSQP()
        prob.setup(P.T@P, -P.T@q,
                   A=scipy.sparse.eye(len(sol0)).tocsc(),
                   l=sol0*0,
                   u=sol0*0+1,
                   verbose=False)
        prob.warm_start(x=sol1)
        prob_res = prob.solve()
        sol1 = prob_res.x
        sol1[sol1 >= 1] = 1.
        sol1[sol1 <= 0] = 0.
        assert not np.any(np.isnan(sol1))
        return sol1
    for it in range(num_iter):
        sol1 = smooth_iter(sol1)
        sol_list.append(sol1)
        np.savez(f'{prefix}{model}_lv{vis_lv}_height.npz',
                 newF=newF, newE=newE, V1=update_V(sol1), operator=sp_dh_D_s, sol_list=sol_list)


if __name__ == '__main__':
    import fire
    fire.Fire()
