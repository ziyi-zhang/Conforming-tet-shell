#!/usr/bin/env python
# coding: utf-8

# In[2]:


## import glob
import os
import plotly.graph_objects as go
import tqdm
from sidecar import Sidecar
from vis_utils import *


# In[3]:


import sys
sys.path.append('../python')
import seism
from curve import assembler, cubic_fitting
from curve import fem_generator


# In[4]:


basis = fem_generator.basis_info(nsd=2,order=3,force_codec='tri10')['basis']
usV, usF = igl.upsample(
    np.array([[0, 0.], [1, 0], [0, 1]]), np.array([[0, 1, 2]]), 5)
bnd0 = igl.boundary_loop(usF)
usE = np.vstack([bnd0, np.roll(bnd0, -1)]).T
u, v = usV[:, :1], usV[:, 1:]
sampled_basis_value = np.hstack(basis(u, v)).T


# In[13]:


refV,refF,B,T,V,F,cp = h5reader('../build_clang/spot_collapse_n95.h5', 'ref.V', 'ref.F', 
                 'mbase', 'mtop', 'mV', 'mF', 'complete_cp')


# In[14]:


F.shape


# In[15]:


sv = np.einsum('fcd,cs->fsd', 
               cp, sampled_basis_value)
all_F = np.vstack([usF+len(usV)*i for i in range(len(sv))])
newV,SVI,SVJ,_ = igl.remove_duplicate_vertices(sv.reshape(-1,3), all_F, 1e-6)
newF = SVJ[all_F]
plt = mplot(newV,newF,wireframe=False,flat=False, width=1000,height=1000)
#plt.add_edges(sv.reshape(-1,3),np.vstack([usE+len(usV)*i for i in range(len(sv))]))


# In[10]:


from curve import assembler, cubic_fitting

newV, newF, newE, SVI, SVJ, values_at_points, all_base, all_dir, all_global_ids = cubic_fitting.upsampled(B, T, F, basis,  cubic_fitting.edge_globals(F), level=4)


# In[11]:


M_p = seism.AABB(refV,refF)
coord, hits =[], []
for b, t in zip(tqdm.tqdm(all_base.reshape(-1,3)),
                (all_base+all_dir).reshape(-1,3)):
    hit = M_p.segment_hit(b[None,:], t[None,:], True)
    assert(hit[0]>=0)
    hits.append(hit)
    coord.append(b + hit[-1]*(t-b))
h_fid, h_u, h_v, h_t = np.array(hits).T


# In[12]:


operator = assembler.spline_operator(
        values_at_points, all_global_ids, SVI, SVJ, all_dir, all_global_ids.max()+1, len(F))


# In[13]:


v_data = np.einsum('fed,fe->fd',
                   refV[refF[h_fid.astype(np.int)]],
                   np.array([1-h_u-h_v, h_u, h_v]).T).reshape(len(F),-1,3)
hitV = v_data.reshape(-1,3)[SVI]
massmat = igl.massmatrix(hitV, newF).data
diagmat = scipy.sparse.diags(massmat)
sol0 = scipy.sparse.linalg.spsolve(
    operator.T@diagmat@operator, operator.T@diagmat@(hitV))


# In[ ]:


def mean_curv_normal(V,F):
    massmat = igl.massmatrix(V, F, igl.MASSMATRIX_TYPE_VORONOI).data
    HN = -scipy.sparse.diags(1/massmat)@(igl.cotmatrix(V,F)@V)
    return HN


# In[14]:


plt = mplot(operator@sol0,newF,
                #c = np.linalg.norm(operator@sol0 - v_data.reshape(-1,3)[SVI],axis=1),
                wireframe=False,width=1000,height=1000)


# In[ ]:




