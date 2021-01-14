#!/usr/bin/env python
# coding: utf-8

# In[1]:


from vis_utils import *
import sys
sys.path.append('../python')


# In[2]:


from curve import fem_generator
import meshio


# In[4]:


tetra_o4 = fem_generator.basis_info(order=4, nsd=3, derivative=True)


# In[7]:


with np.load('../python/curve/data/p4_quniform5_dxyz.npz') as fp:
    dxyz, pts = fp['dxyz'], fp['points']


# In[8]:


deris5 = fem_generator.bernstein_deri_evaluator(*pts[:,1:].T, 
                                       codecs=np.asarray(tetra_o4['codec']))


# In[30]:


gmsh_cod = (fem_generator.codecs()['tetra35'][2]*4).astype(np.int)
reorder = np.lexsort(np.asarray(tetra_o4['codec']).T)[
    fem_generator.invert_permutation(np.lexsort(gmsh_cod.T))
]


# In[16]:


deris = fem_generator.bernstein_deri_evaluator(*np.asarray(tetra_o9['codec'])[:,1:].T/9, 
                                               codecs=np.asarray(tetra_o4['codec']))


# In[4]:


tetra_o9 = fem_generator.basis_info(order=9, nsd=3, derivative=False)


# In[24]:


lagr = np.linalg.det(deris.transpose((2,0,1))@cp)


# In[35]:


import meshio
cp = np.asarray(
    tetra_o4['codec'])[:,1:] - 1e-1*np.sin(np.arange(35*3).reshape(-1,3))
meshio.write('test.msh',
             meshio.Mesh(points=tetra_o4['b2l']@cp, 
                         cells=[('tetra35', reorder[None,:])])
            )


# In[6]:


with h5py.File('../python/curve/data/tetra_o9_l2b.h5','w') as fp:
    fp['b2l'] = tetra_o9['b2l']
    fp['l2b'] = tetra_o9['l2b']


# In[8]:


tetra_o9['codec'][50]

