# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 13:03:24 2020

@author: Mohammad Asif Zaman

Code for testng differentiation matrix operators

Solves 𝛁²u(x,y) = f(x,y) 

May 2, 20200: v1:
    - Dirichlet boundary condition implementation 
    
May 5, 20200: v1_1:
    - Fixed a bug regarding the right-hand function
    - Figure size and font size adjusted
            
July 15, 2021: v1_2:
    - Added capabilities for defining Dirichlet boundary conditions at internal points

July 16, 2021: v1_3
    -Added Neumann boundary conditions
March 26, 2022: v1_4
    - Added ability to simulate multiple inner boundaries easily
    - Note that defining inner boundaries to Neumann leads to issues 
March 26, 2022: v1_4_1
    - Calculated the gradient of the solution, grad(u). |grad(u)| would be the electric field distribution if u is considered the electric potential.
    - Added contour plot and streamplot of the |grad(u)| function
    - Added colorset to plot in custom colors
March 26-30, 2022: v1_4_4
    - Fixed aliasing problem in contour plot export
    - 
May 3, 2022: v1_5
    - Generalized the boundary definition process. A separate function file is used to find all the outer and inner boundary indices.
      All boundary quantities are in list form. This makes it easier to implement all the boundary operations in one go (rather than
      treating each boundary separately).
                                                                                                                         
May 6, 2022: v_1_5_2
    - Created geo_source_definition file. Geometry and source definition functions have been moved there
    - Created utilities_definition file. Moved the plotting (and color) definition functions there    
"""


from __future__ import print_function    

import time
import math
import random
import numpy as np
import pylab as py
import scipy.sparse as sp                 # import sparse matrix library
from scipy.sparse.linalg import spsolve

# Change math font style and general font size
py.rcParams['mathtext.fontset'] = 'stix'
py.rcParams.update({'font.size': 20})





#==============================================================================
# Import functions from other files
#==============================================================================
# import the file where the differentiation matrix operators are defined
from diff_matrices import Diff_mat_1D, Diff_mat_2D   



from utilities_definition import my_contourf, my_scatter, color_distinct
from boundary_definition import boundary_regions
from geo_source_definition import geo_src_def


#==============================================================================




#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////




#===================================================================================================================
# Call geometry + source function. Call boundary region calculator function.
#===================================================================================================================

clr_set = color_distinct()   # load color set

prb_type = 'test_structure_3.bmp'
# prb_type = 'cap'

x,y, ub_o,B_type_o, xb_i,yb_i,ub_i,B_type_i,f =  geo_src_def(prb_type)


# Boundary definitions
if prb_type.find('.bmp') != -1:           # if the problem type is a bmp image file
    B_ind = xb_i
    B_type = B_type_i
    B_val = ub_i
else:                                     # if the problem type is a predefined structure and not a bmp image file
    B_ind, B_type, B_val = boundary_regions(x,y, ub_o,B_type_o, xb_i,yb_i,ub_i,B_type_i)


Nx = len(x)
Ny = len(y)

#===================================================================================================================



dx = x[1] - x[0]                # grid spacing along x direction
dy = y[1] - y[0]                # grid spacing along y direction

X,Y = np.meshgrid(x,y)          # 2D meshgrid

# 1D indexing
Xu = X.ravel()                  # Unravel 2D meshgrid to 1D array
Yu = Y.ravel()




#================================================================================================
# Loading finite difference matrix operators
#================================================================================================
Dx_2d, Dy_2d, D2x_2d, D2y_2d = Diff_mat_2D(Nx,Ny)   # Calling 2D matrix operators from funciton

#================================================================================================



#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////






#==============================================================================
# Plotting geometry and source regions
#==============================================================================

start_time = time.time()

# Find indices, types and value of all boundary regions


N_B = len(B_val)    # Number of different boundary regions. Each region consists of one or more boundary points (usually many many points for a well defined solution grid)

print("Boundary search time = %1.6s" % (time.time()-start_time))


# Plot solution domain (with boundary)
py.close('all')
py.figure(figsize = (9,7))
my_scatter(X,Y,clr_set[0],msize = 4)

for m in range(N_B):
    my_scatter(Xu[B_ind[m]], Yu[B_ind[m]],clr_set[N_B-m],msize = 4)

# Plot source function
py.figure(figsize = (9,7))
my_contourf(x,y,f.reshape(Ny,Nx),r'$f\,(x,y)$','RdBu')

#==============================================================================

#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////





#==============================================================================
# Construction of right hand vector (function of x and y)
#==============================================================================

start_time = time.time()


b = f
# Insert boundary values at the outer boundary points

for m in range(N_B):
    b[B_ind[m]] = B_val[m]



#==============================================================================





#============================================================================================
# Construction of the system matrix and adjust the right hand vector for boundary conditions
#============================================================================================

I_sp = sp.eye(Nx*Ny).tocsr()
L_sys = D2x_2d/dx**2 + D2y_2d/dy**2     # system matrix without boundary conditions

# Boundary operators
BD = I_sp       # Dirichlet boundary operator
BNx = Dx_2d     # Neumann boundary operator for x component
BNy = Dy_2d     # Neumann boundary operator for y component

# Selectively replace the rows of the system matrix that correspond to boundary value points. We replace these rows with 
# those of the boundary operator


for m in range(N_B):
    if B_type[m] == 0:
        L_sys[B_ind[m],:] = BD[B_ind[m],:]
    if B_type[m] == 1:
        L_sys[B_ind[m],:] = BNx[B_ind[m],:]
    if B_type[m] == 2:
        L_sys[B_ind[m],:] = BNy[B_ind[m],:]
        
    
#============================================================================================


#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////





print("System matrix and right hand vector computation time = %1.6s" % (time.time()-start_time)) 



#==============================================================================
# solve
start_time = time.time()
u = spsolve(L_sys,b).reshape(Ny,Nx)
print("spsolve() time = %1.6s" % (time.time()-start_time))
#==============================================================================

#==============================================================================
# Calculating the gradient of the solution
vx = -(Dx_2d*u.ravel()).reshape(Ny,Nx)
vy = -(Dy_2d*u.ravel()).reshape(Ny,Nx)
v = np.sqrt(vx**2 + vy**2)
#==============================================================================


#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////

#==============================================================================
# Plot solution


py.figure(figsize = (12,7))
my_contourf(x,y,u,r'$u\,(x,y)$')


# Plotting the gradient
py.figure(figsize = (12,7))
my_contourf(x,y,v,r'$|-\nabla u\,(x,y)|$','afmhot')
py.streamplot(x,y,vx,vy,color = 'w',density = 1.2, linewidth = 0.4)





# thin_factor = 10
# skip = (slice(None, None, thin_factor), slice(None, None, thin_factor))
# py.quiver(X[skip],Y[skip],vx[skip]/v[skip],vy[skip]/v[skip],color = 'w')


#==============================================================================





