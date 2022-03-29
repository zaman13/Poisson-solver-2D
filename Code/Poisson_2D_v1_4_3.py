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
    

"""


from __future__ import print_function    

import time
import math
import numpy as np
import pylab as py
import scipy.sparse as sp                 # import sparse matrix library
from scipy.sparse.linalg import spsolve

# Change math font style and general font size
py.rcParams['mathtext.fontset'] = 'stix'
py.rcParams.update({'font.size': 20})


# import the file where the differentiation matrix operators are defined
from diff_matrices import Diff_mat_1D, Diff_mat_2D   


# Defining custom plotting functions
def my_contourf(x,y,F,ttl,clrmp = 'inferno'):
    cnt = py.contourf(x,y,F,41,cmap = clrmp)
    cbar = py.colorbar()
    py.xlabel(r'$x$'); py.ylabel(r'$y$'); 
    # py.title(ttl)
    cbar.set_label(ttl)
    return 0
    

def my_scatter(x,y,clr,ttl=''):
    py.plot(x,y,'.',markersize=2,color=clr)
    py.xlabel(r'$x$'); py.ylabel(r'$y$'); py.title(ttl)
    return 0
    



clr_set = ['#eaeee0','#ce9c9d','#adb5be','#57838d','#80adbc','#b4c9c7','#dadadc','#f3bfb3','#ccadb2','#445a67']




#==============================================================================
# Dirichlet/Neumann boundary conditions at outerwalls (boundary condition type is defined through boundary operators)
uL = 0
uR = 0
uT = 0
uB = 0
B_type_o = [1,1,1,1]    # Type of outer boundary conditions. 0 = Dirichlet, 1 = Neumann 
                        # Order element 1 = left, element 2 = right, element 3 = top, element 4 = bottom
#==============================================================================



#==============================================================================
# Dirichlet boundary conditions at an inner rectangular region

# Boundary defintion for the capacitor example
# ub_i = [2,-2]                      # boundary values at inner region

# xb_i = [[-2, 2],[-2, 2]]           # lower and upper limits of x defining the inner boundary region Format:[ [xlow1,xhigh1], [xlow2,xhigh2], .... ]
# yb_i = [[.5,.6],[-.6, -.5]]        # lower and upper limits of y defining the inner boundary region Format:[ [ylow1,yhigh1], [ylow2,yhigh2], .... ]

# Nb_i = len(ub_i)                   # Number of inner boundary regions
# B_type_i = [0, 0]                  # Type of inner boundary conditions. 0 = Dirichlet, 1 = Neumann. Note that setting inner boundaries to Neumann can lead to issues



# Boundary definitions for the diode example
ub_i = [-2,2]                    # boundary values at inner region

xb_i = [[-4,-3.8],[3.8,4]]       # lower and upper limits of x defining the inner boundary region Format:[ [xlow1,xhigh1], [xlow2,xhigh2], .... ]
yb_i = [[-1,1],[-1,1]]           # lower and upper limits of y defining the inner boundary region Format:[ [ylow1,yhigh1], [ylow2,yhigh2], .... ]

Nb_i = len(ub_i)                 # Number of inner boundary regions
B_type_i = [0,0]                 # Type of inner boundary conditions. 0 = Dirichlet, 1 = Neumann. Note that setting inner boundaries to Neumann can lead to issues

#==============================================================================



#==============================================================================
# Define independent variables
Nx = 300                         # No. of grid points along x direction
Ny = 200                         # No. of grid points along y direction
x = np.linspace(-6,6,Nx)        # x variables in 1D
y = np.linspace(-3,3,Ny)        # y variable in 1D
#==============================================================================

#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////

#==============================================================================


dx = x[1] - x[0]                # grid spacing along x direction
dy = y[1] - y[0]                # grid spacing along y direction

X,Y = np.meshgrid(x,y)          # 2D meshgrid

# 1D indexing
Xu = X.ravel()                  # Unravel 2D meshgrid to 1D array
Yu = Y.ravel()


# Source function (right hand side vector)
f = np.zeros(Nx*Ny) 

# Source term for the diode example
for m in range(Nx*Ny):
    if np.abs(Xu[m]) < 3 and np.abs(Yu[m])<1:
        f[m] = -np.tanh(Xu[m])/np.cosh(Xu[m])
        

# Loading finite difference matrix operators

Dx_2d, Dy_2d, D2x_2d, D2y_2d = Diff_mat_2D(Nx,Ny)   # Calling 2D matrix operators from funciton



#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////


#==============================================================================
# Boundary indices
start_time = time.time()
ind_unravel_L = np.squeeze(np.where(Xu==x[0]))          # Left boundary
ind_unravel_R = np.squeeze(np.where(Xu==x[Nx-1]))       # Right boundary
ind_unravel_B = np.squeeze(np.where(Yu==y[0]))          # Bottom boundary
ind_unravel_T = np.squeeze(np.where(Yu==y[Ny-1]))       # Top boundary

ind_boundary_unravel = np.squeeze(np.where((Xu==x[0]) | (Xu==x[Nx-1]) | (Yu==y[0]) | (Yu==y[Ny-1])))  # outer boundaries 1D unravel indices
ind_boundary = np.where((X==x[0]) | (X==x[Nx-1]) | (Y==y[0]) | (Y==y[Ny-1]))    # outer boundary


ind_boundary_i = []           # empty lists (initialization)
ind_boundary_i_unravel = []   # empty lists (initialization)

for m in range(Nb_i):
    temp1 = np.squeeze(np.where((Xu>xb_i[m][0]) & (Xu<xb_i[m][1]) & (Yu>yb_i[m][0]) & (Yu<yb_i[m][1])))  # inner boundaries defined by xb_i and yb_i
    temp2 = np.where((X>xb_i[m][0]) & (X<xb_i[m][1]) & (Y>yb_i[m][0]) & (Y<yb_i[m][1]))    #  inner boundaries
    
    ind_boundary_i_unravel.append(temp1)
    ind_boundary_i.append(temp2)

print("Boundary search time = %1.6s" % (time.time()-start_time))


# Plot solution domain (with boundary)
py.close('all')
py.figure(figsize = (12,7))
my_scatter(X,Y,clr_set[0],'Solution grid')
my_scatter(X[ind_boundary], Y[ind_boundary],clr_set[1],'Solutiohn grid with boundary')
for m in range(Nb_i):
    my_scatter(X[ind_boundary_i[m]], Y[ind_boundary_i[m]],clr_set[2])


# Plot source function
py.figure(figsize = (12,7))
my_contourf(x,y,f.reshape(Ny,Nx),r'$f\,(x,y)$','RdBu')

#==============================================================================

#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////





#==============================================================================
# Construction of right hand vector (function of x and y)
start_time = time.time()

b = f
# Insert boundary values at the outer boundary points
b[ind_unravel_L] = uL
b[ind_unravel_R] = uR
b[ind_unravel_T] = uT
b[ind_unravel_B] = uB


#==============================================================================


#==============================================================================
# Construction of the system matrix and adjust the right hand vector for boundary conditions
I_sp = sp.eye(Nx*Ny).tocsr()
L_sys = D2x_2d/dx**2 + D2y_2d/dy**2     # system matrix without boundary conditions

# Boundary operators
BD = I_sp       # Dirichlet boundary operator
BNx = Dx_2d     # Neumann boundary operator for x component
BNy = Dy_2d     # Neumann boundary operator for y component

# Selectively replace the rows of the system matrix that correspond to boundary value points. We replace these rows with 
# those of the boundary operator

# L_sys[ind_boundary_unravel,:] = BD[ind_boundary_unravel,:]    # Boundaries at the four edges


# Handle inner boundary condtions
for m in range(Nb_i):
    # Boundaries defined by xb_i and yb_i (inside boundaries)
    
    b[ind_boundary_i_unravel[m]] = ub_i[m]   # Fixing the right hand vector
    
    # Fixing the system matrix
    if B_type_i[m] == 0:    # For Dirichlet boundary conditions
        L_sys[ind_boundary_i_unravel[m],:] = BD[ind_boundary_i_unravel[m],:]  
    else:                # For Neumann boundary conditions
        L_sys[ind_boundary_i_unravel[m],:] = BNx[ind_boundary_i_unravel[m],:]  
            
            




L_sys[ind_unravel_L,:] = BD[ind_unravel_L,:] if B_type_o[0] == 0 else BNx[ind_unravel_L,:]    # Boundaries at the left layer
L_sys[ind_unravel_R,:] = BD[ind_unravel_R,:] if B_type_o[1] == 0 else BNx[ind_unravel_R,:]    # Boundaries at the right edges
L_sys[ind_unravel_T,:] = BD[ind_unravel_T,:] if B_type_o[2] == 0 else BNy[ind_unravel_T,:]    # Boundaries at the top layer
L_sys[ind_unravel_B,:] = BD[ind_unravel_B,:] if B_type_o[3] == 0 else BNy[ind_unravel_B,:]    # Boundaries at the bottom layer

#==============================================================================

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





