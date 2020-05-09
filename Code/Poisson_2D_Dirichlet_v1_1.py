# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 13:03:24 2020

@author: Mohammad Asif Zaman

Code for testng differentiation matrix operators

Solves 𝛁²u(x,y) = g(x,y) 

May 2, 20200: v1:
    - Dirichlet boundary condition implementation 
    
May 5, 20200: v1_1:
    - Fixed a bug regarding the right-hand function
    - Figure size and font size adjusted
            


"""


from __future__ import print_function    

import time
import math
import numpy as np
import pylab as py
import scipy.sparse as sp                 # import sparse matrix library
from scipy.sparse.linalg import spsolve

py.rcParams.update({'font.size': 20})

# import the file where the differentiation matrix operators are defined
from diff_matrices import Diff_mat_1D, Diff_mat_2D   


# Defining custom plotting functions
def my_contourf(x,y,F,ttl):
    cnt = py.contourf(x,y,F,41,cmap = 'inferno')
    py.colorbar()
    py.xlabel('x'); py.ylabel('y'); py.title(ttl)
    return 0
    

def my_scatter(x,y,clr,ttl):
    py.plot(x,y,'.',markersize=2,color=clr)
    py.xlabel('x'); py.ylabel('y'); py.title(ttl)
    return 0
    

# Dirichlet boundary conditions
uL = 0.5
uR = 1.2
uT = 1.5
uB = 0



# Define independent variables
Nx = 50                         # No. of grid points along x direction
Ny = 80                         # No. of grid points along y direction
x = np.linspace(-6,6,Nx)        # x variables in 1D
y = np.linspace(-3,3,Ny)        # y variable in 1D

dx = x[1] - x[0]                # grid spacing along x direction
dy = y[1] - y[0]                # grid spacing along y direction

X,Y = np.meshgrid(x,y)          # 2D meshgrid

# 1D indexing
Xu = X.ravel()                  # Unravel 2D meshgrid to 1D array
Yu = Y.ravel()

# Source function (right hand side vector)
g = np.zeros(Nx*Ny) 


 

# Loading finite difference matrix operators

Dx_2d, Dy_2d, D2x_2d, D2y_2d = Diff_mat_2D(Nx,Ny)   # Calling 2D matrix operators from funciton



# Boundary indices
start_time = time.time()
ind_unravel_L = np.squeeze(np.where(Xu==x[0]))          # Left boundary
ind_unravel_R = np.squeeze(np.where(Xu==x[Nx-1]))       # Right boundary
ind_unravel_B = np.squeeze(np.where(Yu==y[0]))          # Bottom boundary
ind_unravel_T = np.squeeze(np.where(Yu==y[Ny-1]))       # Top boundary

ind_boundary_unravel = np.squeeze(np.where((Xu==x[0]) | (Xu==x[Nx-1]) | (Yu==y[0]) | (Yu==y[Ny-1])))  # All boundary
ind_boundary = np.where((X==x[0]) | (X==x[Nx-1]) | (Y==y[0]) | (Y==y[Ny-1]))    # All boundary
print("Boundary search time = %1.6s" % (time.time()-start_time))


# Plot solution domain (with boundary)

py.close('all')
my_scatter(X,Y,'b','Solution grid')
my_scatter(X[ind_boundary], Y[ind_boundary],'r','Solutiohn grid with boundary')



# Construction of the system matrix
start_time = time.time()
I_sp = sp.eye(Nx*Ny).tocsr()
L_sys = D2x_2d/dx**2 + D2y_2d/dy**2     # system matrix without boundary conditions

L_sys[ind_boundary_unravel,:] = I_sp[ind_boundary_unravel,:]


# Construction of right hand vector (function of x and y)
b = g
b[ind_unravel_L] = uL
b[ind_unravel_R] = uR
b[ind_unravel_T] = uT
b[ind_unravel_B] = uB

print("System matrix and right hand vector computation time = %1.6s" % (time.time()-start_time)) 


# solve
start_time = time.time()
u = spsolve(L_sys,b).reshape(Ny,Nx)
print("spsolve() time = %1.6s" % (time.time()-start_time))



# Plot solution
py.figure(figsize = (14,7))
my_contourf(x,y,u,r'$\nabla^2 u = 0$')





