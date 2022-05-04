#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 12:28:10 2022

@author: Mohammad Asif Zaman

This script contains the boundary definition function. The boundary is defined by three quantities:
    1. The indices of the boundary points for each boundary region
    2. The boundary condition type at each boundary region
    3. The boundary value at each boundary region
    
    So, there should be 3 outputs.
    1. A list of numpy arrays for point 1 above
    2. A numpy array for for point 2 above
    3. A numpy array for for point 3 above
    
We may define the whole geometry of the problem by adding a few additional details:
    1. Defining the x,y independent variables
    2. Defining the f(x,y) source function    

This is not implemented in this script. A separate script will be developed for this. 


Input arguments:
    1. x,y                        : Numpy array 1D. 1D independent spatial variables
    2. ub_o                       : Boundary values at the outer boundary regions (order left, right, top, bottom)
    3. B_type_o                   : Boundary type at the outer boundary region
    4. xb_i, yb_i                 : 1D numpy arrays defining internal rectangular regions with boundary conditions enforced
    5. ub_i                       : Boundary values at the inner boundary regions
    6. B_type_i                   : Boundary types at the inner boundary regions


Output arguments:
    1. B_ind                      : Python list. Each element of the list contains indices of a specific boundary region. 
                                    All boundary regions (outer and inner) are covered by this list.
    2. B_type                     : Python list/numpy array 1D. The boundary type of each boundary region.                                      
    3. B_val                      : Python list/numpy array 1D. The boundary value of each boundary region.



"""

from __future__ import print_function    



import numpy as np




def boundary_regions(x,y, ub_o,B_type_o, xb_i,yb_i,ub_i,B_type_i):
    
    boundary_indices = []
    boundary_type = []
    boundary_value = []
    
    
    
    # Create necessary meshgrid and unravel grid from x,y. Define Nx, Ny
    # =========================================================================
    X,Y = np.meshgrid(x,y)          # 2D meshgrid

    # 1D indexing
    Xu = X.ravel()                  # Unravel 2D meshgrid to 1D array
    Yu = Y.ravel()
    
    Nx = len(x)
    Ny = len(y)
    
    # =========================================================================
    
    
    # Finding outer boundary regions
    # =========================================================================
    
    ind_unravel_L = np.squeeze(np.where(Xu==x[0]))          # Left boundary
    ind_unravel_R = np.squeeze(np.where(Xu==x[Nx-1]))       # Right boundary
    ind_unravel_B = np.squeeze(np.where(Yu==y[0]))          # Bottom boundary
    ind_unravel_T = np.squeeze(np.where(Yu==y[Ny-1]))       # Top boundary
    
    # ind_boundary = np.where((X==x[0]) | (X==x[Nx-1]) | (Y==y[0]) | (Y==y[Ny-1]))    # outer boundary
    # ind_boundary_unravel = np.squeeze(np.where((Xu==x[0]) | (Xu==x[Nx-1]) | (Yu==y[0]) | (Yu==y[Ny-1])))  # outer boundaries 1D unravel indices
    
    boundary_indices.append(ind_unravel_L)
    boundary_indices.append(ind_unravel_R)
    boundary_indices.append(ind_unravel_T)
    boundary_indices.append(ind_unravel_B)
    
    boundary_type.extend(B_type_o)
    boundary_value.extend(ub_o)
    
    # =========================================================================
    
    
    
    # Finding inner boundary regions
    # =========================================================================
    Nb_i = len(ub_i)                # Number of inner boundary regions 
    
    
    
    # ind_boundary_i = []           # empty lists (initialization)
    # ind_boundary_i_unravel = []   # empty lists (initialization)
    for m in range(Nb_i):
        temp1 = np.squeeze(np.where((Xu>xb_i[m][0]) & (Xu<xb_i[m][1]) & (Yu>yb_i[m][0]) & (Yu<yb_i[m][1])))  # inner boundaries defined by xb_i and yb_i
        # temp2 = np.where((X>xb_i[m][0]) & (X<xb_i[m][1]) & (Y>yb_i[m][0]) & (Y<yb_i[m][1]))    #  inner boundaries    
        boundary_indices.append(temp1)
        # ind_boundary_i.append(temp2)
    
    
    boundary_type.extend(B_type_i)
    boundary_value.extend(ub_i)
    
    
    
    
    
    return boundary_indices,boundary_type,boundary_value
    
    
    
    
    