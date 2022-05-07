#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 12:07:40 2022

@author: Mohammad Asif Zaman
"""


from __future__ import print_function    



import numpy as np




def geo_src_def(geo_preset):
    
    #==============================================================================
    # Define independent variables
    Nx = 300                         # No. of grid points along x direction
    Ny = 200                         # No. of grid points along y direction
    x = np.linspace(-6,6,Nx)         # x variables in 1D
    y = np.linspace(-3,3,Ny)         # y variable in 1D
    #==============================================================================
    
        
       
    #================================================================================================================================
    # Dirichlet/Neumann boundary conditions at outerwalls (boundary condition type is defined through boundary operators)
    #================================================================================================================================
    uL = 0
    uR = 0
    uT = 0
    uB = 0
    ub_o = [uL, uR, uT, uB]
    B_type_o = [1,1,2,2]    # Type of outer boundary conditions. 0 = Dirichlet, 1 = Neumann x derivative, 2 = Neumann y derivative
                            # Order element 1 = left, element 2 = right, element 3 = top, element 4 = bottom
    
    # We have put Neumann boundaries at the outerwalls. This should be fine for most cases.
    #================================================================================================================================
    
    
    
    
    #================================================================================================================================
    # Variable range and Dirichlet boundary conditions at an inner rectangular region. This depends on the example problem we are
    # trying to solve
    #================================================================================================================================
    
    # Boundary defintion for the capacitor example
    #==============================================================================
    if geo_preset == 'cap':
        
        #==============================================================================
        # Define independent variables
        Nx = 300                         # No. of grid points along x direction
        Ny = 200                         # No. of grid points along y direction
        x = np.linspace(-6,6,Nx)         # x variables in 1D
        y = np.linspace(-3,3,Ny)         # y variable in 1D
        #==============================================================================
        
        
        ub_i = [2,-2]                      # boundary values at inner region
        
        xb_i = [[-2, 2],[-2, 2]]           # lower and upper limits of x defining the inner boundary region Format:[ [xlow1,xhigh1], [xlow2,xhigh2], .... ]
        yb_i = [[.5,.6],[-.6, -.5]]        # lower and upper limits of y defining the inner boundary region Format:[ [ylow1,yhigh1], [ylow2,yhigh2], .... ]
       
        Nb_i = len(ub_i)                   # Number of inner boundary regions
        B_type_i = [0, 0]                  # Type of inner boundary conditions. 0 = Dirichlet, 1 = Neumann x derivative, 2 = Neumann y derivative. Note that setting inner boundaries to Neumann can lead to issues
        
        # Source function (right hand side vector)
        f = np.zeros(Nx*Ny) 
    #==============================================================================
    
    
    
    
    
    # Boundary definitions for the diode example
    #==============================================================================
    if geo_preset == 'diode':
        
        #==============================================================================
        # Define independent variables
        Nx = 300                         # No. of grid points along x direction
        Ny = 200                         # No. of grid points along y direction
        x = np.linspace(-6,6,Nx)         # x variables in 1D
        y = np.linspace(-3,3,Ny)         # y variable in 1D
        
        X,Y = np.meshgrid(x,y)          # 2D meshgrid

        # 1D indexing
        Xu = X.ravel()                  # Unravel 2D meshgrid to 1D array
        Yu = Y.ravel()

        #==============================================================================
        
        ub_i = [-2,2]                    # boundary values at inner region
        
        xb_i = [[-4,-3.8],[3.8,4]]       # lower and upper limits of x defining the inner boundary region Format:[ [xlow1,xhigh1], [xlow2,xhigh2], .... ]
        yb_i = [[-1,1],[-1,1]]           # lower and upper limits of y defining the inner boundary region Format:[ [ylow1,yhigh1], [ylow2,yhigh2], .... ]
        
        Nb_i = len(ub_i)                 # Number of inner boundary regions
        B_type_i = [0,0]                 # Type of inner boundary conditions. 0 = Dirichlet, 1 = Neumann. Note that setting inner boundaries to Neumann can lead to issues
        
        # Source term for the diode example
        f = np.zeros(Nx*Ny) 
        
        for m in range(Nx*Ny):
            if np.abs(Xu[m]) < 3.8 and np.abs(Yu[m])<1:
                f[m] = -np.tanh(Xu[m])/np.cosh(Xu[m])

    #==============================================================================
    
    # Boundary definition for the gravitational problem
    # ub_i = []                   # boundary values at inner region
    
    # xb_i = []                   # lower and upper limits of x defining the inner boundary region Format:[ [xlow1,xhigh1], [xlow2,xhigh2], .... ]
    # yb_i = []                   # lower and upper limits of y defining the inner boundary region Format:[ [ylow1,yhigh1], [ylow2,yhigh2], .... ]
    
    # Nb_i = len(ub_i)            # Number of inner boundary regions
    # B_type_i = []  
    
           
    # Source term for the gravitational potential example
    # for m in range(Nx*Ny):
    #     if (Xu[m]-1.2)**2 + (Yu[m]+0.5)**2 < .3**2:
    #         f[m] = .1
    #     if (Xu[m]+1.2)**2 + (Yu[m]-1.5)**2 < .3**2:
    #         f[m] = .05      



    
    return x,y, ub_o,B_type_o, xb_i,yb_i,ub_i,B_type_i, f
    
    #==============================================================================
    
     
        
        