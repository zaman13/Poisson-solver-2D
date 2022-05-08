#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 13:13:04 2022

@author: Mohammad Asif Zaman


A program to read an image file (pixel art) and convert it to a numpy matrix


*** All output indices are with respect to unravel meshgrid variables

"""


from __future__ import print_function    

import time
import math
import numpy as np
import pylab as py
from PIL import Image
import scipy.sparse as sp                 # import sparse matrix library
from scipy.sparse.linalg import spsolve

py.rcParams.update({'font.size': 20})





def im_str_read(im_file_name = "test_structure_3.bmp"):

    
    
    im = Image.open(im_file_name)
    
    #==============================================================================
    # Define independent variables
    # Nx =  No. of grid points along x direction
    # Ny = No. of grid points along y direction
    
    Nx,Ny = im.size
    
    x = np.linspace(-3,3,Nx)        # x variables in 1D
    y = np.linspace(-3,3,Ny)        # y variable in 1D
    
    # Note that the ranges of x and y are completely arbitrary here. These are local
    # dummy variables here. It's the indices that are important. The indices will
    # be the same for any range of x and y as long as Nx and Ny have the same
    # value. 
    
    #==============================================================================
    
    #//////////////////////////////////////////////////////////////////////////////
    #//////////////////////////////////////////////////////////////////////////////
    #//////////////////////////////////////////////////////////////////////////////
    
    #==============================================================================
    
    

    X,Y = np.meshgrid(x,y)          # 2D meshgrid
    # 1D indexing
    Xu = X.ravel()                  # Unravel 2D meshgrid to 1D array
    Yu = Y.ravel()
    
    #==============================================================================
    
    # *******
    # We are considering 16 color bmp image. So, color values range from 0 to 15
    # For color values less than boundary_region_color_threshold, we consider the 
    # region to be a boundary region. If the color value is larger than that, we 
    # consider that a solution region. Note that there might be multiple solution region
    # depending on how many color values above boundary_region_color_threshold we have 
    # used. These distinct solution regions can be used to define source functions
    # that vary depending on the region. 
    
    
    boundary_region_color_threshold = 10
    # *******
    #==============================================================================
    
    
    
    strct_map = np.flipud(np.array(im))   # Read the image value. The flipud is necessary to match orientation.
    
    
    strct_map_u = strct_map.ravel()
    
    regions = np.unique(strct_map_u)
    n_regions = np.size(regions)
    print('Number of regions identified = %d ' % (n_regions) )
    
    
    # Boundary regions
    # b_regions_indices = []
    b_regions_indices_u = []
    b_regions_size = []
    
    # Other regions
    # o_regions_indices = []
    o_regions_indices_u = []
    o_regions_size = []
    
    
    counter_b = 0
    counter_o = 0
    
    
    
    for m in range(n_regions):
        
        ind_temp = np.where(strct_map==regions[m])
        ind_u_temp = np.squeeze(np.where(strct_map_u==regions[m])) 
        
        if regions[m] <= boundary_region_color_threshold:  # boundary region
            # b_regions_indices.append(ind_temp)    
            b_regions_size.append(len(ind_u_temp))
            b_regions_indices_u.append(ind_u_temp)
            
        else:
            # o_regions_indices.append(ind_temp)
            o_regions_size.append(len(ind_u_temp))
            o_regions_indices_u.append(ind_u_temp)
        
            
    # Note that we read all the non-boundary indices distinctly. This would be useful in future if we want to define a complex
    # source function with different value in different non-boundary regions.         
        
        
    
    n_b_regions = len(b_regions_size)
    n_o_regions = len(o_regions_size)
    print('Number of boundary regions = %d ' % (n_b_regions) )
    
    # clr_set = ['#eaeee0','#ce9c9d','#adb5be','#57838d','#80adbc','#b4c9c7','#dadadc','#f3bfb3','#ccadb2','#445a67']
    
    # my_scatter(X,Y,'k','Solution grid')
    
    # for m in range(n_b_regions):
    #     ind1 = b_regions_indices[m]
    #     clr = clr_set[m] if m < n_regions-1 else "#ffffff"
        
        
    #     my_scatter(X[ind1], Y[ind1],clr,'Boundary regions')
        
    return Nx,Ny, b_regions_indices_u, o_regions_indices_u
    
    # my_scatter(X[ind_boundary2], Y[ind_boundary2],'b','')
    
    
