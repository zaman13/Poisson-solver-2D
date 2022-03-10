#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 19:42:11 2022

@author: Mohammad ASif Zaman

A program to read an image file (pixel art) and convert it to a numpy matrix.
Different color regions of the image file will be segmented as different regions.
The indices corresponding to each region (with respec to a meshgrid array) will
be calculated. The function returns the indices of each region as a list file.


% Output variables


Nx                :    Number of elements in the x direction
Ny                :    Number of elements in the y direction
n_regions         :    Number of identified regions
regions_indices   :    list of indices of each region. (e.g. regions_indices[1] = indices of all elements that fall in region 1)
regiond_indices_u :    list of indices in unravel form

 
"""


import numpy as np
#==============================================================================
def structure_definition(im):

    Nx,Ny = im.size
    
    # Dummy variables having the same size as the acutal spatial variables but spanning from -1 to 1
    # The boundary indices found for this specific dummy space variables would be the same as the boundary indices
    # of the actual spatial variables.
    x = np.linspace(-1,1,Nx)        # x variables in 1D
    y = np.linspace(-1,1,Ny)        # y variable in 1D


    # Meshgrid of the dummy spatial variables
    X,Y = np.meshgrid(x,y)          # 2D meshgrid
    # 1D indexing of the dummy spatial variables
    
    Xu = X.ravel()                  # Unravel 2D meshgrid to 1D array
    Yu = Y.ravel()
    

    # Read the image file in a numpy array    
    strct_map = np.flipud(np.array(im))   # the flipud is necessary to match the numpy array coordinate with the figure pixel coordinates
    
    # Map the structure in 1D (unravel)
    strct_map_u = strct_map.ravel()
    
    # Find the unique elements of structure array. Each unique element referst to a specific color.
    # When drawing the image, different boundary regions should be drawn with different colors. 
    # All regions having the same color would be assigned the same boundary value (or will be assigned as unknown)
    
    regions = np.unique(strct_map_u)
    
    n_regions = np.size(regions)  # The number of unique regions
    
    
    print('Number of regions identified = %d \n' % (n_regions) )
    
    
    # Defining empty lists. They will be populated by indices of different boundary regions
    regions_indices = [];
    regions_indices_u = [];
    
    
    for m in range(n_regions):
        ind_temp = np.where(strct_map==regions[m])       # indices for the meshgrid
        ind_u_temp = np.squeeze(np.where(strct_map_u==regions[m])) # indices for the unravel case
        
        regions_indices.append(ind_temp)         # store meshgrid indices
        regions_indices_u.append(ind_u_temp)     # store unravel indices



    return Nx,Ny,n_regions, regions_indices, regions_indices_u


#==============================================================================