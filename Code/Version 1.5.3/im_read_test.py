#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 21:55:34 2022

@author: asif
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


from image_structure_read import im_str_read
from utilities_definition import my_contourf, my_scatter, color_distinct


Nx,Ny,b_regions_indices, o_regions_indices = im_str_read('test_structure_3.bmp')

x = np.linspace(-3,3,Nx)        # x variables in 1D
y = np.linspace(-3,3,Ny)        # y variable in 1D





X,Y = np.meshgrid(x,y)          # 2D meshgrid
# 1D indexing
Xu = X.ravel()                  # Unravel 2D meshgrid to 1D array
Yu = Y.ravel()


# my_scatter(Xu,Yu,'k','Solution grid')
# clr_set = ['#eaeee0','#ce9c9d','#adb5be','#57838d','#80adbc','#b4c9c7','#dadadc','#f3bfb3','#ccadb2','#445a67']
clr_set = ['r','k','b','g','#80adbc','#b4c9c7','#dadadc','#f3bfb3','#ccadb2','#445a67']

for m in range(len(b_regions_indices)):
    ind1 = b_regions_indices[m]
    clr = clr_set[m] 
    
    
    my_scatter(Xu[ind1], Yu[ind1],clr,'Boundary regions')


# for m in range(len(o_regions_indices)):
#     ind1 = o_regions_indices[m]
#     clr = clr_set[m + len(b_regions_indices) ] 
    
    
#     my_scatter(Xu[ind1], Yu[ind1],clr,'Boundary regions')