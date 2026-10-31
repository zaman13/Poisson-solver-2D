#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 12:07:40 2022

@author: Mohammad Asif Zaman
"""


from __future__ import print_function    



import numpy as np
import pylab as py



# Defining custom plotting functions
def my_contourf(x,y,F,ttl,clrmp = 'inferno'):
    cnt = py.contourf(x,y,F,41,cmap = clrmp)
    
    # Antialiasing block for exporting figure to pdf later
    for c in cnt.collections:
        c.set_edgecolor("face")
    
    cbar = py.colorbar()
    py.xlabel(r'$x$',fontsize=26); py.ylabel(r'$y$',fontsize=26); 
    # py.title(ttl)
    cbar.set_label(ttl,fontsize=26)
    py.xlim([x[0],x[-1]])
    py.ylim([y[0],y[-1]])
    return 0
    

def my_scatter(x,y,clr,ttl='',msize = 2):
    py.plot(x,y,'.',markersize=msize,color=clr)
    py.xlabel(r'$x$',fontsize=26); py.ylabel(r'$y$',fontsize = 26); py.title(ttl)
    return 0


def color_distinct():
    
    # clr_set = ['#eaeee0','#ce9c9d','#adb5be','#57838d','#80adbc','#b4c9c7','#dadadc','#f3bfb3','#ccadb2','#445a67']
    #===================================================================================================================
    # clr_set = []
    # for j in range(20):
    #     temp= ["#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])]
    #     clr_set.extend(temp)
    
    # clr_set[0] = '#eaeee0'
    
    # https://sashamaps.net/docs/resources/20-colors/
    # clr_set = ['#eaeee0', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000','#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6']
    clr_set = ['#eaeee0','#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000']
    #===================================================================================================================

    return clr_set

