#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:09:46 2024

@author: theresabuettner
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#%% INPUT & HEAT DATA

def colorcoding(foldernames):
    
    #nrColors = len(foldernames)
    nrColors = 5
    
    cmap = mpl.colormaps['viridis_r']
    colors = cmap(np.linspace(0,1,nrColors))
    
    return colors

#%% BOUNDARY DATA

def colorcoding_b():

    colormap    = {'Radial Velocity':          'summer', 
                   'Radial Viscosity':         'YlGnBu', 
                   'Radial Temperature':       'magma', 
                   'Surface Heat Flux':        'YlOrRd'}
    
    return colormap