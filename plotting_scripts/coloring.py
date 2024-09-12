#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:09:46 2024

@author: theresabuettner
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

colormap = {}

#%% INPUT DATA

def colorcoding(foldernames):
    
    nrColors = len(foldernames)
    
    cmap = mpl.colormaps['viridis']
    colors = cmap(np.linspace(0,1,nrColors))
    
    return colors

#%% HEAT DATA

colormap['h'] = {}


#%% BOUNDARY DATA

colormap    = {'Radial Velocity':          'summer', 
               'Radial Viscosity':         'YlGnBu', 
               'Radial Temperature':       'magma', 
               'Surface Heat Flux':        'YlOrRd'}