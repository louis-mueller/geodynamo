#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:09:35 2024

@author: theresabuettner
"""

from input_data import loc
from getdata import *
#from location import *
from coloring import *
from plotting import *


#%% GET DATA

def get_data(data_type, local_dir, server_dir, comparison, params_units):
    
    files(data_type)
    
    if comparison == True:
        findfolders(loc, local_dir, server_dir)

    if loc == 'r': 
        getfiles(comparison, local_dir, server_dir)
        
    if loc == 's':
        homedir = server_dir
        
    if loc == 'r' or loc == 'l':
        homedir = local_dir
        
    all_data = importfiles(data_type, homedir, comparison, params_units)
    
    return all_data

#%% PLOT DATA

def plot_data(all_data, data_type, local_dir, server_dir, comparison, params, params_units):
    
    if loc == 's':
        homedir = server_dir
        
    if loc == 'r' or loc == 'l':
        homedir = local_dir
    
    if 'i' in data_type:
        plot_input(all_data, homedir, comparison, params, params_units)
    
    if 'h' in data_type:
        plot_heat(all_data, homedir, comparison, params, params_units)
    
    if 'b' in data_type:
        plot_2D(all_data, homedir, comparison, params, params_units)