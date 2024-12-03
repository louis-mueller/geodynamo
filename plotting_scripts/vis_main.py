# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:09:09 2024

@author: theresabuettner
"""

from input_data import params_units
from vis_functions import *
#from input_data import *

# BEFORE STARTING: install paramiko by pasting the following line into your terminal
# $ pip install paramiko

#%% LOCATION 

# If you want to access a specific folder, specify it in the following variables.
# If you want to access multiple folders, specify the parent folder.

# In case you are working locally, specify the folder you want to access in 'local_dir'.
# In case you are working on the server, specify the folder you want to access in 'server_dir'.  
# 
# curta local: '/scratch/louim95/'
# thor local: '/scratch/Simulations_Louis/'
# ideapad-louis local: 'D:\\Geodynamo_Core_Simulation\\tc_250924_M1-2_Fe30-60\\'  

local_dir    = 'D:\\Geodynamo_Core_Simulation\\tc_251024_M1-5_Fe30-60\\'
server_dir   = 'other'

# DATA TYPE

# What kind of data do you want to plot? You can also choose multiple types.

data_type = 'c'      # i - input profiles           
                     # h - heat stuff               
                     # b - parameters at boundaries
                     # c - core data ! ToDo implement this                

# ACCESS TO MULTIPLE FOLDERS

comparison = True    # Do you want to compare and plot data from more than one folder?

#%% GET DATA

all_data = get_data(data_type, local_dir, server_dir, comparison, params_units)

#%% PLOT DATA

# Put 'True' for all parameters you want to plot, and 'False' for the ones you don't want to plot.

# INPUT DATA  
params = {}

params['i'] = {'Gravity':                               True,
               'Pressure':                              True,
               'Density':                               True,
               'Temperature':                           True,
               'Melt Temperature':                      True,
               'Heat Capacity':                         True,
               'Thermal Expansivity':                   True,
               'Grüneisen Parameter':                   False,
               'Isothermal Bulk Modulus':               False,
               'Adiabatic Bulk Modulus':                False,
               'Shear Modulus':                         False,
               'Electrical Conductivity':               False,
               'Material Phase Number':                 False}
              
# HEAT DATA  

params['h'] = {'Heat Source':                           True, 
               'Mean Core Mantle Boundary Temperature': True, 
               'Core Mantle Boundary Temperature':      True, 
               'Surface Heatflow':                      True, 
               'Core Heatflow':                         True}

# BOUNDARY DATA

params['b'] = {'Radial Velocity':                       True, 
               'Radial Viscosity':                      True, 
               'Radial Temperature':                    True, 
               'Surface Heat Flux':                     True}

# CORE DATA
params['c'] = {'Time':                                  False,    
               'Inner Core Radius':                     True,
               'CMB Temperature':                       True,
               'CMB Heat':                              True,
               'Secular Cooling':                       True,
               'Latent Heat Release':                   True,
               'CMB Temperature Change':                False,
               'ICB Radius Change':                     False,
               'Thermal Buoyancy Flux':                 True,
               'Magnetic Moment':                       True,
               'CMB Magnetic Field Strength':           True,
               'Surface Magnetic Field Strength':       True}

plot_data(all_data, data_type, local_dir, server_dir, comparison, params, params_units)