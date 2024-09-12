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

local_dir    = '\Users\louis\Documents\plots'
server_dir   = '/scratch/Simulations_Louis'


# DATA TYPE

# What kind of data do you want to plot? You can also choose multiple types.

data_type = 'h'    # i - input profiles           
                     # h - heat stuff               
                     # b - parameters at boundaries                

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
               'Grüneisen Parameter':                   True,
               'Isothermal Bulk Modulus':               True,
               'Adiabatic Bulk Modulus':                True,
               'Shear Modulus':                         True,
               'Electrical Conductivity':               True,
               'Material Phase Number':                 True}
              
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

plot_data(all_data, data_type, local_dir, server_dir, comparison, params, params_units)
