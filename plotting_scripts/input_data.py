#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:18:03 2024

@author: theresabuettner
"""

from getpass import getpass

#%% REMOTE ACCESS 

# In case you are working remotely, fill out the login and server information below.
# Please make sure you are in the local FU network or connected by VPN before logging into the server.

# If you are working locally or on the server, you can leave the strings empty.

hostip        = 'thor.planet.fu-berlin.de'    # alternatively '130.133.87.145' for Freya
username      = 'louismueller'
password      = 'hi'

safety        = ['username', 'password', 'hostip']

#%% DEFINE CASES

loc = 'l'            # Where are you working?
                     # l - local  (files are saved locally on the computer you are working locally)
                     # r - remote (files are saved on a remote server and need to be downloaded)
                     # s - server (files are saved on server, and you are working on the server)


#%% FILENAMES

# Define the filenames of the differend output files you want to access.

input_files    =  'profs.res'             # input profiles
heat_files     =  'data_char_val.res'     # heat stuff
boundary_files = ['data_vel_profs.res',   # boundaries:   velocity
                  'data_visc_profs.res',                # viscosity
                  'data_temp_profs.res',                # temperature
                  'data_surf_heatflux_profs.res']       # surface heat flux
core_files     =  'data_core.res'          # core data

#%% DATASETS 

params_units = {}

# INPUT DATA  

params_units['i'] = {'Gravity':                               r'$g \ [\mathrm{m/s^2}]$',
                     'Pressure':                              r'$p \ [\mathrm{GPa}]$',
                     'Density':                               r'$\rho \ [\mathrm{kg/m^3}]$',
                     'Radius':                                r'$r \ [\mathrm{km}]$',
                     'Temperature':                           r'$T \ [\mathrm{K}]$',
                     'Melt Temperature':                      r'$T_m \ [\mathrm{K}]$',
                     'Heat Capacity':                         r'$C_p \ [\mathrm{J/kgK}]$',
                     'Thermal Expansivity':                   r'$\alpha \ [10^{-5}/s]$',
                     'Grüneisen Parameter':                   r'$\gamma$',
                     'Isothermal Bulk Modulus':               r'$K_T \ [\mathrm{GPa}]$',
                     'Adiabatic Bulk Modulus':                r'$K_S \ [\mathrm{GPa}]$',
                     'Shear Modulus':                         r'$G_S \ [\mathrm{GPa}]$',
                     'Electrical Conductivity':               r'$\sigma \ [\mathrm{S/m}]$',
                     'Material Phase Number':                 r'$mat$'}
              
# HEAT DATA  

params_units['h'] = {'Time':                                  r'$t \ [\mathrm{Myr}]$',
                     'Heat Source':                           r'$h \ [\mathrm{pW/kg}] $', 
                     'Mean Core Mantle Boundary Temperature': r'$T \ [\mathrm{K}]$', 
                     'Core Mantle Boundary Temperature':      r'$T \ [\mathrm{K}]$', 
                     'Surface Heatflow':                      r'$q_s \ [\mathrm{mW/m^2}]$', 
                     'Core Heatflow':                         r'$q_c \ [\mathrm{mW/m^2}]$'}

# BOUNDARY DATA

params_units['b'] = {'Radial Velocity':                       r'$v \ [\mathrm{m/s}]$', 
                     'Radial Viscosity':                      r'$\eta \ [\mathrm{mPas}]$', 
                     'Radial Temperature':                    r'$T \ [\mathrm{K}]$', 
                     'Surface Heat Flux':                     r'$q_s \ [\mathrm{mW/m^2}]$'}

# CORE DATA

params_units['c'] = {'Time':                                  r'$t \ [\mathrm{Myr}]$',
                     'Inner Core Volume':                     r'$V_i/V_c \ [\%]$',
                     'CMB Temperature':                       r'$T_c \ [\mathrm{K}]$',
                     'CMB Heat':                              r'$Q_{CMB} \ [\mathrm{TW}]$',
                     'Secular Cooling':                       r'$Q_S \ [\mathrm{TW}]$',
                     'Latent Heat Release':                   r'$Q_L \ [\mathrm{TW}]$',
                     'CMB Temperature Change':                r'$\Delta T \ [\mathrm{K/Myr}]$',
                     'ICB Radius Change':                     r'$\Delta r \ [\mathrm{m/Myr}]$',
                     'Thermal Buoyancy Flux':                 r'$F_T \ [\mathrm{m^2/s^3}]$',
                     'Magnetic Moment':                       r'$m \ [\mathrm{A m^2}]$',
                     'CMB Magnetic Field Strength':           r'$B_{CMB} \ [\mathrm{\mu T}]$',
                     'Surface Magnetic Field Strength':       r'$B_{S} \ [\mathrm{\mu T}]$'}



#%% PLOTTING

smoothing       = False
show_rawdata    = True
spacing         = 5
save_fig        = True

# only used if save_fig is set to True
#--------------------------------------------------------------------------------
output_name     = 'M3-4_Fe30-60'
file_type       = '.svg'        # will only effect i, h, and c output.
plot_dir        = '/scratch/Simulations_Louis/Plots/'
#--------------------------------------------------------------------------------
