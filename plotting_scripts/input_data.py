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
username      = 'tbuettner'
password      = 'Thor-Al-Paka1213'
#password     = getpass(" password: ") 

safety        = ['username', 'password', 'hostip']

#%% DEFINE CASES

loc = 'l'            # Where are you working?
                     # l - local  (files are saved locally on the computer you are working at)
                     # r - remote (files are saved on a remote server and need to be downloaded)
                     # s - server (files are saved on server, and you are working on the server)


#%% FILENAMES

# Define the filenames of the differend output files you want to access.

input_files    =  '/profs.res'             # input profiles
heat_files     =  '/data_char_val.res'     # heat stuff
boundary_files = ['/data_vel_profs.res',   # boundaries:   velocity
                  '/data_visc_profs.res',                # viscosity
                  '/data_temp_profs.res',                # temperature
                  '/data_surf_heatflux_profs.res']       # surface heat flux

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
                     'Thermal Expansivity':                   r'$\alpha \ [10^{-5} \ \mathrm{s^{-1}}]$',
                     'Grüneisen Parameter':                   r'$\gamma$',
                     'Isothermal Bulk Modulus':               r'$K_T \ [\mathrm{GPa}]$',
                     'Adiabatic Bulk Modulus':                r'$K_S \ [\mathrm{GPa}]$',
                     'Shear Modulus':                         r'$G_S \ [\mathrm{GPa}]$',
                     'Electrical Conductivity':               r'$\sigma \ [\mathrm{S/m}]$',
                     'Material Phase Number':                 r'$mat$'}
              
# HEAT DATA  

params_units['h'] = {'Time':                                  r'$t \ [\mathrm{Myr}]$',
                     'Heat Source':                           r'$ ??? $', 
                     'Mean Core Mantle Boundary Temperature': r'$T \ [\mathrm{K}]$', 
                     'Core Mantle Boundary Temperature':      r'$T \ [\mathrm{K}]$', 
                     'Surface Heatflow':                      r'$SHF \ [\mathrm{mW/m^2}]$', 
                     'Core Heatflow':                         r'$CHF \ [\mathrm{mW/m^2}]$'}

# BOUNDARY DATA

params_units['b'] = {'Radial Velocity':                       r'$v \ [\mathrm{m/s}]$', 
                     'Radial Viscosity':                      r'$\eta \ [\mathrm{mPas}]$', 
                     'Radial Temperature':                    r'$T \ [\mathrm{K}]$', 
                     'Surface Heat Flux':                     r'$SHF \ [\mathrm{mW/m^2}]$'}


#%% PLOTTING

smoothing       = False
show_rawdata    = True
spacing         = 5

