#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 10:47:55 2024

@author: theresabuettner
"""

from input_data import *
for name in safety:
    del globals()[name]

import paramiko
import numpy as np

#%% DEFINE FILES

def files(data_type):
    
    global filenames
    filenames = []
    
    if 'i' in data_type:
        filenames.append(input_files)
        
    if 'h' in data_type:
        filenames.append(heat_files)

    if 'b' in data_type:
        for entry in boundary_files:
            filenames.append(entry) 

    if 'c' in data_type:
        filenames.append(core_files) 
        

#%% LOGIN 

def login(username, password, hostip):
    try:
        global ssh_client  # If you need this to be global; consider alternatives
        
        ssh_client = paramiko.SSHClient()
        ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        # Attempt to connect
        ssh_client.connect(hostname=hostip, username=username, password=password)
        print("Connection successful.")
    
    except paramiko.SSHException as e:
        print(f"SSH connection error: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")
    
    return ssh_client

#%% FORMAT FOLDERNAMES

def manyfolders(chosenfolders):
    
    chosenfolders_formatted = chosenfolders.replace(' ', '')
    chosenfolders_split     = chosenfolders_formatted.split(',')
    
    chosenfolders_array     = [int(numeric_string) for numeric_string in chosenfolders_split]
    
    return chosenfolders_array
    
#%% UPLOAD AND TRANSFORM DATA

def uploadtransform(path, file):
    
    if file == boundary_files:
        
        vel  = np.ndarray.tolist(np.loadtxt(path + boundary_files[0], dtype=float))
        # for line in range(len(vel)):
        #     for column in range(len(vel[line])):
        #         vel[line][column] = float(vel[line][column])
        
        visc = np.ndarray.tolist(np.loadtxt(path + boundary_files[1], dtype=float))
        # for line in range(len(visc)):
        #     for column in range(len(visc[line])):
        #         visc[line][column] = float(visc[line][column])
        
        temp = np.ndarray.tolist(np.loadtxt(path + boundary_files[2], dtype=float))
        # for line in range(len(temp)):
        #     for column in range(len(temp[line])):
        #         temp[line][column] = float(temp[line][column])
        
        try:
            shf = np.ndarray.tolist(np.loadtxt(path + boundary_files[3], dtype=str))
            for line in range(len(shf)):
                for column in range(len(shf[line])):
                    
                    E_there     = shf[line][column].find('E+')
                    plus_there  = shf[line][column].find('+')
                    
                    if plus_there != -1 and E_there == -1:
                        shf[line][column] = float(shf[line][column].replace('+', 'E+'))
                    
                    else:
                        shf[line][column] = float(shf[line][column])
        except:
            shf = []
                
                
        data = [vel, visc, temp, shf]
    
    elif file == core_files:

        try:
            data = np.ndarray.tolist(np.loadtxt(path + file, dtype=str, skiprows=1))
            
            for line in range(len(data)):
                for column in range(len(data[line])):
                    
                    E_there     = data[line][column].find('E+')
                    plus_there  = data[line][column].find('+')
                    
                    if plus_there != -1 and E_there == -1:
                        data[line][column] = float(data[line][column].replace('+', 'E+'))
                    
                    else:
                        data[line][column] = float(data[line][column])
        except:
            data = []
            print(file, 'could not be found in this folder.')
            print('didnt work')
    
    else:
        
        try:
            
            data = np.ndarray.tolist(np.loadtxt(path + file, dtype=str))
            
    
            for line in range(len(data)):
                for column in range(len(data[line])):
                    
                    E_there     = data[line][column].find('E+')
                    plus_there  = data[line][column].find('+')
                    
                    if plus_there != -1 and E_there == -1:
                        data[line][column] = float(data[line][column].replace('+', 'E+'))
                    
                    else:
                        data[line][column] = float(data[line][column])
        except:
            data = []
            print(file, 'could not be found in this folder.')
                
    return data

#%% ASSIGN NAMES TO DATA

def assign(data, collected_data, folder, file, params_units):
    
    if file == input_files:
        
        g                   = [item[0] for item in data]
        P                   = [item[1] for item in data]
        rho                 = [item[2] for item in data]
        r                   = [item[3] for item in data]
        T                   = [item[4] for item in data]
        T_melt              = [item[5] for item in data]
        C                   = [item[6] for item in data]
        alpha               = [item[7] for item in data]
        Gamma               = [item[8] for item in data]
        K_iso               = [item[9] for item in data]
        K_adi               = [item[10] for item in data]
        mu                  = [item[11] for item in data]
        sigma               = [item[12] for item in data]
        phasenumber         = [item[13] for item in data]
        
        collected_data[folder][list(params_units['i'].keys())[0]]      = g
        collected_data[folder][list(params_units['i'].keys())[1]]      = P
        collected_data[folder][list(params_units['i'].keys())[2]]      = rho
        collected_data[folder][list(params_units['i'].keys())[3]]      = r
        collected_data[folder][list(params_units['i'].keys())[4]]      = T
        collected_data[folder][list(params_units['i'].keys())[5]]      = T_melt
        collected_data[folder][list(params_units['i'].keys())[6]]      = C
        collected_data[folder][list(params_units['i'].keys())[7]]      = alpha
        collected_data[folder][list(params_units['i'].keys())[8]]      = Gamma
        collected_data[folder][list(params_units['i'].keys())[9]]      = K_iso
        collected_data[folder][list(params_units['i'].keys())[10]]      = K_adi
        collected_data[folder][list(params_units['i'].keys())[11]]     = mu
        collected_data[folder][list(params_units['i'].keys())[12]]     = sigma
        collected_data[folder][list(params_units['i'].keys())[13]]     = phasenumber


    if file == heat_files:
        
        time                = [item[0] for item in data]
        heatsource          = [item[1] for item in data]
        cmb_temp_mean       = [item[2] for item in data]
        cmb_temp            = [item[3] for item in data]
        heatflow_surface    = [item[4] for item in data]
        heatflow_core       = [item[5] for item in data]
        
        collected_data[folder][list(params_units['h'].keys())[0]]      = time
        collected_data[folder][list(params_units['h'].keys())[1]]      = heatsource
        collected_data[folder][list(params_units['h'].keys())[2]]      = cmb_temp_mean
        collected_data[folder][list(params_units['h'].keys())[3]]      = cmb_temp
        collected_data[folder][list(params_units['h'].keys())[4]]      = heatflow_surface
        collected_data[folder][list(params_units['h'].keys())[5]]      = heatflow_core
        
    if file == boundary_files:
        
        vel                 = data[0]
        visc                = data[1]
        temp                = data[2]
        shf                 = data[3]
        
        collected_data[folder][list(params_units['b'].keys())[0]]      = vel
        collected_data[folder][list(params_units['b'].keys())[1]]      = visc
        collected_data[folder][list(params_units['b'].keys())[2]]      = temp
        collected_data[folder][list(params_units['b'].keys())[3]]      = shf

    if file == core_files:

        time                = [item[0] for item in data]
        icb_radius          = [item[1] for item in data]
        cmb_temp            = [item[2] for item in data]
        cmb_heat            = [item[3] for item in data]
        sec_cooling         = [item[4] for item in data]
        latent_heat         = [item[5] for item in data]
        cmb_temp_change     = [item[6] for item in data]
        icb_radius_change   = [item[7] for item in data]
        thermal_buoyancy    = [item[8] for item in data]
        magnetic_moment     = [item[9] for item in data]
        cmb_field_strength  = [item[10] for item in data]
        surf_field_strength = [item[11] for item in data]
        
        
        collected_data[folder][list(params_units['c'].keys())[0]]      = time
        collected_data[folder][list(params_units['c'].keys())[1]]      = icb_radius
        collected_data[folder][list(params_units['c'].keys())[2]]      = cmb_temp
        collected_data[folder][list(params_units['c'].keys())[3]]      = cmb_heat
        collected_data[folder][list(params_units['c'].keys())[4]]      = sec_cooling
        collected_data[folder][list(params_units['c'].keys())[5]]      = latent_heat
        collected_data[folder][list(params_units['c'].keys())[6]]      = cmb_temp_change
        collected_data[folder][list(params_units['c'].keys())[7]]      = icb_radius_change
        collected_data[folder][list(params_units['c'].keys())[8]]      = thermal_buoyancy
        collected_data[folder][list(params_units['c'].keys())[9]]      = magnetic_moment
        collected_data[folder][list(params_units['c'].keys())[10]]     = cmb_field_strength
        collected_data[folder][list(params_units['c'].keys())[11]]     = surf_field_strength
        
    
    return collected_data
    