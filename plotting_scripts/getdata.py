#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 16:02:01 2024

@author: theresabuettner
"""

from input_data import *
from supp_functions import *

for name in safety:
    del globals()[name]

from glob import glob
from stat import S_ISDIR
import os.path
import paramiko
import numpy as np
from getpass import getpass
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl


#%% DEFINE FINDFOLDERS

def findfolders(loc, local_dir, server_dir):
    
    if loc == 'l':
        homedir = local_dir
        
    else:
        homedir = server_dir
    
    testfolders = []

    if loc == 'r':
        
        for name in safety:
            from input_data import name
            
        login(username, password, hostip)
        
        for name in safety:
            del globals()[name]
        
        global sftp_client
        sftp_client = ssh_client.open_sftp()
        sftp_client.chdir(homedir)
    
        for entry in sftp_client.listdir_attr(homedir):
            
            mode = entry.st_mode
            
            if S_ISDIR(mode):
                testfolders.append(entry.filename)
        
    else:
        
        allfolders = glob(homedir+'*/', recursive = True)
        
        for entry in allfolders:
            entry = entry.replace(homedir, '').replace('/','')
            testfolders.append(entry)


    folderlist = {}
    foldernumbers = list(range(1, len(testfolders)+1))
    
    for number in foldernumbers:
        for folder in testfolders:
            
            folderlist[number] = folder
            testfolders.remove(folder)
            break
    
    print('\n')
    for x,y in folderlist.items():
        print(x,y)
    
    global chosenfolders
    chosenfolders = input('\nHere are the folders you can choose from. \n'\
                          'Please enter the number(s) corresponding to your chosen folders, separated by a comma. ' \
                          'If you want to choose all folders, enter 0. \n')
        
    chosenfolders_array = manyfolders(chosenfolders)
        
    global folderstrings
    folderstrings = []
        
    for entry in chosenfolders_array:
        
        value = folderlist.get(entry)
        
        if value is not None:
            folderstrings.append(folderlist.get(entry))
            
        else:
            pass

#%% GET FILES

def getfiles(comparison, local_dir, server_dir):

    # LOGIN
    for name in safety:
        from input_data import name
        
    login(username, password, hostip)
    
    for name in safety:
        del globals()[name]

    # ACCESS REMOTE FILES
    global sftp_client
    sftp_client = ssh_client.open_sftp()
    sftp_client.chdir(server_dir)
        
    if comparison == True:
        
        sftp_client = ssh_client.open_sftp()
        sftp_client.chdir(server_dir)
        
        print('\nThe folders you have chosen are: ')
        
        for folder in folderstrings:
            
            print(folder)
            
            Path(local_dir + folder).mkdir(exist_ok=True)
            
            for filename in filenames:

                local_path = local_dir + folder + filename
                try:
                    sftp_client.get(folder + filename, local_path)
                except:
                    print(filename, 'does not exist in folder', folder)
                
    else:
        
        for filename in filenames:
            Path(local_dir).mkdir(exist_ok=True)
            
            try:
                sftp_client.get(server_dir + filename, local_dir + filename)  
            except:
                print(filename, 'does not exist in this folder.')
     
    sftp_client.close()
    ssh_client.close()
    
    print('\nThe desired files have been downloaded. \n\n')


#%% IMPORT FILES

def importfiles(data_type, homedir, comparison, params_units):

    collected_data = dict()    

    if comparison == True:

        for folder in folderstrings:
            
            collected_data[folder] = {}
            
            path  = homedir + folder
            #files = [path + i for i in filenames]
            
            if 'i' in data_type:
                
                profs = uploadtransform(path, input_files)
                profs_assigned = assign(profs, collected_data, folder, input_files, params_units)
                
            if 'h' in data_type:
                
                data_char_val = uploadtransform(path, heat_files)
                data_char_val_assigned = assign(data_char_val, collected_data, folder, heat_files, params_units)
    
            
            if 'b' in data_type:
                
                data_profs = uploadtransform(path, boundary_files)
                data_profs_assigned = assign(data_profs, collected_data, folder, boundary_files, params_units)
            
            if 'c' in data_type:
                
                data_core = uploadtransform(path, core_files)
                data_core_assigned = assign(data_core, collected_data, folder, core_files, params_units)
                
    else:
        
        path = homedir[:-1]
        folder = path[path.rfind('/')+1:]
        
        collected_data[folder] = {}
        
        if 'i' in data_type:
            
            profs = uploadtransform(path, input_files)
            profs_assigned = assign(profs, collected_data, folder, input_files, params_units)
            
        if 'h' in data_type:
            
            data_char_val = uploadtransform(path, heat_files)
            data_char_val_assigned = assign(data_char_val, collected_data, folder, heat_files, params_units)

        
        if 'b' in data_type:
            
            data_profs = uploadtransform(path, boundary_files)
            data_profs_assigned = assign(data_profs, collected_data, folder, boundary_files, params_units)
        
        if 'c' in data_type:

            data_core = uploadtransform(path, core_files)
            data_core_assigned = assign(data_core, collected_data, folder, core_files, params_units)


    return collected_data









        
        
        
        
        
        
        
        
        
        
        