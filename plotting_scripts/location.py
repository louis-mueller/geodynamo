#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:10:09 2024

@author: theresabuettner
"""

import paramiko
import numpy as np
from getpass import getpass
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl    

#%% GET FILES

#def getfiles(ssh_homedir, ssh_folder, local_homedir, local_folder, variables):
def getfiles(local_dir, server_dir):
    
    global ssh_path
    global local_path
    
    local_path  = local_dir
    ssh_path    = server_dir
    
    sftp_client = ssh_client.open_sftp()
    sftp_client.chdir(ssh_path)
    
    Path(local_path).mkdir(exist_ok=True)
    
    ssh_filename = []
    
    if 'T' in variables:
        
        ssh_filename.append('data_temp_profs.res')
        
    if 'V' in variables:
        
        ssh_filename.append('data_visc_profs.res')
        
    if 'v' in variables:
        
        ssh_filename.append('data_vel_profs.res')    
    
    if 's' in variables:
        
        ssh_filename.append('data_surf_heatflux_profs.res')    

    local_filename = ssh_filename.copy()
    
    for i in range(len(ssh_filename)):
        
        local_filename[i] = local_filename[i].replace('.res', '.txt')
        sftp_client.get(ssh_path + ssh_filename[i], local_path + local_filename[i])

    sftp_client.close()
    
    

