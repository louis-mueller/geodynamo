# 1D Interior Structure (profs.res) Data2Lib Function

import pandas as pd
import plotting_input as pf
import numpy as np
import scipy as sp
import os


def prof_data2lib(profs, prem='PREM500.csv', names=pf.names, colors=pf.cmap_inferno, prem_clr=pf.prem_clr,prem_lbl=pf.prem_lbl
                  ,prem_ls=pf.prem_ls):
    
   # Check if profs is a list; if not, convert it to a list
    if not isinstance(profs, list):
        profs = [profs]
    
    lbl, texture, planets, cores = [], {}, {}, {}

    # Load data for each profile
    #---------------------------------------------------------------------------
    for prof in profs:
        if not isinstance(prof, str):
            raise ValueError("ERROR: prof_data2lib: profs must be a list of strings. Found an invalid element: {}".format(prof))
        
        string = prof[0:16]
        string = string.replace('_', ' ')
        string = string.replace('Fe','CMF')
        string = string.replace('-','.')
        lbl.append(string)
        df = pd.read_csv(prof,delimiter=r'\s+', names=names)

        for cmb in range(len(df['phase'])):
            if df['phase'][cmb] == 8.0:
                break

        i = profs.index(prof)
        df['P'] = df['P']*1.0e9
        df['alpha'] = df['alpha']*1.0e-5
        df['r']=df['r']*1.0e-3
        planets[prof] = df

        df = df.drop(df.index[:cmb])
        ncmb = 1000-cmb
        df.index = range(0,ncmb)
        df['n'] = ncmb
        df['r'] = df['r']-df['r'][ncmb-1]
        clrs = colors(np.linspace(0,1,len(profs)))
        texture[prof] = {'label': lbl[i], 'color': clrs[i], 'linestyle': '-'}
        cores[prof] = df
    #---------------------------------------------------------------------------

    # Load PREM data
    #---------------------------------------------------------------------------
    try:
        data_PREM = np.genfromtxt(prem, delimiter=',', skip_header=1)
        PREM_r = data_PREM[:,0]
        PREM_rho = data_PREM[:, 1]
    except:
        print("ERROR: prof_data2lib: could not load PREM data")
    
    # Define new radius values (for example, n=1000)
    n=len(planets[profs[0]]['r'])
    r = np.linspace(min(PREM_r), max(PREM_r), n)
    G = 6.6743e-11

    try:
        # Interpolate density to the new radius values
        density_interpolator = sp.interpolate.interp1d(PREM_r, PREM_rho, kind='linear', fill_value="extrapolate")
        rho = density_interpolator(r)
    except:
        print("ERROR: prof_data2lib: could not interpolate PREM data")

    # Initialize arrays for mass, gravitational acceleration, and pressure
    M = np.zeros_like(r)
    g = np.zeros_like(r)
    P = np.zeros_like(r)

    try:
        # Compute mass, gravitational acceleration, and pressure
        for i in range(1, len(r)):
            # Integrate mass up to radius r[i] using the trapezoidal rule
            M[i] = 4 * np.pi * np.trapz(rho[:i+1] * r[:i+1]**2, r[:i+1])
            
            # Compute gravitational acceleration at radius r[i]
            g[i] = G * M[i] / r[i]**2
            
            # Compute pressure at radius r[i] using hydrostatic equilibrium
            P[i] = P[i-1] + np.trapz(-rho[i-1:i+1] * g[i-1:i+1], r[i-1:i+1])

        P = (P-P[-1])/1e9 # [GPa]

        # find r core where g is maximum
        rc = r[np.argmax(g)]

        # find the density jump at the inner core boundary where rho suddenly increases
        ri = r[np.argmax(-1*np.gradient(rho[:np.argmax(g)]))]
    except:
        print("ERROR: prof_data2lib: could not compute mass, gravitational acceleration, and pressure for PREM data")
    
    try:
        planets['PREM'] = pd.DataFrame({'r': r/1000, 'rho': rho, 'g': g, 'P': P, 'M': M, 'rc': rc/1000, 'ri': ri/1000})
        texture['PREM'] = {'label': prem_lbl, 'color': prem_clr, 'linestyle': prem_ls}
    except:
        print("ERROR: prof_data2lib: could not create DataFrame for PREM data")
    #---------------------------------------------------------------------------
        
    return planets, cores, texture



def core_data(sim_dir,data_name,getcore=False,clm_names=pf.names,colors=pf.cmap_inferno):
    '''
    sim _dir : directory where any number of Simulation directories are stored
    data_name : name of the data file to be loaded 
    '''

    sim_list = os.listdir(sim_dir)
    sim_list = [d for d in sim_list if os.path.isdir(os.path.join(sim_dir, d))]

    for i, d in enumerate(sim_list):
        print(i, d)

    sim_num = input("Input the numbers of Simulation directories to be stored (comma seperatted): ")

    sim_num = sim_num.split(",") # split the input into a list of strings
    sim_num = [int(i) for i in sim_num]
    sim_list = [sim_list[i] for i in sim_num]
    
    clrs = colors(np.linspace(0,1,len(sim_list)//2))

    core_data, texture = {}, {}
    for j, sim in enumerate(sim_list):
        sim_data = os.path.join(sim_dir, sim, data_name)

        # Check if the file has a header
        first_row = pd.read_csv(sim_data, nrows=1, sep=r'\s+', header=None)
        first_row = first_row.values[0]

        if all(isinstance(x, str) for x in first_row):
            # Assume the first row is a header
            df = pd.read_csv(sim_data, sep=r'\s+', header=0)
        else:
            # Assume the first row is not a header
            df = pd.read_csv(sim_data, sep=r'\s+', header=None)
            
            if len(clm_names) == len(df.columns):
                df.columns = clm_names  # Rename the columns directly
            else:
                raise ValueError("Number of new column names does not match number of columns in DataFrame, COlums are numbered as default")
        
        if data_name == "profs.res" and getcore == True:
            for i in range(len(df['phase'])):
                if df['phase'][i] == 8.0:
                    cmb = i
                    # delete all values after index cmb
                    break
            df = df.drop(df.index[:cmb])
            # convert the index from 0 to 511
            ncmb = 1000-cmb
            df.index = range(0,ncmb)

            # corrections
            df['n'] = ncmb
            df['r'] = df['r']-df['r'][ncmb-1]
            df['P'] = df['P']*1.0e9
            df['alpha'] = df['alpha']*1.0e5
    
        if j % 2 == 0:
            a = 0.5
        else:
            a = 1.0
        texture[sim] = {'color': clrs[j//2], 'linestyle': '-', 'alpha': a}
        core_data[sim] = df
        print(sim, "loaded")

    return core_data, texture