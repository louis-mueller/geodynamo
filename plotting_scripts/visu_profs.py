# PyVis-CHIC_profs_res.py
# -*- coding: utf-8 -*-

'''
Python Visualization of Interior Structure Simulations data of profs.res 
Created by Louis Müller (18.04.2024)

Version (12.08.2024)

Summary:
Run This program in the base directory as >>python3 visu_profs.py <DIR_STRS> <OUTPUT_FILE_NAME> <MAX_COUNT>
Any MAX_COUNT of directories named as listed in DIR_STRS containing a file named INT_STRUCT_DATA (e.g., profs.res), 
with the data stored in columns, are plotted in a pre defined configuration of Subplots. 
The color is chosen by a letter snippet defining the planet mass (M), and three further snippets can 
be defined in ID to differentiate simulations accordingly.
'''

import os
import numpy as np 
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.interpolate as interp
from scipy.integrate import quad

# Input
# -----------------------------------------------------------------------------
INT_STRUCT_DATA = 'profs.res'               # This must have a very specific structure othewise the code would not read it correctly
TEST = 0                                    # in configure_simulation_dictionary() the plot attributes are set easier if TEST == TRUE 
ADD_PREM = 1                                # show PREM500.csv data
n = 1000                                    # data size for PREM500 interpolation (n must equal len(data[:,0]))
ID = ['_p', 'Fe30', 'sFe','M']              # Identifiers in dir-name to distinguish simulations
DATA_TO_LOAD = list(range(10))              # Specify the colums of INT_STRUCT_DATA you would like to load
ROWS, COLUMS = 3, 3                         # Setup how the INT_STRUCT_DATA should be plotted (r x c >= len(DATA_TO_LOAD))
SCALE_FACTOR = 4                            # Scalefactor for the figure size (r*s_f,c*s_f)
MARKER_STYLES = ['solid', 'dashed']         # considers perplex or not ID[0]
OPACITY = [0.5, 1]                          # considers changes in the Iron value ID[1]
GRID = False                                # Grid on or off                                 
# -----------------------------------------------------------------------------

# Color range (normalized) 
cmap = mpl.colormaps['viridis_r']
COLOR_PALETTE = cmap(np.linspace(0,1,5))

#np.array([
#    [255, 195, 0], [255, 87, 51], [199, 0, 57], [144, 12, 63], [88, 24, 69]]) / 255

# from:https://colordesigner.io/color-scheme-builder#5C4B51-8CBEB2-F2EBBF-F3B562-F06060 (23.07.24)
COLOR_PALETTE_TEST = np.array([ [92, 75, 81],[140, 190, 178],[243, 181, 98],[240, 96, 96]]) / 255

# Plot Titles and Labels
TITLES = [
    'Gravity', 'Pressure', 'Density', 'Temperature', 'Melt Temperature', 'Heat Capacity',
    'Thermal Expansivity', 'Grüneisen Parameter', 'Isothermal Bulk Modulus', 
    'Adiabatic Bulk Modulus', 'Shear Modulus', 'Electrical Conductivity', 'Material Phase Number'
]

Y_AXIS_LABELS = [
    r'$g \ [\mathrm{m/s^2}]$', r'$p \ [\mathrm{GPa}]$', r'$\rho \ [\mathrm{kg/m^3}]$',
    r'$T \ [\mathrm{K}]$', r'$T_m \ [\mathrm{K}]$', r'$C_p \ [\mathrm{J/kgK}]$',
    r'$\alpha \ [10^{-5} \ \mathrm{s^{-1}}]$', r'$\gamma$', r'$K_T \ [\mathrm{GPa}]$',
    r'$K_S \ [\mathrm{GPa}]$', r'$G_S \ [\mathrm{GPa}]$', r'$\sigma \ [\mathrm{S/m}]$', 
    r'$mat$'
]

LABELS = [r'1 M  30% Fe', r'1 M  60% Fe', r'2 M  30% Fe', r'2 M  60% Fe', r'3 M  30% Fe',
    r'3 M  60% Fe', r'4 M  30% Fe', r'4 M  60% Fe', r'5 M  30% Fe', r'5 M  60% Fe'
]

def configure_simulation_dictionary(DIR_STRS, MAX_COUNT):
    '''Configure the simulation dictionary with color, marker, and opacity values.'''
    sim_dict = {}
    COUNT = 0
    for dir in DIR_STRS:
        if COUNT > MAX_COUNT:
            break
        if TEST:
            sim_dict[dir] = [COLOR_PALETTE_TEST[DIR_STRS.index(dir)+ADD_PREM], MARKER_STYLES[1], OPACITY[1]]
        else:
            # Define the color index based on the planet mass found in the directory name
            if dir[0] == 'M':
                color_index = int(dir[1]) - 1
            else:
                color_index = int(dir[dir.index(ID[3])+2]) - 1

            # Determine marker style and opacity based on directory identifiers
            if ID[0] in dir:
                opacity = OPACITY[0] if ID[1] in dir else OPACITY[1]
                sim_dict[dir] = [COLOR_PALETTE[color_index], MARKER_STYLES[0], opacity]
            else:
                opacity = OPACITY[0] if ID[1] in dir else OPACITY[1]
                sim_dict[dir] = [COLOR_PALETTE[color_index], MARKER_STYLES[1], opacity]
        COUNT += 1
    return sim_dict

def load_data(directory):
    '''Load and process data from the specified directory.'''
    try:
        data = np.loadtxt(os.path.join(directory, INT_STRUCT_DATA), usecols=DATA_TO_LOAD)
        radius = data[:, 3] / 1000  # Convert to km
        data = np.delete(data, 3, axis=1)  # Remove the radius column from data
        n = len(radius) # define the length for later possible interpolation
        return radius, data, n
    except Exception as e:
        print(f'Error loading data from {directory}: {e}')
        return None, None, None

def plot_data(axs, sim_dict, DIR_STRS, MAX_COUNT, LABELS=LABELS):
    '''Plot the data from all directories.'''
    COUNT = 0
    for dir_name in DIR_STRS:
        if COUNT > MAX_COUNT:
            break

        radius, data, n = load_data(dir_name)
        
        if radius is None or data is None or n is None:
            print(f'Error loading data from {dir_name}')
            continue
        
            
        if ROWS*COLUMS < len(DATA_TO_LOAD)-1:
            return print('Error: Not enough Subplots to plot all the expected data. Change ROWS, COLUMS or DATA_TO_LOAD input')
            # Plot each subplot

        for i in range(ROWS):
            for j in range(COLUMS):
                k = i * COLUMS + j
                if k < len(DATA_TO_LOAD)-1:
                    axs[i, j].plot(radius, data[:, k], color=sim_dict[dir_name][0], linestyle=sim_dict[dir_name][1],
                                alpha=sim_dict[dir_name][2], label=LABELS[COUNT])
                    axs[i, j].set_xlabel('Radius [km]')
                    axs[i, j].set_ylabel(Y_AXIS_LABELS[k])
                    axs[i, j].set_title(TITLES[k])
                    axs[i, j].grid(GRID)
   
        COUNT += 1
        
def remove_empty_subplots(fig):
    '''Remove empty subplots from a matplotlib figure.'''
    # Get all axes from the figure
    all_axes = fig.get_axes()
    
    # Find and delete empty subplots
    index = 0
    for ax in all_axes:
        index += 1
        if not ax.has_data() and not ax.get_images():  # Check if subplot is empty
            ax.remove()
            print(f'Axis {index} was not used and therefore removed from the figure.')

def add_prem(n,axs,sim_dict):
    '''Prem load [columns left to right: radius(m),density(kg/m^3),
    Vpv(m/s),Vsv(m/s),Q-kappa,Q-mu,Vph(m/s),Vsh(m/s),eta]'''
    
    data_PREM = np.genfromtxt('PREM500.csv', delimiter=',', skip_header=1)
    PREM_r = data_PREM[:,0]
    PREM_rho = data_PREM[:, 1]
    # Define new radius values (for example, n=1000)
    r = np.linspace(min(PREM_r), max(PREM_r), n)
    G = 6.6743e-11

    # Interpolate density to the new radius values
    density_interpolator = interp.interp1d(PREM_r, PREM_rho, kind='linear', fill_value="extrapolate")
    rho = density_interpolator(r)

    # Initialize arrays for mass, gravitational acceleration, and pressure
    M = np.zeros_like(r)
    g = np.zeros_like(r)
    P = np.zeros_like(r)

    # Compute mass, gravitational acceleration, and pressure
    for i in range(1, len(r)):
        # Integrate mass up to radius r[i] using the trapezoidal rule
        M[i] = 4 * np.pi * np.trapz(rho[:i+1] * r[:i+1]**2, r[:i+1])
        
        # Compute gravitational acceleration at radius r[i]
        g[i] = G * M[i] / r[i]**2
        
        # Compute pressure at radius r[i] using hydrostatic equilibrium
        P[i] = P[i-1] + np.trapz(-rho[i-1:i+1] * g[i-1:i+1], r[i-1:i+1])

    P = (P-P[-1])/1e9 # [GPa]
    axs[0,0].plot(r/1000,g,color=COLOR_PALETTE_TEST[0],linestyle='dashed',alpha=0.5)
    axs[ROWS-1,2].plot(0,0,color=COLOR_PALETTE_TEST[0],linestyle='dashed',alpha=0.5,label='PREM')
    axs[ROWS-1,2].legend(fontsize='small')

    if COLUMS >= 2:
        axs[0,1].plot(r/1000,P,color=COLOR_PALETTE_TEST[0],linestyle='dashed',alpha=0.5)
        if COLUMS >= 3: 
            axs[0,2].plot(r/1000,rho,color=COLOR_PALETTE_TEST[0],linestyle='dashed',alpha=0.5)
         
def main(DIR_STRS, OUTPUT_FILE_NAME, MAX_COUNT):
    '''Main function to configure and plot data.'''

    sim_dict = configure_simulation_dictionary(DIR_STRS, MAX_COUNT)
    fig, axs = plt.subplots(nrows=ROWS, ncols=COLUMS, figsize=(COLUMS*SCALE_FACTOR, ROWS*SCALE_FACTOR))
    plot_data(axs, sim_dict, DIR_STRS, MAX_COUNT)
    if ADD_PREM:
        add_prem(n,axs,sim_dict)
    remove_empty_subplots(fig)
    plt.tight_layout()
    plt.legend(fontsize='small')
    
    # Ensure the "Plots" directory exists and save the SVG file there
    PLOTS_DIR = os.path.join(os.getcwd(), 'Plots')
    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    OUTPUT_FILE_PATH = os.path.join(PLOTS_DIR, OUTPUT_FILE_NAME)
    plt.savefig(OUTPUT_FILE_PATH)
    print(f'Plot saved as {OUTPUT_FILE_NAME}')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage should be: python3 visu_profs.py <DIR_STRS> <OUTPUT_FILE_NAME> <MAX_COUNT>")
        sys.exit(1)

    DIR_STRS = sys.argv[1]
    DIR_STRS = DIR_STRS.split(',')
    OUTPUT_FILE_NAME = sys.argv[2]
    MAX_COUNT = int(sys.argv[3])
    
    main(DIR_STRS, OUTPUT_FILE_NAME, MAX_COUNT)