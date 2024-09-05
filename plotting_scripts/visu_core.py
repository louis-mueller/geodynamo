import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#Input
#--------------------------------------------------------------------------------
sim_name = 'Debug'              # name it 'Debug' if you want to test the script
file_name = 'data_core.res'
scale_factor = 5
time_steps = 0                 # 0 means all time steps
#titles = ['Time', 'ICB Radius', 'Density at ICB', 'CMB Temperature', 'Bottom TBL mantle Temperature',
#          'CMB Heat', 'Secular Cooling', 'Latent Heat Release', 'Radiogenic Heating', 'Pressure Heating']
ylabels = ['ICB Radius [%]', 'Rayleigh Num.', 'TBL Thickness [m]', 'Temperatures [K]',
          'Core Heat [TW]', r'Magnetic Field Strengnth [$\mu$T]']
clrs = np.array([ [92, 75, 81],[140, 190, 178],[243, 181, 98],[240, 96, 96]]) / 255
#--------------------------------------------------------------------------------

os.chdir('/scratch/Simulations_Louis/')
plots_dir = os.path.join(os.getcwd(), 'Plots')
sim_dir = os.path.join(os.getcwd(), sim_name)
file_dir = os.path.join(sim_dir, file_name)

if not os.path.exists(sim_dir):
    print(f"Error: {sim_dir} does not exist.")
    sys.exit(1)

print(f"Reading data from {file_dir}")

try:
    df = pd.read_csv(file_dir, delim_whitespace=True)
except pd.errors.ParserError as e:
    print(f"Error parsing the file: {e}")
    raise

num_cols = df.shape[1]

# only plot the first n time steps if n > time_steps
if df.shape[0] > time_steps and time_steps > 0:
    df = df.iloc[:time_steps]

ids = df.columns

colums = 3
rows = 2

fig, axs = plt.subplots(rows, colums, figsize=(colums*scale_factor, rows*scale_factor), sharex=True)

axs = axs.flatten()
    
axs[0].plot(df[ids[0]], df[ids[1]], color=clrs[0])
axs[0].set_ylabel(ylabels[0])

axs[1].plot(df[ids[0]], df[ids[2]], color=clrs[0])
axs[1].set_ylabel(ylabels[1])

axs[2].plot(df[ids[0]],df[ids[3]], label=r'$M_{\oplus} = 1$, $w_{Fe} = 30\%$, $n = 5000$', color=clrs[0])
axs[2].set_ylabel(ylabels[2])
axs[2].legend(loc='upper right', fontsize=14)

axs[3].plot(df[ids[0]], df[ids[4]], color=clrs[2], label=r'$T_c$')
axs[3].plot(df[ids[0]], df[ids[4]], color=clrs[1], label=r'$T_b$')
axs[3].plot(df[ids[0]], df[ids[6]], color=clrs[3], label=r'$T_i$')
axs[3].legend(loc='upper right', fontsize=14)
axs[3].set_ylabel(ylabels[3])

axs[4].plot(df[ids[0]], df[ids[7]], color=clrs[2], label=r'$Q_{CMB}$')
axs[4].plot(df[ids[0]], df[ids[8]], color=clrs[1], label=r'$Q_S$')
axs[4].plot(df[ids[0]], df[ids[9]], color=clrs[3], label=r'$Q_L$')
axs[4].plot(df[ids[0]], df[ids[10]], color=clrs[0], label=r'$Q_D$')
axs[4].legend(loc='upper right', fontsize=14)
axs[4].set_ylabel(ylabels[4])

axs[5].plot(df[ids[0]], df[ids[12]]*1e6, color='red', label=r'$B_{CMB}$')
axs[5].plot(df[ids[0]], df[ids[13]]*1e6, color='blue', label=r'$B_{surface}$')
axs[5].legend(loc='upper right', fontsize=14)
axs[5].set_ylabel(ylabels[5])

for ax in axs[-3:]:
    ax.set_xlabel('t [Myr]')

plt.tight_layout()

# Ensure the "Plots" directory exists and save the SVG file there
if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

output_file_name = sim_name+'.svg'

output_file_path = os.path.join(plots_dir, output_file_name)
plt.savefig(output_file_path)
print(f'Plot saved as {output_file_name}')