import os
import pandas as pd
import matplotlib.pyplot as plt

file_name = 'data_core.res'
output_file_name = 'core_plot.svg'
plots_dir = os.path.join(os.getcwd(), 'Plots')
debug_dir = os.path.join(os.getcwd(), 'Debug')
file_dir = os.path.join(debug_dir, file_name)

if not os.path.exists(debug_dir):
    print(f"Error: {debug_dir} does not exist.")
    sys.exit(1)

print(f"Reading data from {file_dir}")

titles = ['Time', 'ICB Radius', 'Density at ICB', 'CMB Temperature', 'Bottom TBL mantle Temperature',
          'CMB Heat', 'Secular Cooling', 'Latent Heat Release', 'Radiogenic Heating', 'Pressure Heating']

try:
    df = pd.read_csv(file_dir, delim_whitespace=True)
except pd.errors.ParserError as e:
    print(f"Error parsing the file: {e}")
    raise

if df.shape[1] != 10:
    raise ValueError("The input file must have exactly 10 columns.")

ids = df.columns

fig, axs = plt.subplots(3, 3, figsize=(15, 15), sharex=True)

axs = axs.flatten()

for i in range(1, 10): 
    axs[i-1].plot(df[ids[0]], df[ids[i]], label=ids[i], color='black')
    axs[i-1].set_ylabel(ids[i])
    #axs[i-1].set_title(titles[i])
    #axs[i-1].legend(loc='upper right')
    #axs[i-1].grid(True)

for ax in axs[-3:]:
    ax.set_xlabel('t [Myr]')

plt.tight_layout()

# Ensure the "Plots" directory exists and save the SVG file there

if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

output_file_path = os.path.join(plots_dir, output_file_name)
plt.savefig(output_file_path)
print(f'Plot saved as {output_file_name}')
