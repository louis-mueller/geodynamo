import pandas as pd
import matplotlib.pyplot as plt

# Define the file path using os.path.join for platform independence
file_path = 'data_core.res'

# Load the .res file into a DataFrame with the first line as header
try:
    df = pd.read_csv(file_path, delim_whitespace=True)
except pd.errors.ParserError as e:
    print(f"Error parsing the file: {e}")
    # Handle the error, perhaps by inspecting or fixing the file
    raise

# Ensure the DataFrame has the correct number of columns
if df.shape[1] != 10:
    raise ValueError("The input file must have exactly 10 columns.")

# Get the column names from the header
column_names = df.columns

# Plotting
fig, axs = plt.subplots(3, 3, figsize=(15, 15), sharex=True)  # 3x3 grid of subplots

# Flatten the 3x3 array of axes to easily loop through them
axs = axs.flatten()

for i in range(1, 10):  # Loop over columns 2 to 10 (index 1 to 9)
    axs[i-1].plot(df[column_names[0]], df[column_names[i]], label=column_names[i])
    axs[i-1].set_ylabel(column_names[i])
    axs[i-1].legend(loc='upper right')
    axs[i-1].grid(True)

# Set the x-axis label for the last row of subplots
for ax in axs[-3:]:
    ax.set_xlabel(column_names[0])

plt.tight_layout()
plt.show()
