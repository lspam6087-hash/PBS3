import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the uploaded file to inspect its structure
file_path = "data/grf_data.dat"

# Try reading as whitespace-separated values
try:
    df = pd.read_csv(file_path, delim_whitespace=True, header=None)
except Exception as e:
    df = str(e)

# Extract x and y from the dataframe
x = df[0].values
y = df[1].values

# Estimate bin edges from midpoints
bin_edges = np.concatenate((
    [x[0] - (x[1] - x[0]) / 2],   # first edge
    (x[:-1] + x[1:]) / 2,         # midpoints
    [x[-1] + (x[-1] - x[-2]) / 2] # last edge
))

# Plot histogram-style bar chart
plt.figure(figsize=(8,5))
plt.bar(x, y, width=np.diff(bin_edges), align='center', edgecolor='black')

plt.xlabel("Distance from center [Angstrom]")
plt.ylabel("Probability [-]")
# plt.title("Probability distribution")
plt.savefig("Figures/B2.png")
plt.show()
