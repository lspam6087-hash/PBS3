# import pandas as pd
# import matplotlib.pyplot as plt

# df = pd.read_csv("data\dens_histogramC2.csv")

# # Assign generic column names

# print(df.head())
# # Plot the three series against col_0
# plt.figure(figsize=(8,6))
# plt.plot(df["x_center"], df["density_A"], label='Rho_A', marker='o')
# # plt.plot(data['col_0'], data['col_2'], label='Series 2', marker='s')
# # plt.plot(data['col_0'], data['col_3'], label='Series 3', marker='^')

# plt.xlabel("col_0")
# plt.ylabel("Values")
# plt.title("Line Plot of Data")
# plt.legend()
# plt.grid()
# plt.show()

# It seems there might be hidden characters or whitespace in the column names.
# Let's inspect the column names to confirm and clean them up before plotting.

# Reload file
import pandas as pd
import matplotlib.pyplot as plt
file_path = "data/dens_histogramC2.csv"
data = pd.read_csv(file_path)

# Display column names
print(data.columns.tolist())

# Clean column names (strip whitespace)
data.columns = data.columns.str.strip()

# Plot again


plt.figure(figsize=(8, 5))
plt.plot(data["x_center"], data["density_A"], marker='o', label='Density A')
plt.plot(data["x_center"], data["density_B"], marker='s', label='Density B')
plt.plot(data["x_center"], data["density_total"], marker='^', label='Density Total')

plt.xlabel("X Center")
plt.ylabel("Density")
plt.title("Density Distribution (A, B, Total)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()