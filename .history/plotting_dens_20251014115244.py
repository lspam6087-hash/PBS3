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
plt.plot(data["x_center"], data["density_A"], label='Density A')
plt.plot(data["x_center"], data["density_B"], label='Density B')
plt.plot(data["x_center"], data["density_total"], label='Density Total')

plt.xlabel("X Center")
# plt.ylim([0,3])
plt.ylabel("Density")
plt.title("Density Distribution (A, B, Total)")
plt.legend()
plt.grid(True)
plt.tight_layout()
# plt.savefig("Figures/D2_N2.png")
plt.show()