import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("data/vel_histogramB3.csv")

# Reset index to get usable columns
df_clean = df.reset_index()

# Rename columns for clarity
df_clean.columns = ["v_center", "density", "count"]

# Convert columns to numeric
df_clean["v_center"] = pd.to_numeric(df_clean["v_center"], errors="coerce")
df_clean["density"] = pd.to_numeric(df_clean["density"], errors="coerce")
df_clean["count"] = pd.to_numeric(df_clean["count"], errors="coerce")


# Drop any NaN rows (in case of header junk)
df_clean = df_clean.dropna()

df_clean["count"] = df_clean["count"]/max(df_clean["count"])

print(sum(df_clean["count"]))

# Speed range
v = np.linspace(0, 10, 100)  # speeds in m/s
m = 1
k_B = 1
T = 1

# Maxwell-Boltzmann distribution
f_v = ( (m / (2 * np.pi * k_B * T))**1.5 ) * 4 * np.pi * v**2 * np.exp(-m * v**2 / (2 * k_B * T))

f_v = f_v/max(f_v)

# print(max(f_v))


plt.figure(figsize=(10,6))
# plt.bar(df_clean["v_center"], df_clean["count"], width=0.4, align='center')
plt.plot(df_clean["v_center"], df_clean["count"], 'o')
plt.plot(v, f_v,'-')

plt.xlim([0, 10])

plt.xlabel("Velocity Center (v_center)")
plt.ylabel("Count")
plt.title("Histogram of Velocity Distribution")
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()