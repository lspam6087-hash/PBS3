import matplotlib.pyplot as plt
import pandas as pd


df = pd.read_csv("data/diagnostics.csv")

plt.figure(figsize=(10,8))

# 1) Energies
plt.plot(df["time"], label="Kinetic Energy", color="red")
plt.plot(df["kinetic_energy"], label="Potential Energy", color="blue")
plt.plot(df["potential_energy"], label="Total Energy", color="green")

plt.ylabel("Energy (kJ/mol)")
plt.xlabel("Time (ps)") 
plt.title("Energies vs Time")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
 