import numpy as np
import matplotlib.pyplot as plt

# ---------------------------
# Load data
# ---------------------------
# Columns: i, dt, t, ke, vcm
data = np.loadtxt("log", skiprows=1)

t   = data[:, 2]   # time
vcm = data[:, 4]   # centre-of-mass velocity (U_d)

# ---------------------------
# Plot
# ---------------------------
plt.figure()
plt.plot(t, vcm)
plt.xlabel("Time")
plt.ylabel(r"$v_{\mathrm{cm}}$")
plt.title("Droplet centre-of-mass velocity vs time")
plt.grid(True)

plt.tight_layout()
plt.show()
