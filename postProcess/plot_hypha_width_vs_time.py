import numpy as np
import matplotlib.pyplot as plt

# File written by event log_hypha_deformation:
# columns: t  y_top_max
data = np.loadtxt("hypha-def-log", skiprows=1)

t = data[:, 0]
y_top_max = data[:, 1]

# Sort by time (safe)
idx = np.argsort(t)
t = t[idx]
y_top_max = y_top_max[idx]

# "Max hypha width" (as logged) = y_top_max
# If you want DIAMETER assuming symmetry about y=0, use: width = 2*y_top_max
width = y_top_max

# Optional: deformation relative to initial value
w0 = width[0]
deformation = width - w0

# Plot width vs time
plt.figure()
plt.plot(t, width, label="max hypha width (logged)")

# Uncomment if you also want deformation on the same plot
# plt.plot(t, deformation, "--", label="deformation = width - width(t=0)")

plt.xlabel("Time")
plt.ylabel("Max hypha width")
plt.title("Max hypha width vs time")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Optional save:
# plt.savefig("hypha_width_vs_time.png", dpi=300)
