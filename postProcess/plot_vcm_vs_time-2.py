import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------
# Read log file, skipping non-numeric header lines
# --------------------------------------------------
t = []
vcm = []

with open("log", "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        # skip header / description lines
        if line.startswith("Level") or line.startswith("i"):
            continue

        parts = line.split()
        if len(parts) < 5:
            continue

        # columns: i dt t ke vcm
        t.append(float(parts[2]))
        vcm.append(float(parts[4]))

t = np.array(t)
vcm = np.array(vcm)

# --------------------------------------------------
# Plot
# --------------------------------------------------
plt.figure()
plt.plot(t, vcm, marker="o")
plt.xlabel("Time")
plt.ylabel(r"$v_{\mathrm{cm}}$")
plt.title("Centre-of-mass velocity vs time")
plt.grid(True)
plt.tight_layout()
plt.show()

# Optional save
# plt.savefig("vcm_vs_time.png", dpi=300)
