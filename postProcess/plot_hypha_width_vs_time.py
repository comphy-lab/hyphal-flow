"""
# plot_hypha_width_vs_time.py

Plot hypha width statistics from `hypha-def-log`.

## Input Format

`hypha-def-log` is expected to contain:
- `t`: time
- `y_if_max`: maximum detected interface height
"""

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    """
    Load hypha deformation log data and plot maximum width vs time.
    """
    data = np.loadtxt("hypha-def-log", skiprows=1)
    t = data[:, 0]
    y_if_max = data[:, 1]

    idx = np.argsort(t)
    t = t[idx]
    y_if_max = y_if_max[idx]

    width = y_if_max

    plt.figure()
    plt.plot(t, width, label="max hypha width (logged)")
    plt.xlabel("Time")
    plt.ylabel("Max hypha width")
    plt.title("Max hypha width vs time")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
