"""
# plot_hypha_width_vs_time.py

Plot hypha width statistics from `hypha-def-log`.

## Dependencies

- `numpy`: load and sort tabular data.
- `matplotlib`: generate the time-series figure.

## Input Format

`hypha-def-log` is expected to contain:
- `t`: time
- `y_if_max`: maximum detected interface height

## Output

Displays an interactive plot of `y_if_max` against time.

#### Example

```bash
python3 postProcess/plot_hypha_width_vs_time.py
```
"""

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    """
    Load `hypha-def-log` and plot maximum hypha width versus time.

    #### Raises

    - `OSError`: The input file cannot be opened.
    - `ValueError`: The file does not contain at least two numeric columns.
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
