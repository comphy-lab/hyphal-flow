"""
# plot_vcm_vs_time.py

Plot droplet center-of-mass velocity from a standard `log` file.

## Dependencies

- `numpy`: tabular data loading.
- `matplotlib`: time-series plotting.

## Input Format

The script expects a text file named `log` with columns:
`i dt t ke vcm`.

## Output

Displays a `vcm(t)` plot in an interactive window.

#### Example

```bash
python3 postProcess/plot_vcm_vs_time.py
```
"""

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    """
    Load time-series data from `log` and generate a velocity-vs-time plot.

    #### Raises

    - `OSError`: The input file cannot be opened.
    - `ValueError`: The file cannot be parsed as a numeric table.
    """
    data = np.loadtxt("log", skiprows=1)
    t = data[:, 2]
    vcm = data[:, 4]

    plt.figure()
    plt.plot(t, vcm)
    plt.xlabel("Time")
    plt.ylabel(r"$v_{\mathrm{cm}}$")
    plt.title("Droplet center-of-mass velocity vs time")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
