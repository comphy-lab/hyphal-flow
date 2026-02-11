"""
# plot_vcm_vs_time-2.py

Robust `log` parser for center-of-mass velocity time series.

Unlike the minimal parser, this script skips descriptive header lines
before plotting `vcm` against time.

## Dependencies

- `numpy`: array conversion for parsed columns.
- `matplotlib`: time-series visualization.

#### Example

```bash
python3 postProcess/plot_vcm_vs_time-2.py
```
"""

import matplotlib.pyplot as plt
import numpy as np


def read_log(path: str = "log") -> tuple[np.ndarray, np.ndarray]:
    """
    Parse a simulation `log` file and return `(t, vcm)` arrays.

    #### Args

    - `path`: Log file path (default: `log`).

    #### Returns

    - `tuple[np.ndarray, np.ndarray]`: Time and center-of-mass velocity arrays.

    #### Raises

    - `OSError`: The log file cannot be opened.
    - `ValueError`: A parsed row contains non-numeric values.
    """
    t_values = []
    vcm_values = []

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("Level") or line.startswith("i"):
                continue

            parts = line.split()
            if len(parts) < 5:
                continue

            t_values.append(float(parts[2]))
            vcm_values.append(float(parts[4]))

    return np.array(t_values), np.array(vcm_values)


def main() -> None:
    """
    Load data with header filtering and plot droplet velocity vs time.

    #### Raises

    - `OSError`: The input log file cannot be opened.
    - `ValueError`: Parsed data is malformed.
    """
    t, vcm = read_log("log")

    plt.figure()
    plt.plot(t, vcm, marker="o")
    plt.xlabel("Time")
    plt.ylabel(r"$v_{\mathrm{cm}}$")
    plt.title("Center-of-mass velocity vs time")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
