import os

import numpy as np
import pypython
from matplotlib import pyplot as plt


def find_closest_element_index(arr: np.ndarray, target: float) -> int:
    """
    Finds the index of the element in the NumPy array with a value closest to the target value.

    Parameters
    ----------
    arr : numpy.ndarray
        A NumPy array of numeric values.
    target : float
        The target value to find the closest element to.

    Returns
    -------
    int
        The index of the element in the array closest to the target value.
    """
    return int(np.argmin(np.abs(arr - target)))


PLOT_ALPHA = 0.75
SMOOTH = 3
wind1d = pypython.Wind(
    "input", "../simulations/1d/Bin-4", distance_units="rg", version="88"
)
wind2d = pypython.Wind(
    "input", "../simulations/2d/reduced", distance_units="rg", version="88"
)

fig, ax = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
ax = ax.flatten()

TARGET_ANGLE = 14
THETA_CELL = 128 / 90 * TARGET_ANGLE

for m, r_cell in enumerate([2, 60, 110, 155]):
    r = wind2d.r[int(r_cell), int(THETA_CELL)]
    rr_cell = find_closest_element_index(wind1d.r, r)
    cell_spec = wind1d.cell_spec[rr_cell]
    if cell_spec is None:
        rr_cell += 1
        if rr_cell >= wind1d.nx:
            rr_cell = wind1d.nx - 1
        cell_spec = wind1d.cell_spec[rr_cell]
    if cell_spec is None:
        continue
    r = wind1d.r[rr_cell, 0]

    # benchmark
    ax[m].plot(
        cell_spec["Freq."],
        pypython.smooth_array(cell_spec["Freq."] * cell_spec["Flux"], SMOOTH),
        alpha=PLOT_ALPHA,
        label="1D",
    )

    # 2d
    cell_spec = wind2d.cell_spec[int(r_cell), int(THETA_CELL)]
    if cell_spec is None:
        print(f"Skipping theta = {THETA_CELL} and r = {r_cell}")
        continue

    r = wind2d.r[int(r_cell), int(THETA_CELL)]

    ax[m].plot(
        cell_spec["Freq."],
        pypython.smooth_array(cell_spec["Freq."] * cell_spec["Flux"], SMOOTH),
        label="2D",
        alpha=PLOT_ALPHA,
        linestyle="-",
    )
    ax[m].set_xlim(1e14, 2e17)
    ax[m].set_ylim(1e10, 3e17)
    ax[m].text(
        0.03,
        0.91,
        r"$\log_{10}(r / r_{g}) = " + f"{np.log10(r):.1f}$",
        transform=ax[m].transAxes,
    )
    ax[m] = pypython.plot.set_axes_scales(ax[m], "loglog")

leg0 = ax[1].legend(loc="upper right")
leg0.set_zorder(0)

fig.text(
    0.02,
    0.5,
    r"$\nu J_{\nu}$ [ergs s$^{-1}$ cm$^{-3}$ sr$^{-1}$]",
    rotation="vertical",
    ha="center",
    va="center",
)
fig.text(0.5, 0.02, "Rest-frame frequency [Hz]", ha="center", va="center")
fig = pypython.plot.finish_figure(fig, hspace=0, wspace=0)
fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
plt.show()
