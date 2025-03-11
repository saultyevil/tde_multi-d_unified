import os

import numpy as np
import pysi.util.array
from astropy.constants import c
from matplotlib import pyplot as plt
from pysi.wind import Wind
from pysi.wind.enum import DistanceUnits


def find_closest_element_index(arr: np.ndarray, target: float) -> int:
    """
    Finds the index of the element in the Numpy array with a value closest to
    the target value.

    """
    return int(np.argmin(np.abs(arr - target)))


def plot_multi_d_cell_seds(
    wind1d: Wind, wind2d: Wind, display: bool
) -> tuple[plt.Figure, plt.Axes]:
    """Create a plot comparing cell SEDs from a 1D and the 2D model."""
    fig, ax = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    ax = ax.flatten()

    plot_alpha = 0.75
    smooth = 3
    target_angle = 14
    theta_cell = 128 / 90 * target_angle

    original_units1d = wind1d.distance_units
    original_units2d = wind2d.distance_units

    wind1d.change_units(DistanceUnits.GRAVITATIONAL_RADIUS)
    wind2d.change_units(DistanceUnits.GRAVITATIONAL_RADIUS)

    for m, r_cell in enumerate([2, 60, 110, 155]):
        r = wind2d["r"][int(r_cell), int(theta_cell)]
        rr_cell = find_closest_element_index(wind1d["r"], r)
        if wind1d["inwind"][rr_cell] < 0:
            rr_cell += 1
            if rr_cell >= wind1d.nx:
                rr_cell = wind1d.nx - 1
        r = wind1d["r"][rr_cell]

        # 1d simulation
        cell_spec_flux = wind1d["spec_flux"][int(rr_cell)]
        cell_spec_freq = wind1d["spec_freq"][int(rr_cell)]

        wavelength = (c.value / cell_spec_freq) / 1e-10
        ax[m].plot(
            wavelength,
            pysi.util.array.smooth_array(cell_spec_freq * cell_spec_flux, smooth),
            alpha=plot_alpha,
            label="1D",
        )

        # 2d simulation
        cell_spec_flux = wind2d["spec_flux"][int(r_cell), int(theta_cell)]
        cell_spec_freq = wind2d["spec_freq"][int(r_cell), int(theta_cell)]

        r = wind2d["r"][int(r_cell), int(theta_cell)]
        wavelength = (c.value / cell_spec_freq) / 1e-10
        ax[m].plot(
            wavelength,
            pysi.util.array.smooth_array(cell_spec_freq * cell_spec_flux, smooth),
            label="2D",
            alpha=plot_alpha,
            linestyle="-",
        )

        # set axes styles and labels
        ax[m].set_xlim(c.value / 2e17 / 1e-10, c.value / 1e14 / 1e-10)
        ax[m].set_ylim(1e10, 5e17)
        ax[m].text(
            0.03,
            0.91,
            r"$\log_{10}(r / r_{g}) = " + f"{np.log10(r):.1f}$",
            transform=ax[m].transAxes,
        )
        ax[m] = pysi.util.plot.set_axes_scales(ax[m], "loglog")

    leg0 = ax[1].legend(loc="upper right")
    leg0.set_zorder(0)

    fig.text(
        0.02,
        0.5,
        r"$\lambda J_{\lambda}$ [ergs s$^{-1}$ cm$^{-3}$ sr$^{-1}$]",
        rotation="vertical",
        ha="center",
        va="center",
    )
    fig.text(0.5, 0.02, "Rest-frame wavelength [\AA]", ha="center", va="center")
    fig = pysi.util.plot.finish_figure(fig, hspace=0, wspace=0)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.jpeg", dpi=300)

    wind1d.change_units(original_units1d)
    wind2d.change_units(original_units2d)

    if display:
        plt.show()
    else:
        plt.close()

    return fig, ax
