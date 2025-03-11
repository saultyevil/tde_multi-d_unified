import os

import numpy as np
import pysi
import pysi.util.array
from astropy.constants import c
from matplotlib import pyplot as plt
from pysi.spec import labels
from pysi.util.array import find_where_target_in_array
from pysi.wind import Wind
from pysi.wind.enum import DistanceUnits


def plot_model_cell_seds(wind: Wind, display: bool) -> tuple[plt.Figure, plt.Axes]:
    """Create a plot showing how the cell SEDs evolve with radius and angle."""
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.5), sharex=True, sharey=True)

    original_units = wind.distance_units
    wind.change_units(DistanceUnits.GRAVITATIONAL_RADIUS)

    theta_target = [14, 77]

    input_bb = np.loadtxt("input_bb_nu_jnu.txt")

    for n, theta_value in enumerate(theta_target):
        theta_cell = find_where_target_in_array(wind["theta"][0, :], theta_value)
        theta_target[n] = theta_cell  # not good!
        rg_target = [112, 274, 797, 1355, 2547, 3600]
        for rg_value in rg_target:
            r_cell = find_where_target_in_array(wind["r"][:, int(theta_cell)], rg_value)
            flux = wind["spec_flux"][int(r_cell), int(theta_cell)]
            freq = wind["spec_freq"][int(r_cell), int(theta_cell)]

            wavelength = (c / freq) / 1e-10

            r = wind["r"][int(r_cell), int(theta_cell)]
            ax[n].plot(
                wavelength,
                pysi.util.array.smooth_array(freq * flux, 4),
                label=f"${np.log10(r):.1f}$",
                alpha=0.75,
            )

        ax[n].set_xlim(20, 20000)
        ax[n].set_ylim(1e10, 1e18)

        ax[n] = pysi.util.plot.set_axes_scales(ax[n], "loglog")

        ax[n].plot(
            (c / input_bb[:, 0]) / 1e-10,
            pysi.util.array.smooth_array(input_bb[:, 1], 100),
            linestyle="--",
            color="k",
            zorder=0,
        )

    leg = ax[1].legend(
        loc="upper right",
        fontsize=10,
        ncol=2,
    )
    ax[0].text(
        0.92,
        0.92,
        r"$" + f"{wind['theta'][int(2), int(theta_target[0])]:.0f}" + r"^{\circ}$",
        transform=ax[0].transAxes,
    )
    ax[1].text(
        0.6,
        0.92,
        r"$" + f"{wind['theta'][int(2), int(theta_target[1])]:.0f}" + r"^{\circ}$",
        transform=ax[1].transAxes,
    )
    leg.set_zorder(0)
    fig.text(
        0.02,
        0.5,
        r"$\lambda J_{\lambda}$ [ergs s$^{-1}$ cm$^{-3}$ sr$^{-1}$]",
        rotation="vertical",
        ha="center",
        va="center",
    )

    edges = [
        [r"O \textsc{vii}", 16],
        [r"C \textsc{v}", 32],
        [r"O \textsc{vi }", 93],
        [r"He \textsc{ii}", 227],
        [r"He \textsc{i}", 504],
    ]

    ax[0] = labels.add_transition_labels_to_ax(
        ax[0], edges, label_linestyle="none", label_offset=0
    )
    ax[1] = labels.add_transition_labels_to_ax(
        ax[1], edges, label_linestyle="none", label_offset=0
    )

    fig.text(0.5, 0.03, r"Rest-frame wavelength [\AA]", ha="center", va="center")
    fig = pysi.util.plot.finish_figure(fig, hspace=0, wspace=0)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.jpeg", dpi=300)

    if display:
        plt.show()
    else:
        plt.close()

    wind.change_units(original_units)

    return fig, ax
