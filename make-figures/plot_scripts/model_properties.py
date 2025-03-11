"""Create a plot of the physical properties of the wind."""

import os

import numpy as np
import pysi
from matplotlib import pyplot as plt
from pysi.wind import Wind
from pysi.wind.enum import DistanceUnits


def plot_model_properties(wind: Wind, show_figure: bool) -> tuple[plt.Figure, plt.Axes]:
    """Create a plot of the physical properties of the wind."""
    original_units = wind.distance_units
    wind.change_units(DistanceUnits.GRAVITATIONAL_RADIUS)
    panel_parameters = ["t_e", "h_density", "t_r", "ip", "H_i01_frac", "He_i02_frac"]
    panel_names = [
        r"Electron temperature",
        r"Hydrogen density",
        r"Radiation temperature",
        r"Ionization parameter",
        r"H~\textsc{i} fraction",
        r"He~\textsc{ii} fraction",
    ]
    panel_units = [" [K]", r" [cm$^{-3}$]", " [K]", "", "", ""]
    fig, ax = plt.subplots(3, 2, figsize=(12, 15), subplot_kw={"projection": "polar"})

    # In each panel, plot sight lines and markers for where cells in seds are
    sight_lines_plotted = []
    outlines_plotted = []
    for inclination, marker, linestyle in zip(
        ["14", "34", "56", "77"], ["o", None, None, "d"], ["-", "--", ":", "-"]
    ):
        r = [45, 112, 274, 797, 1355, 2547, 3600, 4500]  # r_g
        theta = np.ones_like(r) * float(inclination)
        # r = np.logspace(np.log10(45), np.log10(wind["r"].max()), n_points)
        for n in range(ax.size):
            i, j = np.unravel_index(n, ax.shape)
            color = "silver"
            plotted_lines = ax[i, j].plot(
                np.deg2rad(theta),
                np.log10(r),
                linestyle=linestyle,
                marker=marker,
                color=color,
                markerfacecolor=color,
                markeredgecolor="black",
                linewidth=1.5,
                markersize=7,
                # markeredgewidth=1.5,
                # alpha=0.7,
                label=r"$i = " + f"{inclination}" + r"^{\circ}$",
            )
            plotted_lines[0].set_markevery(slice(1, None))
            sight_lines_plotted.extend(plotted_lines)

    ax[0, 0].legend(loc="upper right", fontsize=10).set_zorder(0)

    # Plot each parameter on a sub-panel
    for n, (parameter, name, unit) in enumerate(
        zip(panel_parameters, panel_names, panel_units)
    ):
        i, j = np.unravel_index(n, (3, 2))
        if parameter == "He_i02_frac":
            fig, ax = wind.plot_parameter(
                parameter, fig=fig, ax=ax, a_idx=i, a_jdx=j, vmin=-8, vmax=-2
            )
        elif parameter == "h_density":
            thing = wind["H_i01_den"] + wind["H_i02_den"]
            fig, ax = wind.plot_parameter(thing, fig=fig, ax=ax, a_idx=i, a_jdx=j)
        else:
            fig, ax = wind.plot_parameter(parameter, fig=fig, ax=ax, a_idx=i, a_jdx=j)
        ax[i, j].grid(False)
        ax[i, j].set_title(r"$\log_{10}($" + name + r"$)$" + unit)
        # ax[i, j].set_rlim(np.log10(45), np.log10(wind.x_coords[-3]))
        ax[i, j].set_rlim(np.log10(45), np.log10(4000))
        ax[i, j].set_ylabel("")
        ax[i, j].set_thetamin(5)
        ax[i, j].set_thetamax(88)
        # ax[i, j].tick_params(labeltop=False, labelright=True)

    for line in sight_lines_plotted:
        line.set_zorder(10)
    for line in outlines_plotted:
        line.set_zorder(9)

    fig.text(0.0075, 0.5, r"$\log_{10}(r/r_{g})$", rotation="vertical")
    # fig = pysi.util.plot.finish_figure(fig)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.jpeg", dpi=300)

    if show_figure:
        plt.show()
    else:
        plt.close()

    wind.change_units(original_units)

    return fig, ax
