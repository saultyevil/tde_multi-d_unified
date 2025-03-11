import os

import numpy as np
import pysi
from matplotlib import pyplot as plt
from pysi.wind import Wind
from pysi.wind.enum import DistanceUnits, VelocityUnits

sight_lines = [
    7 * np.pi / 16,
    5 * np.pi / 16,
    np.pi / 4,
    3 * np.pi / 16,
    np.pi / 8,
    np.pi / 16,
]


def extract_from_original_grid(grid, what, angle, *, ignore_plot_start=False):
    """Get something from the unstructured grid."""

    thing = []

    mdims = 128
    ndims = 64

    plot_start = -1

    for i in range(mdims):
        start = i * ndims
        end = (i + 1) * ndims
        t = grid[start:end]
        j = 0
        while t["theta"][j] < angle:
            j += 1
        thing.append(t[what][j])
        if t["r"][j] > 100 and plot_start < 0:
            plot_start = i

    radius = grid["r"].data[::ndims]

    if ignore_plot_start:
        plot_start = 0

    return radius[plot_start:], thing[plot_start:]


def plot_regrid_comparison(
    wind: Wind, original, dai_rho_models, dai_vel_models, display, *, alpha=1
):
    fig, ax = plt.subplots(2, 1, figsize=(7, 9), sharex=True)
    original_distance_units = wind.distance_units
    original_velocity_units = wind.velocity_units
    wind.change_units(DistanceUnits.GRAVITATIONAL_RADIUS)
    wind.change_units(VelocityUnits.SPEED_OF_LIGHT)

    for n, sight_line in enumerate(sight_lines):
        x, z, w = wind.get_variable_along_sight_line(np.rad2deg(sight_line), "rho")
        r = np.unique(wind["r"])
        r, w = pysi.util.array.get_subset_in_second_array(r, w, 0, 3500)
        ax[0].plot(
            r[1:-1],
            w[1:-1],
            linewidth=3,
            zorder=0,
            label=f"{np.rad2deg(sight_line):.1f}" + r"$^{\circ}$",
        )

        x, z, w = wind.get_variable_along_sight_line(np.rad2deg(sight_line), "v_l")
        rrr = np.unique(wind["r"])
        r, w = pysi.util.array.get_subset_in_second_array(rrr, w, 0, 3500)

        ax[1].plot(
            r[1:-1],
            w[1:-1],
            linewidth=3,
            zorder=0,
            label=f"{np.rad2deg(sight_line):.1f}" + r"$^{\circ}$",
        )

    ax[0].set_ylabel(r"$\log_{10}(\rho)$ [g cm$^{-3}$]")
    ax[1].set_ylabel(r"$v_{r} / c$")
    ax[1].set_xlabel(r"$r / r_{g}$")
    ax[1].set_xscale("log")

    for n, model in enumerate(dai_rho_models):
        data = np.loadtxt(model, delimiter=",")
        ax[0].plot(data[:, 0], data[:, 1], "--", color=f"C{n}", alpha=alpha)

    for n, model in enumerate(dai_vel_models):
        if "6_vel_pi_16.txt" in model:
            r, v_r = extract_from_original_grid(original, "v_r", sight_line)
            v_r = pysi.util.array.smooth_array(v_r, 5)
            r, v_r = r[2:], v_r[2:]
            r, v_r = pysi.util.array.get_subset_in_second_array(r[2:], v_r[2:], 0, 3000)
            data = np.zeros((len(v_r), 2))
            data[:, 0] = r
            data[:, 1] = v_r
        else:
            data = np.loadtxt(model, delimiter=",")

        ax[1].plot(data[:, 0], data[:, 1], "--", color=f"C{n}", alpha=alpha)

    ax[0] = pysi.util.plot.set_axes_scales(ax[0], "logy")
    ax[0].legend(loc="upper right", ncol=1)

    fig = pysi.util.plot.finish_figure(fig, hspace=0)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.jpeg", dpi=300)

    if display:
        plt.show()
    else:
        plt.close()

    wind.change_units(original_distance_units)
    wind.change_units(original_velocity_units)

    return fig, ax
