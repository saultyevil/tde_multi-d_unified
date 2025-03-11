# TODO(EP): this needs cleaning up more


import os
from subprocess import CalledProcessError

import numpy as np
import pysi
from matplotlib import pyplot as plt

REUSE_DATA = True


def read_in_photosphere_locations(
    tau: str, root: str, fp: str
) -> tuple[float, np.ndarray]:
    """Read in the photosphere file."""

    with open(f"{fp}/{root}_{tau}.photosphere", "r", encoding="utf-8") as file_in:
        lines = file_in.readlines()

    tau_from_file = 0
    surfaces_from_file = []

    for line in lines:
        if line.startswith("# Electron scatter photosphere"):
            tau_from_file = float(line.split()[-1])
        if line.startswith("#"):
            continue
        surfaces_from_file.append(line.split())
    surfaces_from_file = np.array(surfaces_from_file[1:], dtype=np.float64)

    return tau_from_file, surfaces_from_file


def plot_optical_depth_surfaces(wind: pysi.Wind, display):
    original_units = wind.distance_units
    wind.change_units(pysi.wind.DistanceUnits.GRAVITATIONAL_RADIUS)

    fig, axes = plt.subplots(
        1, 2, figsize=(12.5, 5), subplot_kw={"projection": "polar"}
    )

    im_ax0 = axes[0].pcolormesh(
        np.deg2rad(wind["theta"]),
        np.log10(wind["r"]),  # in units of r / rg
        np.log10(wind["rho"]),
        shading="auto",
        alpha=1,
        linewidth=0,
        rasterized=True,
    )
    im_ax1 = axes[1].pcolormesh(
        np.deg2rad(wind["theta"]),
        np.log10(wind["r"]),  # in units of r / rg
        np.log10(wind["t_e"]),
        shading="auto",
        alpha=1,
        linewidth=0,
        rasterized=True,
    )

    fig.colorbar(im_ax0, ax=axes[0])
    fig.colorbar(im_ax1, ax=axes[1])
    axes[0].set_title(r"$\log_{10}(\rho)$ [g cm$^{-3}$]")
    axes[1].set_title(r"$\log_{10}(T_{\mathrm{e}})$ [K]")
    surface_colours = [1, 3, 4, 5, 6, 7, 8]

    # In each panel, plot sight lines and markers for where cells in seds are
    for inclination, marker, linestyle in zip(
        ["14", "34", "56", "77"], ["o", None, None, "d"], ["-", "--", ":", "-"]
    ):
        r = [45, 112, 274, 797, 1355, 2547, 3600, 4500]  # r_g
        theta = np.ones_like(r) * float(inclination)
        for n in range(axes.size):
            color = "silver"
            plotted_lines = axes[n].plot(
                np.deg2rad(theta),
                np.log10(r),
                linestyle=linestyle,
                marker=marker,
                color=color,
                markerfacecolor=color,
                markeredgecolor="black",
                linewidth=1.5,
                markersize=7,
                label=r"$i = " + f"{inclination}" + r"^{\circ}$" if n == 0 else "",
            )
            plotted_lines[0].set_markevery(slice(1, None))

    axes[0].legend(loc="upper right", fontsize=11, ncol=1).set_zorder(0)

    for n, tau in enumerate([1, 5, 10, 50, 100]):
        if REUSE_DATA:
            try:
                filename = f"../simulations/2d/full/input_{tau}.photosphere"
                locations = np.loadtxt(filename, comments="#", skiprows=7)
            except OSError:
                try:
                    pysi.util.run.run_py_optical_depth(
                        wind.root, wind.directory, scatter_surface=tau
                    )
                except CalledProcessError:
                    print(f"Failed to run py_optical_depth for tau_es = {tau}")
                    continue
        else:
            try:
                pysi.util.run.run_py_optical_depth(
                    wind.root, wind.directory, scatter_surface=tau
                )
            except CalledProcessError:
                print(f"Failed to run py_optical_depth for tau_es = {tau}")
                continue

            err = pysi.util.run.run_py_optical_depth(
                wind.root, wind.directory, scatter_surface=tau
            )
            if err:
                print(f"error return of {err} from py_optical_depth for tau_es = {tau}")
                continue

        tau_es, locations = read_in_photosphere_locations(
            tau, wind.root, wind.directory
        )
        r = np.sqrt(locations[:, 0] ** 2 + locations[:, 2] ** 2)
        theta = np.arctan2(locations[:, 2], locations[:, 0])

        sm = 5
        r, theta = (
            pysi.util.array.smooth_array(r, sm),
            pysi.util.array.smooth_array(theta, sm),
        )
        for ax in axes:
            ax.plot(
                np.pi / 2 - theta,
                np.log10(r / wind.grav_radius),
                label=r"$\tau_{\rm es} =" + f"{tau_es:.0f}" + "$",
                color=f"C{surface_colours[n]}",
            )

    for ax in axes:
        ax.grid(False)
        ax.set_rlim(np.log10(45), np.log10(wind.x_coords[-3]))
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_thetamin(5)
        ax.set_thetamax(88)

    axes[0].set_ylabel(r"$\log_{10}(r / r_{g})$")
    axes[0].set_rlabel_position(90)
    axes[1].legend(loc="upper right", fontsize=11, ncol=1).set_zorder(0)

    fig = pysi.util.plot.finish_figure(fig)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.jpeg", dpi=300)

    wind.change_units(original_units)

    if display:
        plt.show()
    else:
        plt.close()

    return fig, ax
