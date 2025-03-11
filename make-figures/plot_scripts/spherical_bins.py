import os

import astropy
import astropy.table
import numpy as np
import pysi
from matplotlib import pyplot as plt

ALPHA = 0.8


def plot_spherical_bins(
    dai2d_grid: astropy.table.Table, display: bool
) -> tuple[plt.Figure, plt.Axes]:
    nx, nz = 128, 64

    bins = [
        [67.5, 87.4, "Bin 1"],
        [45, 67.5, "Bin 2"],
        [22.5, 45, "Bin 3"],
        [5.7, 22.5, "Bin 4"],
    ]

    fig, ax = plt.subplots(figsize=(8.5, 6), subplot_kw={"projection": "polar"})

    ax.grid(False)
    im = ax.pcolormesh(
        dai2d_grid["theta"].reshape(nx, nz),
        dai2d_grid["r"].reshape(nx, nz),
        dai2d_grid["rho"].reshape(nx, nz),
        linewidth=0,
        rasterized=True,
    )
    ax.set_xlabel(r"$r / r_{g}$")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    ax.set_rlabel_position(90)
    ax.set_rlim(0, 2000)

    cbar = plt.colorbar(im, ax=ax, location="left")
    cbar.set_label(r"$\log_{10}(\rho)$" + " [g cm$^{-3}$]")

    for limits in bins:
        theta1, theta2, name = limits[0], limits[1], limits[2]
        x_coords = np.logspace(0, np.log10(3000), 50)
        theta_coords1, theta_coords2 = (
            np.ones_like(x_coords) * theta1,
            np.ones_like(x_coords) * theta2,
        )
        ax.plot(np.deg2rad(theta_coords1), x_coords, "k--", alpha=ALPHA)
        ax.plot(np.deg2rad(theta_coords2), x_coords, "k--", alpha=ALPHA)
        theta_cen = 0.5 * (theta2 + theta1)
        ax.text(
            np.deg2rad(theta_cen),
            1100,
            name,
            color="k",
            rotation=90 - theta_cen,
            va="center",
            ha="center",
        )

    for inclination, marker, linestyle in zip(
        ["14", "34", "56", "77"], ["o", None, None, "d"], ["-", "--", ":", "-"]
    ):
        r = [45, 112, 274, 797, 1355, 2547, 3600, 4500]  # r_g
        theta = np.ones_like(r) * float(inclination)
        color = "silver"
        plotted_lines = ax.plot(
            np.deg2rad(theta),
            r,
            linestyle=linestyle,
            marker=marker,
            color=color,
            markerfacecolor=color,
            markeredgecolor="black",
            linewidth=1.5,
            markersize=7,
            alpha=0.7,
            label=r"$i = " + f"{inclination}" + r"^{\circ}$",
        )
        plotted_lines[0].set_markevery(slice(1, None))

    ax.legend(loc="upper right", fontsize=12, ncol=1)

    fig = pysi.util.plot.finish_figure(fig)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.jpeg", dpi=300)

    if display:
        plt.show()
    else:
        plt.close()

    return fig, ax
