import os

import numpy as np
import pypython
from astropy.table import Table
from matplotlib import pyplot as plt

sight_lines = [
    7 * np.pi / 16,
    5 * np.pi / 16,
    np.pi / 4,
    3 * np.pi / 16,
    np.pi / 8,
    np.pi / 16,
]


def extract_from_original_grid(grid, what, angle):
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

    return radius[plot_start:], thing[plot_start:]


alpha = 1

dai_rho_models = pypython.find("*rho_*", "../other-data")
dai_vel_models = pypython.find("*vel_*", "../other-data")
original = Table.read(
    "../other-data/quantities2.txt",
    format="ascii",
    names=("r", "theta", "rho", "t_gas", "t_rad", "v_r", "v_theta"),
)
wind = pypython.Wind(
    "input",
    "../simulations/2d/full",
    distance_units="rg",
    velocity_units="c",
    masked=False,
)

fig, ax = plt.subplots(2, 1, figsize=(7, 9), sharex=True)

for n, sight_line in enumerate(sight_lines):
    x, z, w = wind.get_variable_along_sight_line(np.rad2deg(sight_line), "rho")

    r, w = pypython.get_xy_subset(np.sqrt(x**2 + z**2), w, 0, 3500)
    ax[0].plot(
        r[1:-1],
        w[1:-1],
        linewidth=3,
        zorder=0,
        label=f"{np.rad2deg(sight_line):.1f}" + r"$^{\circ}$",
    )

    x, z, w = wind.get_variable_along_sight_line(np.rad2deg(sight_line), "v_l")
    r, w = pypython.get_xy_subset(np.sqrt(x**2 + z**2), w, 0, 3500)
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
        v_r = pypython.smooth_array(v_r, 5)
        r, v_r = r[2:], v_r[2:]
        r, v_r = pypython.get_xy_subset(r[2:], v_r[2:], 0, 3000)
        data = np.zeros((len(v_r), 2))
        data[:, 0] = r
        data[:, 1] = v_r
    else:
        data = np.loadtxt(model, delimiter=",")

    ax[1].plot(data[:, 0], data[:, 1], "--", color=f"C{n}", alpha=alpha)

ax[0] = pypython.plot.set_axes_scales(ax[0], "logy")
ax[0].legend(loc="upper right", ncol=1)

fig = pypython.plot.finish_figure(fig, hspace=0)
fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)

plt.show()
