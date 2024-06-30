import os
import subprocess

import numpy as np
import pypython
from matplotlib import pyplot as plt

REUSE_DATA = False


def run_py_optical_depth(root: str, fp: str, tau: str):
    """Run py_optical_depth to get the photosphere."""

    command = (
        f"cd {fp}; py_optd -p {tau} --smax 0.1 {root}; cp {root}.photosphere {root}_{tau}.photosphere",
    )

    print(command[0])
    sh = subprocess.run(command, shell=True, check=False, capture_output=True)

    return sh.returncode


def read_in_photosphere_locations(tau: str, root: str, fp: str):
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


rg = pypython.physics.blackhole.gravitational_radius(5e6)
wind = pypython.Wind("input", "../simulations/2d/full", distance_units="rg")

fig, axes = plt.subplots(1, 2, figsize=(12, 5), subplot_kw={"projection": "polar"})

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

cbar_ax0 = fig.colorbar(im_ax0, ax=axes[0])
cbar_ax1 = fig.colorbar(im_ax1, ax=axes[1])

axes[0].set_title(r"$\log_{10}(\rho)$ [g cm$^{-3}$]")
axes[1].set_title(r"$\log_{10}(T_{\mathrm{e}})$ [K]")

surface_colours = [1, 3, 4, 5, 6, 7, 8]

for n, tau in enumerate([1, 5, 10, 50, 100]):
    if REUSE_DATA:
        try:
            filename = f"../simulations/2d/full/input_{tau}.photosphere"
            locations = np.loadtxt(filename, comments="#", skiprows=7)
            print("found", filename)
        except IOError:
            err = run_py_optical_depth(wind.root, wind.fp, tau)
            if err:
                print(f"error return of {err} from py_optical_depth for tau_es = {tau}")
                continue
    else:
        err = run_py_optical_depth(wind.root, wind.fp, tau)
        if err:
            print(f"error return of {err} from py_optical_depth for tau_es = {tau}")
            continue

    tau_es, locations = read_in_photosphere_locations(tau, wind.root, wind.fp)
    r = np.sqrt(locations[:, 0] ** 2 + locations[:, 2] ** 2)
    theta = np.arctan2(locations[:, 2], locations[:, 0])

    sm = 5
    r, theta = pypython.smooth_array(r, sm), pypython.smooth_array(theta, sm)
    for ax in axes:
        ax.plot(
            np.pi / 2 - theta,
            np.log10(r / rg),
            label=r"$\tau_{\rm es} =" + f"{tau_es:.0f}" + "$",
            color=f"C{surface_colours[n]}",
        )

for ax in axes:
    ax.set_rlim(np.log10(45), np.log10(wind.x_axis_coords[-3]))
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_thetamin(5)
    ax.set_thetamax(88)

axes[0].set_ylabel(r"$\log_{10}(r / r_{g})$")
axes[0].set_rlabel_position(90)
axes[1].legend(loc="upper right", fontsize=11, ncol=1).set_zorder(0)

fig = pypython.plot.finish_figure(fig)
fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
plt.show()
