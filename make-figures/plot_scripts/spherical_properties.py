import os

import pandas as pd
import pysi
from astropy.constants import c
from astropy.table import Table
from matplotlib import pyplot as plt


def plot_spherical_properties(display: bool) -> tuple[plt.Figure, plt.Axes]:
    rg = pysi.math.blackhole.gravitational_radius(5e6)
    files = {
        "Bin 1": [
            "../dai-models/Envelope_bin1.txt",
            "../simulations/benchmark/Bin-1/Bin-1.txt",
        ],
        "Bin 2": [
            "../dai-models/Envelope_bin2.txt",
            "../simulations/benchmark/Bin-2/Bin-2.txt",
        ],
        "Bin 3": [
            "../dai-models/Envelope_bin3.txt",
            "../simulations/benchmark/Bin-3/Bin-3.txt",
        ],
        "Bin 4": [
            "../dai-models/Envelope_bin4.txt",
            "../simulations/benchmark/Bin-4/Bin-4.txt",
        ],
    }
    roth_model_headers = ["r_in", "r_out", "rho", "T_init", "v_in", "v_out"]

    fig, ax = plt.subplots(2, 1, figsize=(7, 9), sharex=True)
    ax = ax.flatten()

    for n, (theta, (roth_model, me)) in enumerate(files.items()):
        me = Table.read(me, format="ascii.fixed_width_two_line")
        roth_model = pd.read_csv(
            roth_model, delim_whitespace=True, skiprows=2, names=roth_model_headers
        )

        ax[0].plot(
            me["r"][1:-2] / rg,
            me["rho"][1:-2],
            "-",
            color=f"C{n}",
            label=f"Bin {n + 1}",
        )
        ax[0].plot(roth_model["r_in"] / rg, roth_model["rho"], "--", color=f"C{n}")

        ax[1].plot(
            me["r"][1:-2] / rg,
            me["v_r"][1:-2] / c.cgs.value,
            "-",
            color=f"C{n}",
        )
        ax[1].plot(
            roth_model["r_in"] / rg,
            roth_model["v_in"] / c.cgs.value,
            "--",
            color=f"C{n}",
        )

    ax[0].legend()
    ax[0] = pysi.util.plot.set_axes_scales(ax[0], "loglog")
    ax[1] = pysi.util.plot.set_axes_scales(ax[1], "logx")

    ax[0].set_ylabel(r"$\log_{10}(\rho)$ [g cm$^{-3}$]")
    ax[1].set_ylabel(r"$v_{r} / c$")
    ax[1].set_xlabel(r"$r / r_{g}$")

    fig = pysi.util.plot.finish_figure(fig, hspace=0)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.jpeg", dpi=300)

    if display:
        plt.show()
    else:
        plt.close()

    return fig, ax
