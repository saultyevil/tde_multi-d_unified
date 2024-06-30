import os

import pypython
from matplotlib import pyplot as plt

alpha = 0.75
smooth = 50

spec_full = pypython.Spectrum(
    "../simulations/2d/full/input.pf",
    distance=100 * 1e6,
    smooth=smooth,
    log=True,
)
spec_redu = pypython.Spectrum(
    "../simulations/2d/reduced/input.pf",
    distance=100 * 1e6,
    smooth=smooth,
    log=True,
)

created = pypython.Spectrum(
    "../simulations/2d/created/input.pf",
    distance=100 * 1e6,
    smooth=500,
    log=True,
)

spec_full.convert_flux_to_luminosity()
spec_redu.convert_flux_to_luminosity()
created.convert_flux_to_luminosity()

# just full spectra
for xrange, yrange, name in zip(
    [(18, 8e4)],
    [(7e41, 9e45)],
    ["sed"],
):
    fig, ax = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    ax = ax.flatten()

    for n, inclination in enumerate(["14", "34", "56", "77"][::-1]):
        y = spec_full[inclination] * spec_full["Lambda"]
        ax[n].plot(spec_full["Lambda"], y, label="Full", alpha=alpha, zorder=2)
        ax[n].text(
            0.02, 0.92, f"{inclination}" + r"$^{\circ}$", transform=ax[n].transAxes
        )
        ax[n].set_xlim(xrange[0], xrange[1])
        ax[n].set_ylim(yrange[0], yrange[1])
        ax[n] = pypython.plot.set_axes_scales(ax[n], "loglog")

    for my_ax in ax:
        my_ax.plot(
            created["Lambda"],
            created["Created"] * created["Lambda"],
            color="k",
            linestyle="--",
            alpha=alpha,
            linewidth=3,
            zorder=0,
        )

    fig.text(
        0.5,
        0.02,
        r"Rest-frame wavelength [\AA]",
        ha="center",
        va="center",
    )
    fig.text(
        0.02,
        0.5,
        "Luminosity $\lambda L_{\lambda}$ [ergs s$^{-1}$]",
        rotation="vertical",
        ha="center",
        va="center",
    )
    fig = pypython.plot.finish_figure(fig, wspace=0, hspace=0)
    fig.savefig(
        f"{os.path.splitext(os.path.basename(__file__))[0]}-{name}.pdf", dpi=300
    )
    plt.show()

# Reduced vs. full
for xrange, yrange, name in zip([(18, 9e4)], [(1e42, 3e46)], ["abun_sed"]):
    fig, ax = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    ax = ax.flatten()

    for n, inclination in enumerate(["14", "34", "56", "77"][::-1]):
        y = spec_full[inclination] * spec_full["Lambda"]
        ax[n].plot(spec_full["Lambda"], y, label="Full", alpha=alpha, zorder=2)

        y = spec_redu[inclination] * spec_redu["Lambda"]
        ax[n].plot(spec_redu["Lambda"], y, label="Reduced", alpha=alpha, zorder=1)

        ax[n].text(
            0.02, 0.92, f"{inclination}" + r"$^{\circ}$", transform=ax[n].transAxes
        )
        ax[n].set_xlim(xrange[0], xrange[1])
        ax[n].set_ylim(yrange[0], yrange[1])
        ax[n] = pypython.plot.set_axes_scales(ax[n], "loglog")

        ax[n].plot(
            created["Lambda"],
            created["Created"] * created["Lambda"],
            color="k",
            linestyle="--",
            alpha=alpha,
            linewidth=3,
            zorder=0,
        )

    ax[1].legend(loc="upper right")
    fig.text(
        0.5,
        0.02,
        r"Rest-frame wavelength [\AA]",
        ha="center",
        va="center",
    )
    fig.text(
        0.02,
        0.5,
        "Luminosity $\lambda L_{\lambda}$ [ergs s$^{-1}$]",
        rotation="vertical",
        ha="center",
        va="center",
    )
    fig = pypython.plot.finish_figure(fig, wspace=0, hspace=0)
    fig.savefig(
        f"{os.path.splitext(os.path.basename(__file__))[0]}-{name}.pdf", dpi=300
    )
    plt.show()
