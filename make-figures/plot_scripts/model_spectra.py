import os

import pysi
from matplotlib import pyplot as plt
from pysi.spec import Spectrum


def plot_model_spectra(
    spec_full: Spectrum,
    spec_reduced: Spectrum,
    spec_created: Spectrum,
    show_figure: bool,
) -> tuple[plt.Figure, plt.Axes]:
    """Create plots of the spectra from the 2D models."""

    alpha = 0.75
    spec_created.convert_flux_to_luminosity()

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
            ax[n] = pysi.util.plot.set_axes_scales(ax[n], "loglog")

        for my_ax in ax:
            my_ax.plot(
                spec_created["Lambda"],
                spec_created["Created"] * spec_created["Lambda"],
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
        fig = pysi.util.plot.finish_figure(fig, wspace=0, hspace=0)
        fig.savefig(
            f"{os.path.splitext(os.path.basename(__file__))[0]}-{name}.pdf", dpi=300
        )

        if show_figure:
            plt.show()
        else:
            plt.close()

    # Reduced vs. full
    for xrange, yrange, name in zip([(18, 9e4)], [(1e42, 3e46)], ["abun_sed"]):
        fig, ax = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
        ax = ax.flatten()

        for n, inclination in enumerate(["14", "34", "56", "77"][::-1]):
            y = spec_full[inclination] * spec_full["Lambda"]
            ax[n].plot(spec_full["Lambda"], y, label="Full", alpha=alpha, zorder=2)

            y = spec_reduced[inclination] * spec_reduced["Lambda"]
            ax[n].plot(
                spec_reduced["Lambda"], y, label="Reduced", alpha=alpha, zorder=1
            )

            ax[n].text(
                0.02, 0.92, f"{inclination}" + r"$^{\circ}$", transform=ax[n].transAxes
            )
            ax[n].set_xlim(xrange[0], xrange[1])
            ax[n].set_ylim(yrange[0], yrange[1])
            ax[n] = pysi.util.plot.set_axes_scales(ax[n], "loglog")

            ax[n].plot(
                spec_created["Lambda"],
                spec_created["Created"] * spec_created["Lambda"],
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
        fig = pysi.util.plot.finish_figure(fig, wspace=0, hspace=0)
        fig.savefig(
            f"{os.path.splitext(os.path.basename(__file__))[0]}-{name}.pdf", dpi=300
        )
        fig.savefig(
            f"{os.path.splitext(os.path.basename(__file__))[0]}-{name}.jpeg", dpi=300
        )
        if show_figure:
            plt.show()
        else:
            plt.close()

    spec_created.convert_luminosity_to_flux()

    return fig, ax
