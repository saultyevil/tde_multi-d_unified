import os

import pysi
from matplotlib import pyplot as plt

ALPHA = 0.75


def plot_multi_d_spectra(
    spec2d_reduced: pysi.Spectrum,
    spec2d_created: pysi.Spectrum,
    spectra1d: list[pysi.Spectrum],
    display: bool,
):
    inclinations = [
        "77",
        "56",
        "34",
        "14",
    ]

    names = ["Bin 1", "Bin 2", "Bin 3", "Bin 4"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    axes = axes.flatten()

    for n, spec1d in enumerate(spectra1d):
        # 1D model
        spec1d.convert_flux_to_luminosity()
        y = spec1d["45"] * spec1d["Lambda"]
        axes[n].plot(
            spec1d["Lambda"],
            y,
            color="C0",
            label="1D",
            alpha=ALPHA,
        )
        spec1d.convert_luminosity_to_flux()

        # 2D model
        spec2d_reduced.convert_flux_to_luminosity()

        y = spec2d_reduced[inclinations[n]] * spec2d_reduced["Lambda"]
        axes[n].plot(
            spec2d_reduced["Lambda"],
            y,
            color="C1",
            label="2D",
            alpha=ALPHA,
        )

        # Names
        axes[n].text(0.03, 0.93, names[n], transform=axes[n].transAxes)
        axes[n].text(
            0.03,
            0.86,
            f"{inclinations[n]}" + r"$^{\circ}$",
            transform=axes[n].transAxes,
        )
        axes[n] = pysi.util.plot.set_axes_scales(axes[n], "loglog")

    spec2d_reduced.convert_luminosity_to_flux()

    # spec_created is already in luminosity, so we can leave that for now
    spec2d_created.convert_flux_to_luminosity()
    for ax in axes:
        ax.plot(
            spec2d_created["Lambda"],
            spec2d_created["Created"] * spec2d_created["Lambda"],
            color="k",
            linestyle="--",
            alpha=ALPHA,
            linewidth=3,
            zorder=0,
        )
    spec2d_created.convert_luminosity_to_flux()

    axes[1].legend()
    axes[0].set_ylim(2e42, 9e45)
    axes[0].set_xlim(10, 5e4)

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
        r"Luminosity $\lambda L_{\lambda}$ [ergs s$^{-1}$]",
        ha="center",
        va="center",
        rotation="vertical",
    )

    fig = pysi.util.plot.finish_figure(fig, wspace=0, hspace=0)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.jpeg", dpi=300)

    if display:
        plt.show()
    else:
        plt.close()

    return fig, ax
