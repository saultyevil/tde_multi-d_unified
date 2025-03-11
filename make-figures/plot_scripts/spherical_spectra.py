import os

import numpy
import pysi
import pysi.math
import pysi.math.blackbody
from pysi.spec import Spectrum
from matplotlib import pyplot as plt


def plot_spectra_spherical_models(
    my_spec: list[Spectrum],
    roth_model_spec: list[numpy.ndarray],
    panel_names: list[str],
    display: bool,
) -> tuple[plt.Figure, plt.Axes]:
    """Create a comparison plot between the spherical simulations from D18 and
    this work.
    """
    fig, axes = plt.subplots(2, 2, figsize=(9, 9), sharex=True, sharey=True)
    axes = axes.flatten()

    for n, (d18, p24, label) in enumerate(zip(roth_model_spec, my_spec, panel_names)):
        axes[n] = pysi.util.plot.set_axes_scales(axes[n], "loglog")
        axes[n].set_ylim(3e-4, 2.5)
        axes[n].text(0.03, 0.93, label, transform=axes[n].transAxes)

        # lambda * L_lambda
        axes[n].plot(d18[:, 0], d18[:, 1], label=r"\textsc{sedona}", alpha=0.75)

        # lambda * L_lambda
        p24.convert_flux_to_luminosity()
        y = p24["45"] * p24["Lambda"]
        axes[n].plot(p24["Lambda"], y / y.max(), label=r"\textsc{python}", alpha=0.75)

    wmin, wmax = 0.5, 4.7
    f_bb = pysi.math.blackbody.planck_lambda(1e6, numpy.logspace(wmin, wmax, 500))
    l_bb = 4 * numpy.pi * (my_spec[0]["distance"] * 3.086e18) ** 2 * f_bb
    y = numpy.logspace(wmin, wmax, 500) * l_bb

    for _, ax in enumerate(axes):
        ax.plot(
            numpy.logspace(wmin, wmax, 500),
            y / y.max(),
            color="k",
            linestyle="--",
            linewidth=3,
        )

    axes[0].set_xlim(5, 8e4)
    axes[1].legend(loc="upper right")
    fig.text(
        0.5,
        0.02,
        "Rest-frame wavelength [\AA]",
        ha="center",
        va="center",
    )
    fig.text(
        0.02,
        0.5,
        "$\lambda L_{\lambda}$ [erg s$^{-1}$]",
        ha="center",
        va="center",
        rotation="vertical",
    )
    fig = pysi.util.plot.finish_figure(fig, hspace=0, wspace=0)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}-sed.pdf", dpi=300)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}-sed.jpeg", dpi=300)

    if display:
        plt.show()
    else:
        plt.close()

    return fig, axes
