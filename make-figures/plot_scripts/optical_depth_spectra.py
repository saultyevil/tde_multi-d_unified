import os

import numpy as np
import pysi
from astropy.constants import c
from matplotlib import pyplot as plt
from pysi.spec import Spectrum, labels


def plot_optical_depth_spectra(
    spec_full: Spectrum, spec_full_no_es: Spectrum, display: bool
) -> tuple[plt.Figure, plt.Axes]:
    """Create plots of the optical depth from the 2D models."""

    current1 = spec_full.current
    current2 = spec_full_no_es.current

    spec_full.set_spectrum("spec_tau")
    spec_full_no_es.set_spectrum("spec_tau")
    plot_xmax, plot_xmin = c.value / 3e14 / 1e-10, c.value / 1e18 / 1e-10

    fig, ax = plt.subplots(figsize=(7, 6))
    ax = pysi.util.plot.set_axes_scales(ax, "loglog")

    for n, inclination in enumerate(spec_full["inclinations"][::-1]):
        y = spec_full[inclination]
        if np.count_nonzero(y) == 0:
            continue
        ax.plot(
            spec_full["Lambda"],
            y,
            label=str(inclination) + r"$^{\circ}$",
            alpha=1,
            color=f"C{n}",
        )

    for n, inclination in enumerate(spec_full["inclinations"][::-1]):
        x = spec_full_no_es["Lambda"]
        y = spec_full_no_es[inclination]
        if np.count_nonzero(y) == 0:
            continue
        ax.plot(x, y, linestyle="--", alpha=0.5, color=f"C{n}", zorder=0)

    ax.legend(loc="lower left", fontsize=12)
    ax.set_xlim(plot_xmin, plot_xmax)
    ax.set_ylim(1e-5, 1e9)
    ax.set_xlabel("Rest-frame wavelength [\AA]")
    ax.set_ylabel(r"Continuum optical depth $\tau$")

    edges = [
        [r"O \textsc{vii}", 16],
        [r"C \textsc{v}", 32],
        [r"O \textsc{vi }", 93],
        [r"He \textsc{ii}", 227],
        [r"He \textsc{i}", 504],
        ["Lyman", 912],
        ["Balmer", 3646],
    ]

    # for i in range(len(edges)):
    #     edges[i][1] = pysi.constants.C / (edges[i][1] * pysi.constants.ANGSTROM)

    ax = labels.add_transition_labels_to_ax(
        ax, edges, label_linestyle="none", label_offset=0
    )
    fig = pysi.util.plot.finish_figure(fig)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.jpeg", dpi=300)

    if display:
        plt.show()
    else:
        plt.close()

    spec_full.set_spectrum(current1)
    spec_full.set_spectrum(current2)

    return fig, ax
