import os

import pysi
from matplotlib import pyplot as plt
from scipy import integrate


def integrate_spectrum(spectrum, sight, start, end):
    x = spectrum["Lambda"]
    y = spectrum[sight]  # * spectrum["Lambda"]
    x, y = pysi.util.array.get_subset_in_second_array(x, y, start, end)
    return integrate.simpson(y, x)


def plot_oxr(
    spec2d_full: pysi.Spectrum,
    spec2d_reduced: pysi.Spectrum,
    spec_created: pysi.Spectrum,
    display: bool,
) -> tuple[plt.Figure, plt.Axes]:
    spectra_wrapper = [
        {"Lambda": spec2d_full["Lambda"], "Flux": spec2d_full["87"]},
        {"Lambda": spec2d_full["Lambda"], "Flux": spec2d_full["77"]},
        {"Lambda": spec2d_full["Lambda"], "Flux": spec2d_full["56"]},
        {"Lambda": spec2d_full["Lambda"], "Flux": spec2d_full["34"]},
        {"Lambda": spec2d_full["Lambda"], "Flux": spec2d_full["14"]},
    ]

    inclination_angles = ["87", "77", "56", "34", "14"]

    fig, ax = plt.subplots(1, 1, figsize=(7, 6), sharex=True, sharey=True)

    for n, (spectrum, inclination) in enumerate(
        zip(spectra_wrapper, inclination_angles)
    ):
        y = spectrum["Flux"] * spectrum["Lambda"]
        ax.plot(spectrum["Lambda"], y, label=inclination + r"$^{\circ}$", alpha=0.75)

    print(
        spec_created["Created"] * spec_created["Lambda"],
    )

    spec_created.convert_flux_to_luminosity()
    ax.plot(
        spec_created["Lambda"],
        spec_created["Created"] * spec_created["Lambda"],
        color="k",
        linestyle="--",
        alpha=0.75,
        linewidth=3,
        zorder=1,
    )
    spec_created.convert_luminosity_to_flux()

    ax.legend(loc="upper right")
    ax = pysi.util.plot.set_axes_scales(ax, "loglog")
    ax.fill_between([1700, 6500], 1e41, 1e46, zorder=0, alpha=0.5, color="coral")
    ax.fill_between([10, 62], 1e41, 1e46, zorder=0, alpha=0.5, color="violet")
    ax.text(
        25,
        9e44,
        "X-ray band",
        rotation="vertical",
        ha="center",
        va="center",
        alpha=0.5,
    )
    ax.text(
        3500,
        9e44,
        "UV/Optical band",
        rotation="vertical",
        ha="center",
        va="center",
        alpha=0.5,
    )

    ax.set_ylim(1e42, 1e46)
    ax.set_xlim(10, 4e4)
    ax.set_xlabel(r"Rest-frame wavelength [\AA]")
    ax.set_ylabel(r"Luminosity $\lambda L_{\lambda}$ [ergs s$^{-1}$]")
    fig = pysi.util.plot.finish_figure(fig, hspace=0, wspace=0)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}-sed.pdf", dpi=300)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}-sed.jpeg", dpi=300)

    if display:
        plt.show()
    else:
        plt.close()

    fig, ax = plt.subplots(figsize=(7, 6))

    oxr_full = []
    oxr_reduced = []

    oxr_full_band2 = []
    oxr_reduced_band2 = []

    inclination_angles.reverse()

    optical_start = 1700
    optical_end = 6500

    for sight in inclination_angles:
        # full model
        lum_xray = integrate_spectrum(spec2d_full, sight, 1, 62)
        lum_optical = integrate_spectrum(spec2d_full, sight, optical_start, optical_end)
        oxr_full.append(lum_optical / lum_xray)
        # reduced model
        lum_xray = integrate_spectrum(spec2d_reduced, sight, 1, 62)
        lum_optical = integrate_spectrum(
            spec2d_reduced, sight, optical_start, optical_end
        )
        oxr_reduced.append(lum_optical / lum_xray)
        # full model - band 2
        lum_xray = integrate_spectrum(spec2d_full, sight, 1, 124)
        lum_optical = integrate_spectrum(spec2d_full, sight, optical_start, optical_end)
        oxr_full_band2.append(lum_optical / lum_xray)
        # reduced model - band 2
        lum_xray = integrate_spectrum(spec2d_reduced, sight, 1, 124)
        lum_optical = integrate_spectrum(
            spec2d_reduced, sight, optical_start, optical_end
        )
        oxr_reduced_band2.append(lum_optical / lum_xray)

    # Plot the first two sets of data
    (line1,) = ax.plot(
        inclination_angles,
        oxr_full,
        "o-",
        color="C0",
        markeredgecolor="k",
        label="Full [0.2 - 10 keV]",
    )
    (line2,) = ax.plot(
        inclination_angles,
        oxr_reduced,
        "D-",
        color="C1",
        markeredgecolor="k",
        label="Reduced [0.2 - 10 keV]",
    )

    # Plot the second two sets of data
    (line3,) = ax.plot(
        inclination_angles,
        oxr_full_band2,
        "o--",
        color="C0",
        markeredgecolor="k",
        label="Full [0.1 - 10 keV]",
    )
    (line4,) = ax.plot(
        inclination_angles,
        oxr_reduced_band2,
        "D--",
        color="C1",
        markeredgecolor="k",
        label="Reduced [0.1 - 10 keV]",
    )

    # Set the scales and labels
    ax = pysi.util.plot.set_axes_scales(ax, "logy")
    ax.set_xlabel("Inclination [$^{\circ}$]")
    ax.set_ylabel(r"$L_{\rm Opt} / L_{\rm X}$")

    # Create the first legend and place it in the upper left
    legend1 = ax.legend(handles=[line1, line2], loc="upper left")

    # Create the second legend and place it in the bottom right
    _legend2 = ax.legend(handles=[line3, line4], loc="lower right")

    # Add the first legend back to the plot
    ax.add_artist(legend1)

    # Finish the figure
    fig = pysi.util.plot.finish_figure(fig)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.jpeg", dpi=300)

    if display:
        plt.show()
    else:
        plt.close()

    return fig, ax
