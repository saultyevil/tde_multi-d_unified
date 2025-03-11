import os

import pysi
from matplotlib import pyplot as plt
from scipy import integrate


def integrate_spectrum(spectrum, sight, start, end):
    x = spectrum["Lambda"]
    y = spectrum[sight]  # * spectrum["Lambda"]
    x, y = pysi.util.array.get_subset_in_second_array(x, y, start, end)
    return integrate.simpson(y, x=x)


def plot_oxr_multi_d(
    spectrum2d_full: pysi.Spectrum,
    spectrum2d_reduced: pysi.Spectrum,
    spectrum1d_bin1: pysi.Spectrum,
    spectrum1d_bin2: pysi.Spectrum,
    spectrum1d_bin3: pysi.Spectrum,
    spectrum1d_bin4: pysi.Spectrum,
    display: bool,
) -> tuple[plt.Figure, plt.Axes]:
    spherical_wrapper = {
        "77": spectrum1d_bin1,
        "56": spectrum1d_bin2,
        "34": spectrum1d_bin3,
        "14": spectrum1d_bin4,
    }

    spectrum2d_full.convert_flux_to_luminosity()
    spectrum2d_reduced.convert_flux_to_luminosity()

    inclination_angles = ["87", "77", "56", "34", "14"]

    fig, ax = plt.subplots(figsize=(7, 6))

    oxr_spherical = []
    oxr_reduced = []
    oxr_full = []

    xray_start = 0
    xray_end = 62

    optical_start = 1700
    optical_end = 6500

    inclination_angles = ["77", "56", "34", "14"]
    inclination_angles.reverse()

    for sight in inclination_angles:
        # # 1d model
        spherical_wrapper[sight].convert_flux_to_luminosity()
        lum_xray = integrate_spectrum(spherical_wrapper[sight], "45", 1, 62)
        lum_optical = integrate_spectrum(
            spherical_wrapper[sight], "45", optical_start, optical_end
        )
        spherical_wrapper[sight].convert_luminosity_to_flux()
        oxr_spherical.append(lum_optical / lum_xray)
        # 2d model - reduced
        lum_xray = integrate_spectrum(spectrum2d_reduced, sight, 1, 62)
        lum_optical = integrate_spectrum(
            spectrum2d_reduced, sight, optical_start, optical_end
        )
        oxr_reduced.append(lum_optical / lum_xray)
        # 2d model - full
        lum_xray = integrate_spectrum(spectrum2d_full, sight, 1, 62)
        lum_optical = integrate_spectrum(
            spectrum2d_full, sight, optical_start, optical_end
        )
        oxr_full.append(lum_optical / lum_xray)

    # Plot the first two sets of data

    (line1,) = ax.plot(
        inclination_angles,
        oxr_full,
        "D-",
        color="C0",
        markeredgecolor="k",
        label="Full",
    )
    (line2,) = ax.plot(
        inclination_angles,
        oxr_reduced,
        "D-",
        color="C1",
        markeredgecolor="k",
        label="Reduced",
    )
    (line3,) = ax.plot(
        inclination_angles,
        oxr_spherical,
        "o-",
        color="C2",
        markeredgecolor="k",
        label="D18",
    )

    # Set the scales and labels
    ax = pysi.util.plot.set_axes_scales(ax, "logy")
    ax.set_xlabel("Inclination [$^{\circ}$]")
    ax.set_ylabel(r"$L_{\rm Opt} / L_{\rm X}$")

    legend1 = ax.legend(handles=[line1, line2, line3], loc="upper left")
    ax.add_artist(legend1)

    # Finish the figure
    fig = pysi.util.plot.finish_figure(fig)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
    fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.jpeg", dpi=300)

    spectrum2d_full.convert_luminosity_to_flux()
    spectrum2d_reduced.convert_luminosity_to_flux()

    if display:
        plt.show()
    else:
        plt.close()

    return fig, ax
