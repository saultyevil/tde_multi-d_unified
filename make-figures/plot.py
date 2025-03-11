"""Create all plots for the paper.

All wind saves and spectra are loaded in this script and passed to the each
plotting function.
"""

import os
import shutil

import numpy
from astropy.table import Table
from pysi.spec import Spectrum
from pysi.util.plot import set_figure_style
from pysi.util.shell import find_file_with_pattern
from pysi.wind import Wind

from plot_scripts.cell_seds import plot_model_cell_seds
from plot_scripts.model_properties import plot_model_properties
from plot_scripts.model_spectra import plot_model_spectra
from plot_scripts.multi_d_cell_seds import plot_multi_d_cell_seds
from plot_scripts.multi_d_spectra import plot_multi_d_spectra
from plot_scripts.optical_depth_spectra import plot_optical_depth_spectra
from plot_scripts.optical_depth_surfaces import plot_optical_depth_surfaces
from plot_scripts.oxr import plot_oxr
from plot_scripts.oxr_multi_d import plot_oxr_multi_d
from plot_scripts.regrid_comparison import plot_regrid_comparison
from plot_scripts.spherical_bins import plot_spherical_bins
from plot_scripts.spherical_properties import plot_spherical_properties
from plot_scripts.spherical_spectra import plot_spectra_spherical_models


def collect_plots() -> None:
    """Move all created plots to separate directories."""
    for directory in ["jpeg", "pdf"]:
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.makedirs(directory)
    for file in os.listdir("."):
        if file.endswith(".jpeg"):
            shutil.move(file, "jpeg")
        elif file.endswith(".pdf"):
            shutil.move(file, "pdf")


DEFAULT_SCALE = "log"
SMOOTH_WIDTH = 50
SHOW_FIGURES = False

wind1d_bin1 = Wind("input", "../simulations/1d/Bin-1", version="88")
wind1d_bin2 = Wind("input", "../simulations/1d/Bin-2", version="88")
wind1d_bin3 = Wind("input", "../simulations/1d/Bin-3", version="88")
wind1d_bin4 = Wind("input", "../simulations/1d/Bin-4", version="88")
wind2d_full = Wind("input", "../simulations/2d/full", version="88")
wind2d_reduced = Wind("input", "../simulations/2d/reduced", version="88")

spectrum2d_full = Spectrum("input", "../simulations/2d/full")
spectrum2d_reduced = Spectrum("input", "../simulations/2d/reduced")
spectrum2d_full_no_es = Spectrum("input_no_es", "../simulations/2d/full")
spectrum2d_created = Spectrum("input", "../simulations/2d/reduced", smooth_width=500)
spectrum2d_full.convert_flux_to_luminosity()
spectrum2d_reduced.convert_flux_to_luminosity()
spectrum1d_bin1 = Spectrum(
    "input",
    "../simulations/1d/Bin-1",
    smooth_width=SMOOTH_WIDTH,
    default_scale=DEFAULT_SCALE,
)
spectrum1d_bin2 = Spectrum(
    "input",
    "../simulations/1d/Bin-2",
    smooth_width=SMOOTH_WIDTH,
    default_scale=DEFAULT_SCALE,
)
spectrum1d_bin3 = Spectrum(
    "input",
    "../simulations/1d/Bin-3",
    smooth_width=SMOOTH_WIDTH,
    default_scale=DEFAULT_SCALE,
)
spectrum1d_bin4 = Spectrum(
    "input",
    "../simulations/1d/Bin-4",
    smooth_width=SMOOTH_WIDTH,
    default_scale=DEFAULT_SCALE,
)
spectrum_benchmark_bin1 = Spectrum(
    "input",
    "../simulations/benchmark/Bin-1",
    smooth_width=SMOOTH_WIDTH,
    default_scale=DEFAULT_SCALE,
)
spectrum_benchmark_bin2 = Spectrum(
    "input",
    "../simulations/benchmark/Bin-2",
    smooth_width=SMOOTH_WIDTH,
    default_scale=DEFAULT_SCALE,
)
spectrum_benchmark_bin3 = Spectrum(
    "input",
    "../simulations/benchmark/Bin-3",
    smooth_width=SMOOTH_WIDTH,
    default_scale=DEFAULT_SCALE,
)
spectrum_benchmark_bin4 = Spectrum(
    "input",
    "../simulations/benchmark/Bin-4",
    smooth_width=SMOOTH_WIDTH,
    default_scale=DEFAULT_SCALE,
)

d18_spectra = [
    numpy.loadtxt(
        "../scripts/m_original/s_dai67-87.txt", delimiter=","
    ),  # The file Nathan sent was empty, so use something extracted from Dai et al. 2018
    numpy.loadtxt("../dai-models/bin2.txt"),
    numpy.loadtxt("../dai-models/bin3.txt"),
    numpy.loadtxt("../dai-models/bin4.txt"),
]

d18_rho = find_file_with_pattern("*rho_*", "../scripts/m_original")
d18_vel = find_file_with_pattern("*vel_*", "../scripts/m_original")
d18_2d_grid = Table.read(
    "../scripts/m_original/quantities2.txt",
    format="ascii",
    names=("r", "theta", "rho", "t_gas", "t_rad", "v_r", "v_theta"),
)

set_figure_style()

try:
    plot_model_properties(wind2d_full, SHOW_FIGURES)
    plot_model_spectra(
        spectrum2d_full, spectrum2d_reduced, spectrum2d_created, SHOW_FIGURES
    )
    plot_optical_depth_spectra(spectrum2d_full, spectrum2d_full_no_es, SHOW_FIGURES)
    plot_model_cell_seds(wind2d_full, SHOW_FIGURES)
    plot_multi_d_cell_seds(wind1d_bin4, wind2d_reduced, SHOW_FIGURES)
    plot_spectra_spherical_models(
        (
            spectrum_benchmark_bin1,
            spectrum_benchmark_bin2,
            spectrum_benchmark_bin3,
            spectrum_benchmark_bin4,
        ),
        d18_spectra,
        ("Bin 1", "Bin 2", "Bin 3", "Bin 4"),
        SHOW_FIGURES,
    )
    plot_regrid_comparison(wind2d_full, d18_2d_grid, d18_rho, d18_vel, SHOW_FIGURES)
    plot_optical_depth_surfaces(wind2d_full, SHOW_FIGURES)
    plot_spherical_bins(d18_2d_grid, SHOW_FIGURES)
    plot_spherical_properties(SHOW_FIGURES)
    plot_multi_d_spectra(
        spectrum2d_reduced,
        spectrum2d_created,
        [spectrum1d_bin1, spectrum1d_bin2, spectrum1d_bin3, spectrum1d_bin4],
        SHOW_FIGURES,
    )
    plot_oxr(spectrum2d_full, spectrum2d_reduced, spectrum2d_created, SHOW_FIGURES)
    plot_oxr_multi_d(
        spectrum2d_full,
        spectrum2d_reduced,
        spectrum_benchmark_bin1,
        spectrum_benchmark_bin2,
        spectrum_benchmark_bin3,
        spectrum_benchmark_bin4,
        SHOW_FIGURES,
    )
finally:
    collect_plots()
