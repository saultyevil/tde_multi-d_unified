import os

import numpy as np
import pypython
from matplotlib import pyplot as plt

roth_model_spec = [
    np.loadtxt(
        "../other-data/s_dai67-87.txt", delimiter=","
    ),  # The file Nathan sent was empty, so use something extracted from Dai et al. 2018
    np.loadtxt("../other-data/bin2.txt"),
    np.loadtxt("../other-data/bin3.txt"),
    np.loadtxt("../other-data/bin4.txt"),
]

log_spec = True
smooth = 20
my_spec = [
    pypython.Spectrum(
        "../simulations/benchmark/Bin-1/input.pf", log=log_spec, smooth=smooth
    ),
    pypython.Spectrum(
        "../simulations/benchmark/Bin-2/input.pf", log=log_spec, smooth=smooth
    ),
    pypython.Spectrum(
        "../simulations/benchmark/Bin-3/input.pf", log=log_spec, smooth=smooth
    ),
    pypython.Spectrum(
        "../simulations/benchmark/Bin-4/input.pf", log=log_spec, smooth=smooth
    ),
]


names = ["Bin 1", "Bin 2", "Bin 3", "Bin 4"]

lines = [
    [r"P \textsc{v}", 1118],
    [r"Ly$\alpha$ / N \textsc{v}", 1216],
    [r"O \textsc{v} / Si \textsc{iv}", 1371],
    [r"C \textsc{iv}", 1548],
    [r"H$_{\delta}$", 4101],
    [r"H$_{\gamma}$", 4340],
    [r"He \textsc{ii}", 4686],
    [r"H$_{\beta}$", 4861],
    [r"Na \textsc{i}", 5891],
    [r"H$_{\alpha}$", 6564],
]

#
# Each model in its own panel
#

fig, axes = plt.subplots(2, 2, figsize=(9, 9), sharex=True, sharey=True)
axes = axes.flatten()

for n, (nathan, mine, name) in enumerate(zip(roth_model_spec, my_spec, names)):
    axes[n] = pypython.plot.set_axes_scales(axes[n], "loglog")
    axes[n].set_ylim(3e-4, 2.5)
    axes[n].text(0.03, 0.93, name, transform=axes[n].transAxes)

    # lambda * L_lambda
    axes[n].plot(nathan[:, 0], nathan[:, 1], label=r"\textsc{sedona}", alpha=0.75)

    # lambda * L_lambda
    mine.convert_flux_to_luminosity()
    y = mine["45"] * mine["Lambda"]
    axes[n].plot(mine["Lambda"], y / y.max(), label=r"\textsc{python}", alpha=0.75)

wmin, wmax = 0.5, 4.7
f_bb = pypython.physics.blackbody.planck_lambda(1e6, np.logspace(wmin, wmax, 500))
l_bb = 4 * np.pi * (my_spec[0].distance * 3.086e18) ** 2 * f_bb
y = np.logspace(wmin, wmax, 500) * l_bb
for i, ax in enumerate(axes):
    ax.plot(
        np.logspace(wmin, wmax, 500),
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
fig = pypython.plot.finish_figure(fig, hspace=0, wspace=0)
fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}-sed.pdf", dpi=300)

plt.show()
