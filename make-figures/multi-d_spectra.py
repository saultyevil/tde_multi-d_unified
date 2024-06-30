import os

import pypython
from matplotlib import pyplot as plt

alpha = 0.75

model2d_reduced = pypython.Spectrum(
    "../simulations/2d/reduced/input.pf",
    distance=100 * 1e6,
    smooth=50,
    log=True,
)

spec_bl = pypython.Spectrum(
    "../simulations/2d/created/input.pf",
    distance=100 * 1e6,
    smooth=500,
    log=True,
)
spec_bl.convert_flux_to_luminosity()

smooth = 50

models_benchmark = [  # Models with the 2D SED
    pypython.Spectrum(
        "../simulations/1d/Bin-1/input.pf",
        log=True,
        smooth=smooth,
        distance=100 * 1e6,
    ),
    pypython.Spectrum(
        "../simulations/1d/Bin-2/input.pf",
        log=True,
        smooth=smooth,
        distance=100 * 1e6,
    ),
    pypython.Spectrum(
        "../simulations/1d/Bin-3/input.pf",
        log=True,
        smooth=smooth,
        distance=100 * 1e6,
    ),
    pypython.Spectrum(
        "../simulations/1d/Bin-4/input.pf",
        log=True,
        smooth=smooth,
        distance=100 * 1e6,
    ),
]

inclinations = [
    "77",
    "56",
    "34",
    "14",
]

names = ["Bin 1", "Bin 2", "Bin 3", "Bin 4"]

fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
axes = axes.flatten()

for n, model in enumerate(models_benchmark):
    # 1D model
    model.convert_flux_to_luminosity()
    y = model["45"] * model["Lambda"]
    axes[n].plot(
        model["Lambda"],
        y,
        color="C0",
        label="1D",
        alpha=alpha,
    )
    # 2D model
    model2d_reduced.convert_flux_to_luminosity()
    y = model2d_reduced[inclinations[n]] * model2d_reduced["Lambda"]
    axes[n].plot(
        model2d_reduced["Lambda"],
        y,
        color="C1",
        label="2D",
        alpha=alpha,
    )

    # Names
    axes[n].text(0.03, 0.93, names[n], transform=axes[n].transAxes)
    axes[n].text(
        0.03, 0.86, f"{inclinations[n]}" + r"$^{\circ}$", transform=axes[n].transAxes
    )
    axes[n] = pypython.plot.set_axes_scales(axes[n], "loglog")

for ax in axes:
    ax.plot(
        spec_bl["Lambda"],
        spec_bl["Created"] * spec_bl["Lambda"],
        color="k",
        linestyle="--",
        alpha=alpha,
        linewidth=3,
        zorder=0,
    )

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

fig = pypython.plot.finish_figure(fig, wspace=0, hspace=0)
fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
plt.show()
