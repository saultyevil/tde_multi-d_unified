import os

import numpy as np
import pypython
from matplotlib import pyplot as plt

PLOT_ALPHA = 0.75
SMOOTH = 4
wind = pypython.Wind(
    "input", "../simulations/2d/full", distance_units="rg", version="88"
)

fig, ax = plt.subplots(1, 2, figsize=(12, 5.5), sharex=True, sharey=True)

for n, theta_cell in enumerate([20, 108]):
    for r_cell in [2, 30, 60, 100, 148, 168, 200, 220]:
        cell_spec = wind.cell_spec[int(r_cell), int(theta_cell)]
        if cell_spec is None:
            print(f"Skipping theta = {theta_cell} and r = {r_cell}")
            continue

        r = wind.r[int(r_cell), int(theta_cell)]
        ax[n].plot(
            cell_spec["Freq."],
            pypython.smooth_array(cell_spec["Freq."] * cell_spec["Flux"], SMOOTH),
            label=f"${np.log10(r):.1f}$",
            alpha=PLOT_ALPHA,
        )

    ax[n].set_xlim(1e14, 2e17)
    ax[n].set_ylim(1e8, 1e18)
    ax[n].text(
        0.92,
        0.92,
        r"$" + f"{wind.theta[int(r_cell), int(theta_cell)]:.0f}" + r"^{\circ}$",
        transform=ax[n].transAxes,
    )
    ax[n] = pypython.plot.set_axes_scales(ax[n], "loglog")

leg = ax[0].legend(loc="upper left", fontsize=10, ncol=2)
leg.set_zorder(0)
fig.text(
    0.02,
    0.5,
    r"$\nu J_{\nu}$ [ergs s$^{-1}$ cm$^{-3}$ sr$^{-1}$]",
    rotation="vertical",
    ha="center",
    va="center",
)
fig.text(0.5, 0.02, "Rest-frame frequency [Hz]", ha="center", va="center")
fig = pypython.plot.finish_figure(fig, hspace=0, wspace=0)
fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)

plt.show()
