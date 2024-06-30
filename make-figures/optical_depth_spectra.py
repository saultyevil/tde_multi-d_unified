import os

import numpy as np
import pypython
from matplotlib import pyplot as plt

xmin, xmax = 3e14, 1e18

alpha = 0.75
spec_full = pypython.Spectrum("input", "../simulations/2d/full", log=True)
spec_full_no_es = pypython.Spectrum("input_no_es", "../simulations/2d/full")
spec_full.set("spec_tau")

fig, ax = plt.subplots(figsize=(7, 6))
ax = pypython.plot.set_axes_scales(ax, "loglog")


for n, inclination in enumerate(spec_full.inclinations[::-1]):
    y = spec_full[inclination]
    if np.count_nonzero(y) == 0:
        continue
    ax.plot(
        spec_full["Freq."],
        y,
        label=str(inclination) + r"$^{\circ}$",
        alpha=1,
        color=f"C{n}",
    )

for n, inclination in enumerate(spec_full.inclinations[::-1]):
    x, y = pypython.get_xy_subset(
        spec_full_no_es["Freq."], spec_full_no_es[inclination], 1e14, 2e17
    )
    if np.count_nonzero(y) == 0:
        continue
    ax.plot(x, y, linestyle="--", alpha=0.5, color=f"C{n}", zorder=0)

ax.legend(loc="lower right", fontsize=12)
ax.set_xlim(xmin, xmax)
ax.set_ylim(1e-5, 1e9)
ax.set_xlabel("Rest-frame frequency [Hz]")
ax.set_ylabel("Continuum optical depth $\tau$")

edges = [
    [r"O \textsc{vii}", 16],
    [r"C \textsc{v}", 32],
    [r"O \textsc{vi }", 93],
    [r"He \textsc{ii}", 227],
    [r"He \textsc{i}", 504],
    ["Lyman", 912],
    ["Balmer", 3646],
]

for i in range(len(edges)):
    edges[i][1] = pypython.constants.C / (edges[i][1] * pypython.constants.ANGSTROM)

ax = pypython.spectrum.plot.add_line_ids(ax, edges, "none", offset=0)

fig = pypython.plot.finish_figure(fig)
fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
plt.show()
