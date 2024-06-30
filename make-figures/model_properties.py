import os

import numpy as np
import pypython
from matplotlib import pyplot as plt

wind = pypython.Wind(
    "input", "../simulations/2d/full", distance_units="rg", version="88"
)
spectrum = pypython.Spectrum("input", "../simulations/2d/full")

parameters = [
    "t_e",
    "ne",
    "t_r",
    "ip",
    "H_i01f",
    "He_i02f",
]

names = [
    r"Electron temperature",
    r"Hydrogen density",
    r"Radiation temperature",
    r"Ionization parameter",
    r"H~\textsc{i} fraction",
    r"He~\textsc{ii} fraction",
]

units = [" [K]", " [K]", "", r" [cm$^{-3}$]", "", ""]
fig, ax = plt.subplots(3, 2, figsize=(12, 15), subplot_kw={"projection": "polar"})

for n, (parameter, name, unit) in enumerate(zip(parameters, names, units)):
    i, j = np.unravel_index(n, (3, 2))

    if parameter == "He_i02f":
        fig, ax = pypython.wind.plot.wind(
            wind, parameter, fig=fig, ax=ax, i=i, j=j, vmin=-8, vmax=-2
        )
    elif parameter == "ne":
        fig, ax = pypython.wind.plot.wind(
            wind,
            np.log10(wind.get("H_i01d") + wind.get("H_i02d")),
            fig=fig,
            ax=ax,
            i=i,
            j=j,
        )
    else:
        fig, ax = pypython.wind.plot.wind(wind, parameter, fig=fig, ax=ax, i=i, j=j)

    ax[i, j].set_title(r"$\log_{10}($" + name + r"$)$" + unit)
    ax[i, j].set_rlim(np.log10(45), np.log10(wind.x_axis_coords[-3]))
    ax[i, j].set_ylabel("")
    ax[i, j].set_thetamin(5)
    ax[i, j].set_thetamax(88)

n_points = 15
for inclination, linestyle in zip(
    spectrum.inclinations[1:-1][::2] + ("87",), ["-", "--", ":", "D-", "s-"]
):
    theta = np.ones(n_points) * float(inclination)
    r = np.logspace(np.log10(45), np.log10(wind["r"].max()), n_points)
    ax[-1, 0].plot(
        np.deg2rad(theta),
        np.log10(r),
        f"k{linestyle}",
        label=r"$i = " + f"{inclination}" + r"^{\circ}$",
        alpha=0.7,
    )

ax[-1, 0].legend(loc="upper right", fontsize=10).set_zorder(0)

fig.text(0.020, 0.5, r"$\log_{10}(r/r_{g})$", rotation="vertical")
fig = pypython.plot.finish_figure(fig)
fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
plt.show()
