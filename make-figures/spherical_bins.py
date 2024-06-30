import os

import numpy as np
import pypython
from matplotlib import pyplot as plt

ALPHA = 0.8
radius, theta, rho, t_gas, t_rad, v_r, v_theta = np.genfromtxt(
    "../other-data/quantities2.txt", unpack=True
)
nx, nz = 128, 64

grid = {
    "r": radius.reshape(nx, nz),
    "theta": theta.reshape(nx, nz),
    "rho": rho.reshape(nx, nz),
    "t_gas": t_gas.reshape(nx, nz),
    "t_rad": t_rad.reshape(nx, nz),
    "v_r": v_r.reshape(nx, nz),
    "v_theta": v_theta.reshape(nx, nz),
}
bins = [
    [67.5, 87.4, "Bin 1"],
    [45, 67.5, "Bin 2"],
    [22.5, 45, "Bin 3"],
    [5.7, 22.5, "Bin 4"],
]

fig, ax = plt.subplots(figsize=(8.5, 6), subplot_kw={"projection": "polar"})

im = ax.pcolormesh(grid["theta"], grid["r"], grid["rho"], linewidth=0, rasterized=True)
ax.set_ylabel(r"$r / r_{g}$")
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_thetamin(0)
ax.set_thetamax(90)
ax.set_rlabel_position(90)
ax.set_rlim(0, 2000)

cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r"$\log_{10}(\rho)$" + " [g cm$^{-3}$]")

for limits in bins:
    theta1, theta2, name = limits[0], limits[1], limits[2]
    x_coords = np.logspace(0, np.log10(3000), 50)
    theta_coords1, theta_coords2 = (
        np.ones_like(x_coords) * theta1,
        np.ones_like(x_coords) * theta2,
    )
    ax.plot(np.deg2rad(theta_coords1), x_coords, "k--", alpha=ALPHA)
    ax.plot(np.deg2rad(theta_coords2), x_coords, "k--", alpha=ALPHA)
    theta_cen = 0.5 * (theta2 + theta1)
    ax.text(
        np.deg2rad(theta_cen),
        1250,
        name,
        color="k",
        rotation=90 - theta_cen,
        va="center",
        ha="center",
    )

fig = pypython.plot.finish_figure(fig)
fig.savefig(f"{os.path.splitext(os.path.basename(__file__))[0]}.pdf", dpi=300)
plt.show()
