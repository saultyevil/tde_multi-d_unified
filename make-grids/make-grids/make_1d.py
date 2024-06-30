import copy

import make_2d
import numpy as np
import pypython
from astropy.table import Table

make_2d.OUTPUT_RADIAL_VELOCTY = True
make_2d.OUTPUT_TEMPERATURE = True
make_2d.INTERPOLATION_LEVEL = "cubic"
original, regrid = make_2d.regrid_model(inner_boundary=1.23391)
grid = make_2d.setup_grid_boundaries(regrid)
r_g = pypython.physics.blackhole.gravitational_radius(5e6)

print(f"R_g = {r_g:e} cm -> 1000 R_g = {1000*r_g:e} cm")

# N. Roth models
# the_1d_models = [
#     [67.5, 87.4, 3.221062e12, 4.0e47, "Dai.Bin1"],
#     [45.0, 67.5, 3.549527e12, 3.0e46, "Dai.Bin2"],
#     [22.5, 45.0, 4.310358e12, 3.0e46, "Dai.Bin3"],
#     [5.70, 22.5, 7.004391e12, 3.0e46, "Dai.Bin4"],
# ]

# Python comparisons
the_1d_models = [
    [67.5, 87.4, 3.221062e12, "Bin-1"],
    [45.0, 67.5, 3.549527e12, "Bin-2"],
    [22.5, 45.0, 4.310358e12, "Bin-3"],
    [5.70, 22.5, 7.004391e12, "Bin-4"],
]

nx, nz = np.max(grid["i"]) + 1, np.max(grid["j"]) + 1
theta_points = grid["theta"][: nz - 1]
r_points = grid["r"][::nz]
max_r = r_points.max()

new_models = []

for theta1, theta2, r_in, name in the_1d_models:
    print(theta1, theta2, r_in / r_g, name)

    theta_idx1 = pypython.get_array_index(theta_points, theta1)
    theta_idx2 = pypython.get_array_index(theta_points, theta2)

    r = []
    v_r = []
    rho = []
    t_r = []

    rho = np.zeros_like(r_points)
    v_r = np.zeros_like(r_points)
    t_r = np.zeros_like(r_points)

    # For each cell, take volume weighted theta average
    for n in range(1, nx - 1):  # loop over each radius
        r_idx = n * nz
        r_min = grid["r"][r_idx]
        r_max = grid["r"][((n + 1) * nz)]

        if r_min < r_in:
            continue

        total_volume = 0

        # sweep over theta cells
        for i in range(r_idx + theta_idx1, r_idx + theta_idx2):
            theta_min = np.deg2rad(grid["theta"][i])
            theta_max = np.deg2rad(grid["theta"][i + 1])

            cell_volume = (
                (2.0 / 3.0)
                * np.pi
                * (r_max**3 - r_min**3)
                * (np.cos(theta_min) - np.cos(theta_max))
                * 2
                * np.pi
            )
            total_volume += cell_volume

            # ignore inflow cells
            if grid["v_r"][i] < 0:
                grid["v_r"][i] = 0
                # don't need to zero these, we'll do it later
                # grid["rho"][i] = 0
                # grid["t_r"][i] = 0

            rho[n] += grid["rho"][i] * cell_volume
            v_r[n] += grid["v_r"][i] * cell_volume
            t_r[n] += grid["t_r"][i] * cell_volume

        rho[n] /= total_volume
        v_r[n] /= total_volume
        t_r[n] /= total_volume

    # find where first non-zero v_r is...
    first_non_zero = np.nonzero(v_r)[0][0]

    # then create the table excluding all the zero v_r cells
    table = Table()
    table["i"] = np.arange(0, len(r_points[first_non_zero:]), 1)
    table["r"] = r_points[first_non_zero:]
    table["v_r"] = v_r[first_non_zero:]
    table["rho"] = rho[first_non_zero:]
    table["t_r"] = t_r[first_non_zero:]
    table["t_e"] = 0.9 * t_r[first_non_zero:]

    table["r"].format = "%1.6e"
    table["v_r"].format = "%1.6e"
    table["rho"].format = "%1.6e"
    table["t_e"].format = "%1.6e"
    table["t_r"].format = "%1.6e"

    # Set up ghost cells, this does the first and the the last two cells
    for i in range(3):
        table[-i]["rho"] = 0
        table[-i]["t_e"] = 0
        table[-i]["t_r"] = 0

    # The CODE WITH NO NAME doesn't really like increasing velocity in the
    # ghost cells
    for i in range(1):
        table[-(i + 1)]["v_r"] = table[-2]["v_r"]

    # take a copy, because we're going to modify the table now
    new_models.append(copy.deepcopy(table))

    # remove everything past 1000 rg, to stay consistent with N. Roth
    table = table[table["r"] < 1000 * r_g]
    # Re-index
    table["i"] = np.arange(0, len(table))

    table.write(
        f"output/{name}.txt", format="ascii.fixed_width_two_line", overwrite=True
    )

# Now I want to create models where the inner edge has been truncated, like the
# 2D model

inner_edge = 30 * r_g

for model, name in zip(new_models, ["Bin-1", "Bin-2", "Bin-3", "Bin-4"]):
    # extract cells which are beyond the inner edge
    new_table = model[model["r"] > inner_edge]
    # Re-index what's left
    new_table["i"] = np.arange(0, len(new_table))
    # Set new first cell as ghost
    new_table[0]["rho"] = 0
    new_table[0]["t_e"] = 0
    new_table[0]["t_r"] = 0
    # Save to disk
    new_table.write(
        f"output/2d_bc.{name}.txt",
        format="ascii.fixed_width_two_line",
        overwrite=True,
    )
