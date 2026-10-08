"""
Efield_viewer.py

Plot a 3-D surface of the normalized electric field inside the dielectric
sample + vial using functions from:

    cavity_sample_field_neighbor_normalized.py

Put this file in the SAME DIRECTORY as
cavity_sample_field_neighbor_normalized.py.

The horizontal axes are local cavity coordinates x and z, with
x = 0, z = 0 at the sample/cavity center.

The vertical axis is normalized E or |E|.
"""

import numpy as np
import matplotlib.pyplot as plt

from cavity_sample_field_neighbor_normalized import (
    Cavity,
    Sample,
    solve_loaded_cavity,
    loaded_field,
)


# ============================================================
# USER SETTINGS
# ============================================================

# Cavity dimensions [mm]
A_MM = 6.0 * 25.4
ELL_MM = 18.0 * 25.4

# TE_10n mode number
MODE_N = 3

# Sample
SAMPLE_RADIUS_MM = 6.78/2
EPS_SAMPLE = 3.08

# Optional vial
# Set VIAL_OUTER_RADIUS_MM = None for no vial.
VIAL_OUTER_RADIUS_MM = 9.08/2
EPS_VIAL = 2.54

# Eigenmode basis size
M_MAX = 21
P_MAX = 41

# Plot resolution
N_GRID = 201

# False -> plot signed E
# True  -> plot |E|
PLOT_MAGNITUDE = True

# 3-D viewing angle
ELEVATION_DEG = 30
AZIMUTH_DEG = -55


# ============================================================
# BUILD PROBLEM
# ============================================================

def make_problem():
    """Construct the Cavity and Sample objects."""

    cavity = Cavity(
        a_mm=A_MM,
        ell_mm=ELL_MM,
    )

    sample = Sample(
        radius_mm=SAMPLE_RADIUS_MM,
        eps_r=EPS_SAMPLE,
        vial_outer_radius_mm=VIAL_OUTER_RADIUS_MM,
        eps_vial=EPS_VIAL,
    )

    return cavity, sample


# ============================================================
# SOLVE CAVITY
# ============================================================

def solve_problem():
    """Solve the loaded-cavity eigenmode once."""

    cavity, sample = make_problem()

    solution = solve_loaded_cavity(
        cavity=cavity,
        sample=sample,
        n=MODE_N,
        m_max=M_MAX,
        p_max=P_MAX,
    )

    return solution


# ============================================================
# SAMPLE + VIAL GRID
# ============================================================

def make_sample_vial_grid(solution, n_grid=N_GRID):
    """
    Build an x-z grid covering the whole sample + vial region.

    If no vial is present, the domain extends only to the sample radius.

    Returns
    -------
    X, Z : 2-D arrays
        Local coordinates [mm].

    E : 2-D array
        Normalized loaded electric field.

    sample_mask : 2-D bool array
        True inside the sample.

    vial_mask : 2-D bool array
        True only in the vial shell.

    outer_radius_mm : float
        Radius of the plotted circular domain.
    """

    if solution.sample.vial_outer_radius_mm is None:
        outer_radius_mm = float(
            solution.sample.radius_mm
        )
    else:
        outer_radius_mm = float(
            solution.sample.vial_outer_radius_mm
        )

    coords = np.linspace(
        -outer_radius_mm,
        +outer_radius_mm,
        int(n_grid),
    )

    X, Z = np.meshgrid(
        coords,
        coords,
    )

    R = np.sqrt(
        X**2 + Z**2
    )

    E = np.asarray(
        loaded_field(
            solution,
            X,
            Z,
        ),
        dtype=float,
    )

    sample_mask = (
        R <= solution.sample.radius_mm
    )

    if solution.sample.vial_outer_radius_mm is None:

        vial_mask = np.zeros_like(
            sample_mask,
            dtype=bool,
        )

    else:

        vial_mask = (
            (R > solution.sample.radius_mm)
            &
            (R <= solution.sample.vial_outer_radius_mm)
        )

    outer_mask = (
        R <= outer_radius_mm
    )

    E = np.where(
        outer_mask,
        E,
        np.nan,
    )

    return (
        X,
        Z,
        E,
        sample_mask,
        vial_mask,
        outer_radius_mm,
    )


# ============================================================
# FIELD ON CIRCULAR BOUNDARY
# ============================================================

def circular_boundary_curve(
    solution,
    radius_mm,
    n_points=500,
    magnitude=False,
):
    """
    Evaluate the field around a circular boundary.

    Used to draw the sample and vial boundaries on the 3-D surface.
    """

    theta = np.linspace(
        0.0,
        2.0 * np.pi,
        int(n_points),
    )

    x = (
        radius_mm
        * np.cos(theta)
    )

    z = (
        radius_mm
        * np.sin(theta)
    )

    E = np.asarray(
        loaded_field(
            solution,
            x,
            z,
        ),
        dtype=float,
    )

    if magnitude:
        E = np.abs(E)

    return x, z, E


# ============================================================
# 3-D PLOT
# ============================================================

def plot_3d_field(
    solution,
    n_grid=N_GRID,
    magnitude=PLOT_MAGNITUDE,
    elevation_deg=ELEVATION_DEG,
    azimuth_deg=AZIMUTH_DEG,
):
    """
    Plot a 3-D surface of the electric field inside the sample + vial.
    """

    (
        X,
        Z,
        E,
        sample_mask,
        vial_mask,
        outer_radius_mm,
    ) = make_sample_vial_grid(
        solution,
        n_grid=n_grid,
    )

    if magnitude:

        E_plot = np.abs(E)
        vertical_label = "Normalized |E|"

        title = (
            f"TE_10{solution.n} electric-field magnitude "
            "inside sample + vial"
        )

    else:

        E_plot = E
        vertical_label = "Normalized E"

        title = (
            f"TE_10{solution.n} electric field "
            "inside sample + vial"
        )

    fig = plt.figure(
        figsize=(10, 8)
    )

    ax = fig.add_subplot(
        111,
        projection="3d",
    )

    surface = ax.plot_surface(
        X,
        Z,
        E_plot,
        cmap="viridis",
        linewidth=0,
        antialiased=True,
    )

    # Draw sample boundary directly on the field surface.
    x_sample, z_sample, E_sample = circular_boundary_curve(
        solution,
        radius_mm=solution.sample.radius_mm,
        magnitude=magnitude,
    )

    ax.plot(
        x_sample,
        z_sample,
        E_sample,
        linewidth=2.0,
        label="Sample boundary",
    )

    # Draw outer vial boundary, if present.
    if solution.sample.vial_outer_radius_mm is not None:

        x_vial, z_vial, E_vial = circular_boundary_curve(
            solution,
            radius_mm=solution.sample.vial_outer_radius_mm,
            magnitude=magnitude,
        )

        ax.plot(
            x_vial,
            z_vial,
            E_vial,
            linewidth=2.0,
            linestyle="--",
            label="Vial outer boundary",
        )

    cbar = fig.colorbar(
        surface,
        ax=ax,
        shrink=0.70,
        pad=0.10,
    )

    cbar.set_label(
        vertical_label
    )

    ax.set_xlabel(
        "x from sample center [mm]"
    )

    ax.set_ylabel(
        "z from sample center [mm]"
    )

    ax.set_zlabel(
        vertical_label
    )

    ax.set_title(
        title
    )

    ax.set_xlim(
        -outer_radius_mm,
        +outer_radius_mm,
    )

    ax.set_ylim(
        -outer_radius_mm,
        +outer_radius_mm,
    )

    ax.set_box_aspect(
        (
            1.0,
            1.0,
            0.75,
        )
    )

    ax.view_init(
        elev=elevation_deg,
        azim=azimuth_deg,
    )

    ax.legend()

    fig.tight_layout()

    return fig, ax


# ============================================================
# PRINT DIAGNOSTICS
# ============================================================

def print_solution_summary(solution):
    """Print a few useful properties of the solved field."""

    print()
    print(
        f"Mode: TE_10{solution.n}"
    )

    print(
        "Empty-cavity frequency:",
        solution.f_empty_Hz / 1e9,
        "GHz",
    )

    print(
        "Loaded-cavity frequency:",
        solution.f_loaded_Hz / 1e9,
        "GHz",
    )

    print(
        "Frequency shift:",
        (
            solution.f_empty_Hz
            -
            solution.f_loaded_Hz
        ) / 1e6,
        "MHz",
    )

    E_center = float(
        np.asarray(
            loaded_field(
                solution,
                0.0,
                0.0,
            )
        )
    )

    print(
        "Normalized center field:",
        E_center,
    )

    if solution.n == 1:

        print(
            "Normalization: target coefficient "
            "(TE_101)"
        )

    else:

        z_neighbor = (
            solution.cavity.ell_mm
            /
            solution.n
        )

        E_plus = float(
            np.asarray(
                loaded_field(
                    solution,
                    0.0,
                    +z_neighbor,
                )
            )
        )

        E_minus = float(
            np.asarray(
                loaded_field(
                    solution,
                    0.0,
                    -z_neighbor,
                )
            )
        )

        mean_abs = 0.5 * (
            abs(E_plus)
            +
            abs(E_minus)
        )

        print(
            "Neighboring antinode z:",
            z_neighbor,
            "mm",
        )

        print(
            "Mean neighboring |E|:",
            mean_abs,
        )


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":

    solution = solve_problem()

    print_solution_summary(
        solution
    )

    plot_3d_field(
        solution,
        n_grid=N_GRID,
        magnitude=PLOT_MAGNITUDE,
        elevation_deg=ELEVATION_DEG,
        azimuth_deg=AZIMUTH_DEG,
    )

    plt.show()
