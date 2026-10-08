import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.interpolate import griddata

from scipy.special import j1
from scipy.linalg import eigh


# ============================================================
# Disk Fourier integral
# ============================================================

def disk_fourier_integral(Q, R):
    """
    Integral over a centered disk:

        F(Q) = ∫_disk cos(qx X + qz Z) dA

    where Q = sqrt(qx^2 + qz^2).

    All lengths are in mm, so F has units mm^2.
    """

    if np.isclose(Q, 0.0):
        return np.pi * R**2

    return (
        2.0 * np.pi * R
        * j1(Q * R)
        / Q
    )


# ============================================================
# Overlap between two rectangular-cavity modes over a
# centered circular disk
# ============================================================

def disk_mode_overlap(
    m1, p1,
    m2, p2,
    a, ell, R,
):
    """
    Calculate

        S_ij(R) = ∫_disk phi_i phi_j dA

    for normalized cavity eigenfunctions

        phi_mp =
            2/sqrt(a*ell)
            sin(m*pi*x/a)
            sin(p*pi*z/ell)

    and a circular disk centered at x=a/2, z=ell/2.

    This implementation assumes odd m and p, which is the
    symmetry sector relevant to centered TE_10n modes with odd n.
    """

    # Center signs
    s_m1 = np.sin(m1 * np.pi / 2.0)
    s_m2 = np.sin(m2 * np.pi / 2.0)

    s_p1 = np.sin(p1 * np.pi / 2.0)
    s_p2 = np.sin(p2 * np.pi / 2.0)

    sign = (
        s_m1 * s_m2
        * s_p1 * s_p2
    )

    total = 0.0

    # Product-to-sum terms:
    #
    # cos(m1 X) cos(m2 X)
    # cos(p1 Z) cos(p2 Z)
    #
    # generate m1 +/- m2 and p1 +/- p2.
    for sx in [-1, 1]:

        qx = (
            (m1 + sx*m2)
            * np.pi / a
        )

        for sz in [-1, 1]:

            qz = (
                (p1 + sz*p2)
                * np.pi / ell
            )

            Q = np.sqrt(
                qx**2 + qz**2
            )

            total += disk_fourier_integral(
                Q,
                R,
            )

    return (
        sign
        * total
        / (a * ell)
    )


# ============================================================
# Build the odd-odd cavity basis
# ============================================================

def make_basis(
    a,
    ell,
    m_max,
    p_max,
):
    """
    Build the odd-m, odd-p symmetry sector.

    This is sufficient for a centered circular sample/vial
    perturbing an odd TE_10n mode.
    """

    modes = []

    for m in range(1, m_max + 1, 2):

        for p in range(1, p_max + 1, 2):

            k2 = (
                (m * np.pi / a)**2
                +
                (p * np.pi / ell)**2
            )

            modes.append(
                (m, p, k2)
            )

    return modes


# ============================================================
# Build overlap matrix S(R)
# ============================================================

def build_disk_overlap_matrix(
    modes,
    a,
    ell,
    R,
):
    """
    Matrix

        S_ij(R) = ∫_disk phi_i phi_j dA
    """

    N = len(modes)

    S = np.zeros(
        (N, N),
        dtype=float,
    )

    for i in range(N):

        m1, p1, _ = modes[i]

        for j in range(i, N):

            m2, p2, _ = modes[j]

            value = disk_mode_overlap(
                m1, p1,
                m2, p2,
                a, ell, R,
            )

            S[i, j] = value
            S[j, i] = value

    return S


# ============================================================
# Solve loaded cavity
# ============================================================

def solve_loaded_cavity(
    a_mm,
    ell_mm,
    target_n,
    sample_radius_mm,
    eps_sample,
    vial_outer_radius_mm=None,
    eps_vial=1.0,
    m_max=21,
    p_max=61,
    normalization="target_coefficient",
):
    """
    Solve the 2-D loaded rectangular cavity eigenproblem.

    All lengths are in mm.

    Parameters
    ----------
    a_mm
        Cavity width in x [mm].

    ell_mm
        Cavity length in z [mm].

    target_n
        Target TE_10n mode.

        Example:
            target_n = 7  -> TE_107

    sample_radius_mm
        Sample radius [mm].

    eps_sample
        Relative permittivity of sample.

    vial_outer_radius_mm
        Outer vial radius [mm].

        If None, a solid sample with no vial is assumed.

    eps_vial
        Relative permittivity of vial.

    m_max, p_max
        Maximum rectangular-cavity basis indices.

    normalization
        "target_coefficient":
            coefficient of the original TE_10n basis mode
            is set equal to its empty-cavity value.

            This is useful for directly examining field
            redistribution.

        "energy":
            electric energy integral

                ∫ epsilon |E|^2 dA

            is set equal to that of an empty TE_10n mode
            having center amplitude 1.

    Returns
    -------
    Dictionary containing eigenvalue, coefficients, basis,
    frequency, and helper field function.
    """

    if target_n % 2 == 0:
        raise ValueError(
            "This centered-sample implementation assumes "
            "an odd target n."
        )

    a = float(a_mm)
    ell = float(ell_mm)
    Rs = float(sample_radius_mm)

    # --------------------------------------------------------
    # Basis
    # --------------------------------------------------------

    modes = make_basis(
        a,
        ell,
        m_max,
        p_max,
    )

    Nbasis = len(modes)

    #print(
    #    f"Number of cavity basis functions: {Nbasis}"
    #)

    # --------------------------------------------------------
    # K matrix
    #
    # -laplacian phi_i = k_i^2 phi_i
    # --------------------------------------------------------

    k2_empty = np.array(
        [mode[2] for mode in modes]
    )

    K = np.diag(k2_empty)

    # --------------------------------------------------------
    # Permittivity matrix M
    #
    # M_ij = ∫ epsilon(x,z) phi_i phi_j dA
    # --------------------------------------------------------

    I = np.eye(Nbasis)

    S_sample = build_disk_overlap_matrix(
        modes,
        a,
        ell,
        Rs,
    )

    if vial_outer_radius_mm is None:

        # Solid dielectric cylinder
        M = (
            I
            +
            (eps_sample - 1.0)
            * S_sample
        )

    else:

        Rh = float(
            vial_outer_radius_mm
        )

        if Rh <= Rs:
            raise ValueError(
                "vial_outer_radius_mm must be "
                "larger than sample_radius_mm."
            )

        S_holder = build_disk_overlap_matrix(
            modes,
            a,
            ell,
            Rh,
        )

        # Fill entire outer disk with vial, then replace
        # the inner disk by sample material.
        M = (
            I
            +
            (eps_vial - 1.0)
            * S_holder
            +
            (eps_sample - eps_vial)
            * S_sample
        )

    # --------------------------------------------------------
    # Solve
    #
    #       K c = k^2 M c
    #
    # eigh returns M-orthonormal eigenvectors.
    # --------------------------------------------------------

    eigvals, eigvecs = eigh(
        K,
        M,
    )

    # --------------------------------------------------------
    # Locate target TE_10n basis function
    # --------------------------------------------------------

    target_index = None

    for i, (m, p, _) in enumerate(modes):

        if m == 1 and p == target_n:

            target_index = i
            break

    if target_index is None:
        raise ValueError(
            "Target TE_10n mode is outside the "
            "chosen m_max/p_max basis."
        )

    # --------------------------------------------------------
    # Identify loaded eigenmode having the strongest
    # projection onto empty TE_10n.
    # --------------------------------------------------------

    target_components = np.abs(
        eigvecs[target_index, :]
    )

    loaded_index = np.argmax(
        target_components
    )

    k2_loaded = eigvals[
        loaded_index
    ]

    k_loaded = np.sqrt(
        k2_loaded
    )

    coeff = eigvecs[
        :,
        loaded_index
    ].copy()

    # --------------------------------------------------------
    # Normalize field
    # --------------------------------------------------------

    if normalization == "target_coefficient":

        # Make the TE_10n modal coefficient equal to 1.
        coeff /= coeff[
            target_index
        ]

        # Our phi_mp basis has center magnitude
        #
        #     2 / sqrt(a ell)
        #
        # for odd modes.
        #
        # Multiply so an empty TE_10n field would have
        # center amplitude +1.
        center_sign = np.sin(
            target_n * np.pi / 2.0
        )

        coeff *= (
            center_sign
            * np.sqrt(a * ell)
            / 2.0
        )

    elif normalization == "energy":

        # eigh gives c^T M c = 1.
        #
        # An empty TE_10n field normalized to center
        # amplitude 1 has
        #
        # ∫ |E0|^2 dA = a*ell/4.
        #
        coeff *= np.sqrt(
            a * ell / 4.0
        )

        # Fix arbitrary sign
        empty_sign = np.sin(
            target_n * np.pi / 2.0
        )

        center_value = 0.0

        for c, (m, p, _) in zip(
            coeff,
            modes,
        ):

            center_value += (
                c
                * 2.0
                / np.sqrt(a * ell)
                * np.sin(m*np.pi/2)
                * np.sin(p*np.pi/2)
            )

        if np.sign(center_value) != np.sign(
            empty_sign
        ):
            coeff *= -1.0

    else:

        raise ValueError(
            "normalization must be "
            "'target_coefficient' or 'energy'."
        )

    # --------------------------------------------------------
    # Empty target wave number
    # --------------------------------------------------------

    k_target_empty = np.sqrt(
        (np.pi / a)**2
        +
        (target_n*np.pi / ell)**2
    )

    # --------------------------------------------------------
    # Frequency
    #
    # k is currently mm^-1.
    # Convert to m^-1.
    # --------------------------------------------------------

    c_light = 299792458.0

    f_loaded = (
        c_light
        * (k_loaded * 1000.0)
        / (2.0*np.pi)
    )

    f_empty = (
        c_light
        * (k_target_empty * 1000.0)
        / (2.0*np.pi)
    )

    # --------------------------------------------------------
    # Field evaluator
    #
    # X,Z are LOCAL coordinates:
    #
    #       X = 0, Z = 0
    #
    # is sample/cavity center.
    # --------------------------------------------------------

    def E_field(X, Z):

        X = np.asarray(X)
        Z = np.asarray(Z)

        # Convert local -> rectangular cavity coordinates
        x_global = X + a/2.0
        z_global = Z + ell/2.0

        E = np.zeros(
            np.broadcast(
                x_global,
                z_global,
            ).shape,
            dtype=float,
        )

        for c, (m, p, _) in zip(
            coeff,
            modes,
        ):

            phi = (
                2.0
                / np.sqrt(a * ell)
                * np.sin(
                    m*np.pi*x_global/a
                )
                * np.sin(
                    p*np.pi*z_global/ell
                )
            )

            E += c * phi

        return E

    # --------------------------------------------------------
    # Empty TE_10n field normalized to center amplitude 1
    # --------------------------------------------------------

    def E_empty(X, Z):

        return (
            np.cos(
                np.pi * X / a
            )
            *
            np.cos(
                target_n
                * np.pi
                * Z
                / ell
            )
        )

    return {
        "modes": modes,
        "coefficients": coeff,

        "K": K,
        "M": M,

        "loaded_mode_index":
            loaded_index,

        "target_basis_index":
            target_index,

        "k_loaded_mm_inv":
            k_loaded,

        "k_empty_mm_inv":
            k_target_empty,

        "f_loaded_Hz":
            f_loaded,

        "f_empty_Hz":
            f_empty,

        "frequency_shift_Hz":
            f_empty - f_loaded,

        "fractional_shift":
            (f_empty - f_loaded)
            / f_empty,

        "E":
            E_field,

        "E_empty":
            E_empty,
    }




# ============================================================
# Read XFdtd planar sensor CSV
# ============================================================

def read_xfdtd_plane(
    filename,
    x_column="vertex_X (mm)",
    z_column="vertex_Z (mm)",
    field_column="Total E y (V/m)",
    time_column="Time (ns)",
    time_ns=None,
    x_center_mm=0.0,
    z_center_mm=0.0,
    average_duplicates=True,
):
    """
    Read an XFdtd planar-sensor CSV.

    Parameters
    ----------
    filename : str
        Path to CSV.

    x_column, z_column : str
        The two XF coordinate columns corresponding to the
        x and z coordinates of the theoretical model.

        For example, if XF's X-Y plane corresponds to your
        theoretical x-z plane, use:

            x_column = "vertex_X (mm)"
            z_column = "vertex_Y (mm)"

    field_column : str
        Electric-field component to compare.

        Examples:
            "Total E y (V/m)"
            "Total E z (V/m)"

    time_column : str
        Time column in the XF export.

    time_ns : float or None
        If the CSV contains multiple times, select the time
        nearest this value.

        If None and only one time is present, that time is used.
        If None and multiple times are present, the last time is used.

    x_center_mm, z_center_mm : float
        XF coordinates of the sample/cavity center.

        Returned coordinates are

            x_local = XF_x - x_center_mm
            z_local = XF_z - z_center_mm

        so that the theoretical sample center is (0,0).

    average_duplicates : bool
        Average field values if XF contains duplicate coordinate
        points.

    Returns
    -------
    dataframe with columns:

        x_mm
        z_mm
        E_xf
    """

    df = pd.read_csv(filename)

    # Remove accidental spaces from headers
    df.columns = df.columns.str.strip()

    required = [
        x_column,
        z_column,
        field_column,
    ]

    for col in required:
        if col not in df.columns:
            raise KeyError(
                f"Column '{col}' not found.\n"
                f"Available columns are:\n{list(df.columns)}"
            )

    # --------------------------------------------------------
    # Select one time if time information exists
    # --------------------------------------------------------

    if time_column in df.columns:

        times = np.sort(
            df[time_column].dropna().unique()
        )

        if len(times) > 1:

            if time_ns is None:
                selected_time = times[-1]
            else:
                selected_time = times[
                    np.argmin(
                        np.abs(times - time_ns)
                    )
                ]

            print(
                f"Using XF time = "
                f"{selected_time:.8g} ns"
            )

            df = df[
                np.isclose(
                    df[time_column],
                    selected_time,
                )
            ].copy()

        elif len(times) == 1:

            print(
                f"Using XF time = "
                f"{times[0]:.8g} ns"
            )

    # --------------------------------------------------------
    # Keep only needed data
    # --------------------------------------------------------

    out = pd.DataFrame({
        "x_mm":
            pd.to_numeric(
                df[x_column],
                errors="coerce",
            )
            - x_center_mm,

        "z_mm":
            pd.to_numeric(
                df[z_column],
                errors="coerce",
            )
            - z_center_mm,

        "E_xf":
            pd.to_numeric(
                df[field_column],
                errors="coerce",
            ),
    })

    out = out.dropna()

    # --------------------------------------------------------
    # Average duplicate vertices if present
    # --------------------------------------------------------

    if average_duplicates:

        out = (
            out
            .groupby(
                ["x_mm", "z_mm"],
                as_index=False,
            )
            ["E_xf"]
            .mean()
        )

    print(
        f"Loaded {len(out)} unique field points."
    )

    print(
        f"x range: "
        f"{out['x_mm'].min():.3f} to "
        f"{out['x_mm'].max():.3f} mm"
    )

    print(
        f"z range: "
        f"{out['z_mm'].min():.3f} to "
        f"{out['z_mm'].max():.3f} mm"
    )

    return out


# ============================================================
# Plot scattered XF field on regular grid
# ============================================================

def plot_field(
    x,
    z,
    E,
    title="Electric field",
    nx=400,
    nz=800,
    cmap="RdBu_r",
    symmetric=True,
):
    """
    Plot scattered field samples as a 2-D interpolated map.
    """

    x = np.asarray(x)
    z = np.asarray(z)
    E = np.asarray(E)

    xi = np.linspace(
        x.min(),
        x.max(),
        nx,
    )

    zi = np.linspace(
        z.min(),
        z.max(),
        nz,
    )

    XX, ZZ = np.meshgrid(
        xi,
        zi,
    )

    EE = griddata(
        (x, z),
        E,
        (XX, ZZ),
        method="linear",
    )

    fig, ax = plt.subplots(
        figsize=(7, 10)
    )

    if symmetric:

        vmax = np.nanmax(
            np.abs(EE)
        )

        vmin = -vmax

    else:

        vmin = np.nanmin(EE)
        vmax = np.nanmax(EE)

    pcm = ax.pcolormesh(
        XX,
        ZZ,
        EE,
        shading="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )

    fig.colorbar(
        pcm,
        ax=ax,
        label="Electric field",
    )

    ax.set_xlabel(
        "x from sample center [mm]"
    )

    ax.set_ylabel(
        "z from sample center [mm]"
    )

    ax.set_title(title)

    ax.set_aspect(
        "equal"
    )

    plt.tight_layout()

    return fig, ax


# ============================================================
# Compare XF field directly to theoretical eigenmode solution
# ============================================================

def compare_xfdtd_to_model(
    xf_data,
    E_model,
    a_mm,
    ell_mm,
    normalization="least_squares",
    wall_margin_mm=5.0,
    energy_weighted=False,
    percent_cutoff=0.01,
    sample_radius_mm=None,
    eps_sample=1.0,
    vial_outer_radius_mm=None,
    eps_vial=1.0,
):
    """
    Compare XFdtd field to the eigenmode model.

    The model is evaluated at the exact XF sample locations.

    If energy_weighted=False, the plotted residual is the ordinary
    signed point-by-point percent residual

        100 * (E_XF - E_model) / |E_XF|.

    For that ordinary percent residual, points with

        |E_XF| < percent_cutoff * max_interior(|E_XF|)

    are masked to avoid divergences at field nodes.

    If energy_weighted=True, the plotted residual is instead weighted
    by the local relative electric-energy density

        w = eps_r |E_XF|^2 / max_interior(eps_r |E_XF|^2),

    clipped to 0 <= w <= 1, and

        weighted residual = w * 100 * (E_XF - E_model) / |E_XF|.

    This weighted form naturally goes to zero at electric-field nodes,
    so the percent_cutoff is not applied to the weighted residual.

    The interior reference region excludes wall_margin_mm from each wall,
    preventing strong feed fields near the walls from setting the model
    normalization or the energy-density reference.
    """

    data = xf_data.copy()

    x = data["x_mm"].to_numpy()
    z = data["z_mm"].to_numpy()
    E_xf = data["E_xf"].to_numpy()

    # ========================================================
    # Evaluate model at EXACT XF coordinates
    # ========================================================

    E_theory_raw = np.asarray(E_model(x, z))
    E_theory_raw = np.real_if_close(E_theory_raw)

    # ========================================================
    # Interior region used for normalization/reference values
    # ========================================================

    interior_mask = (
        (np.abs(x) <= a_mm / 2.0 - wall_margin_mm)
        &
        (np.abs(z) <= ell_mm / 2.0 - wall_margin_mm)
    )

    if np.sum(interior_mask) == 0:
        raise ValueError(
            "wall_margin_mm is too large; no interior points remain."
        )

    # ========================================================
    # Normalize model amplitude using only the interior region
    # ========================================================

    if normalization == "none":
        scale = 1.0
        E_theory = E_theory_raw.copy()

    elif normalization == "peak":
        xf_peak_interior = np.max(np.abs(E_xf[interior_mask]))
        model_peak_interior = np.max(np.abs(E_theory_raw[interior_mask]))

        if model_peak_interior == 0:
            raise ValueError(
                "The theoretical field is zero throughout the interior region."
            )

        scale = xf_peak_interior / model_peak_interior
        E_theory = scale * E_theory_raw

    elif normalization == "least_squares":
        xf_fit = E_xf[interior_mask]
        model_fit = E_theory_raw[interior_mask]

        denominator = np.sum(np.abs(model_fit) ** 2)

        if denominator == 0:
            raise ValueError(
                "The theoretical field is zero throughout the fitting region."
            )

        scale = (
            np.sum(xf_fit * np.conj(model_fit))
            / denominator
        )
        scale = np.real_if_close(scale)
        E_theory = scale * E_theory_raw

    else:
        raise ValueError(
            "normalization must be 'none', 'peak', or 'least_squares'."
        )

    # ========================================================
    # Signed field residual
    # Positive means XF > model
    # ========================================================

    residual = E_xf - E_theory

    # ========================================================
    # Local relative permittivity
    # ========================================================

    eps_local = np.ones_like(E_xf, dtype=float)
    r = np.sqrt(x**2 + z**2)

    # Fill vial region first.
    if vial_outer_radius_mm is not None:
        vial_mask = r <= vial_outer_radius_mm
        eps_local[vial_mask] = eps_vial

    # Replace the inner vial region with sample material.
    if sample_radius_mm is not None:
        sample_mask = r <= sample_radius_mm
        eps_local[sample_mask] = eps_sample

    # ========================================================
    # Electric-energy-density weighting
    # u_E is proportional to eps_r |E|^2.
    # ========================================================

    energy_density = eps_local * np.abs(E_xf) ** 2

    energy_reference = np.max(energy_density[interior_mask])

    if energy_reference <= 0:
        raise ValueError(
            "The interior XF electric-energy density is zero."
        )

    energy_weight = energy_density / energy_reference

    # A strong feed field outside the interior region should not receive
    # a weight larger than the cavity-field reference.
    energy_weight = np.clip(energy_weight, 0.0, 1.0)

    # ========================================================
    # Ordinary signed point-by-point percent residual
    # ========================================================

    xf_peak_interior = np.max(np.abs(E_xf[interior_mask]))
    field_threshold = percent_cutoff * xf_peak_interior

    percent_residual = np.full(E_xf.shape, np.nan, dtype=float)

    valid_percent = np.abs(E_xf) >= field_threshold

    percent_residual[valid_percent] = (
        100.0
        * residual[valid_percent]
        / np.abs(E_xf[valid_percent])
    )

    # ========================================================
    # Energy-weighted signed percent residual
    #
    # This is evaluated safely so an exact node gives zero rather
    # than inf/NaN. Near a node, energy_weight ~ |E_XF|^2, which
    # suppresses the 1/|E_XF| behavior of the ordinary percent error.
    # ========================================================

    safe_fractional_residual = np.zeros_like(E_xf, dtype=float)
    nonzero = np.abs(E_xf) > 0.0

    safe_fractional_residual[nonzero] = (
        residual[nonzero]
        / np.abs(E_xf[nonzero])
    )

    energy_weighted_percent = (
        100.0
        * energy_weight
        * safe_fractional_residual
    )

    # ========================================================
    # Select which residual the plotting function should use
    # ========================================================

    if energy_weighted:
        residual_for_plot = energy_weighted_percent
    else:
        residual_for_plot = percent_residual

    # ========================================================
    # Store results
    # ========================================================

    data["E_model_raw"] = E_theory_raw
    data["E_model"] = E_theory
    data["residual"] = residual
    data["percent_residual"] = percent_residual
    data["eps_r"] = eps_local
    data["energy_density_relative"] = energy_density
    data["energy_weight"] = energy_weight
    data["energy_weighted_percent_residual"] = energy_weighted_percent
    data["residual_for_plot"] = residual_for_plot
    data["interior_mask"] = interior_mask

    # Save plotting metadata.
    data.attrs["energy_weighted"] = bool(energy_weighted)
    data.attrs["percent_cutoff"] = float(percent_cutoff)
    data.attrs["wall_margin_mm"] = float(wall_margin_mm)

    # ========================================================
    # Useful statistics
    # ========================================================

    print()
    print("Model amplitude scale =", scale)
    print("Interior XF peak |E| =", xf_peak_interior, "V/m")
    print("Wall exclusion distance =", wall_margin_mm, "mm")
    print("Maximum interior relative energy density =", energy_reference)

    if energy_weighted:
        print("Residual plotted: energy-weighted percent residual")
        print(
            "Energy-weighted residual range =",
            np.nanmin(energy_weighted_percent),
            "to",
            np.nanmax(energy_weighted_percent),
            "%",
        )
    else:
        print("Residual plotted: ordinary percent residual")
        print(
            "Ordinary-percent cutoff =",
            100.0 * percent_cutoff,
            "% of the interior XF peak field",
        )

    return data


# ============================================================
# Plot XF, model, and residual
# ============================================================

def plot_comparison(
    comparison,
    nx=400,
    nz=800,
    xlim=6*25.4,
    ylim=18*25.4,
    percent_limit=None,
):
    """
    Produce three figures:

        1. XF electric field
        2. Eigenmode-model electric field
        3. Residual map

    If compare_xfdtd_to_model(..., energy_weighted=True) was used,
    the third figure is the energy-weighted signed percent residual.
    Otherwise it is the ordinary signed percent residual.

    percent_limit can be set to a number such as 10 to force the
    residual colorbar to span -10% to +10%. If None, the finite
    residual range is used symmetrically about zero.
    """

    x = comparison["x_mm"].to_numpy()
    z = comparison["z_mm"].to_numpy()

    E_xf = comparison["E_xf"].to_numpy()
    E_model = comparison["E_model"].to_numpy()
    residual_for_plot = comparison["residual_for_plot"].to_numpy()

    energy_weighted = comparison.attrs.get("energy_weighted", False)

    # Common field scale for XF and model.
    common_peak = max(
        np.max(np.abs(E_xf)),
        np.max(np.abs(E_model)),
    )

    # --------------------------------------------------------
    # Common interpolation grid
    # --------------------------------------------------------

    xi = np.linspace(x.min(), x.max(), nx)
    zi = np.linspace(z.min(), z.max(), nz)
    XX, ZZ = np.meshgrid(xi, zi)

    XF_grid = griddata(
        (x, z),
        E_xf,
        (XX, ZZ),
        method="linear",
    )

    model_grid = griddata(
        (x, z),
        E_model,
        (XX, ZZ),
        method="linear",
    )

    # Only interpolate finite residual values.
    valid = np.isfinite(residual_for_plot)

    if np.sum(valid) < 3:
        raise ValueError(
            "Too few valid points remain to interpolate the residual."
        )

    residual_grid = griddata(
        (x[valid], z[valid]),
        residual_for_plot[valid],
        (XX, ZZ),
        method="linear",
    )

    # --------------------------------------------------------
    # XF field
    # --------------------------------------------------------

    fig1, ax1 = plt.subplots(figsize=(7, 10))

    p1 = ax1.pcolormesh(
        XX,
        ZZ,
        XF_grid,
        shading="auto",
        cmap="RdBu_r",
        vmin=-common_peak,
        vmax=common_peak,
    )

    fig1.colorbar(p1, ax=ax1, label="E [V/m]")

    ax1.set_title("XFdtd electric field")
    ax1.set_xlabel("x from sample center [mm]")
    ax1.set_ylabel("z from sample center [mm]")
    ax1.set_aspect("equal")
    fig1.tight_layout()

    # --------------------------------------------------------
    # Model field
    # --------------------------------------------------------

    fig2, ax2 = plt.subplots(figsize=(7, 10))

    p2 = ax2.pcolormesh(
        XX,
        ZZ,
        model_grid,
        shading="auto",
        cmap="RdBu_r",
        vmin=-common_peak,
        vmax=common_peak,
    )

    fig2.colorbar(p2, ax=ax2, label="E [V/m]")

    ax2.set_title("Eigenmode model electric field")
    ax2.set_xlabel("x from sample center [mm]")
    ax2.set_ylabel("z from sample center [mm]")
    ax2.set_aspect("equal")
    fig2.tight_layout()

    # --------------------------------------------------------
    # Residual
    # --------------------------------------------------------

    if percent_limit is None:
        percent_peak = np.nanmax(np.abs(residual_grid))
    else:
        percent_peak = float(percent_limit)

    if percent_peak == 0:
        percent_peak = 1.0

    fig3, ax3 = plt.subplots(figsize=(7, 10))

    p3 = ax3.pcolormesh(
        XX,
        ZZ,
        residual_grid,
        shading="auto",
        cmap="RdBu_r",
        vmin=-percent_peak,
        vmax=percent_peak,
    )

    if energy_weighted:
        colorbar_label = "Energy-weighted (XF - model) / |XF| [%]"
        title = "Energy-weighted electric-field percent residual"
    else:
        colorbar_label = "(XF - model) / |XF| [%]"
        title = "Point-by-point electric-field percent residual"

    fig3.colorbar(
        p3,
        ax=ax3,
        label=colorbar_label,
    )

    ax3.set_ylim(-ylim, ylim)
    ax3.set_xlim(-xlim, xlim)
    ax3.set_title(title)
    ax3.set_xlabel("x from sample center [mm]")
    ax3.set_ylabel("z from sample center [mm]")
    ax3.set_aspect("equal")
    fig3.tight_layout()

    plt.show()




a = 58.16
ell = 304.8
n = 9
Rs = 9/2
Rh = 9.1/2
eps_sample = 2.54
eps_vial = 1.0


result = solve_loaded_cavity(
    a_mm=a,
    ell_mm=ell,
    target_n=n,
    sample_radius_mm=Rs,
    eps_sample=eps_sample,
    vial_outer_radius_mm=Rh,
    eps_vial=eps_vial,
    m_max=21,
    p_max=61,
    normalization=
        "target_coefficient",
)

## "../../Downloads/UHF_1to3_PET_TE107_xz_centerfield.csv"

xf = read_xfdtd_plane(
    "../../Downloads/WR229_Rexolite_TE109_xz_centerfield.csv",

    ## In the WR229 XF model the long axis is along z, width is x
    ## In the UHF1to3 model the long axis is along y, width is x
    x_column="vertex_X (mm)",
    z_column="vertex_Z (mm)",

    field_column="Total E y (V/m)",

    x_center_mm=0.0,
    z_center_mm=210.0,
    time_ns=158.615
)

comparison = compare_xfdtd_to_model(
    xf,
    E_model=result["E"],
    a_mm=a,
    ell_mm=ell,
    normalization="least_squares",
    wall_margin_mm=5.0,
    energy_weighted=True,
    percent_cutoff=0.01,
    sample_radius_mm=Rs,
    eps_sample=eps_sample,
    vial_outer_radius_mm=Rh,
    eps_vial=eps_vial,
)

plot_comparison(
    comparison,
    xlim=58.16/2,
    ylim=304.8/2/9,
    percent_limit=5,
)