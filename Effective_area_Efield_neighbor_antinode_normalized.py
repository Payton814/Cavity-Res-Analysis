import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.interpolate import griddata

from scipy.special import j1
from scipy.linalg import eigh

#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable


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
# Reference-region normalization helpers
# ============================================================

def make_reference_region_mask(
    x,
    z,
    a_mm,
    ell_mm,
    wall_margin_mm=5.0,
    reference_r_min_mm=None,
    reference_r_max_mm=None,
):
    """
    Build the mask used to establish the common cavity amplitude.

    The reference region is an annulus centered on the sample/cavity center,
    restricted to points that are also at least wall_margin_mm from the
    rectangular cavity walls.

    This is intended to use field points that are:
        - outside the dielectric sample/vial,
        - away from the feed and conducting walls,
        - still close enough to the central antinode that the same cavity mode
          dominates.
    """

    x = np.asarray(x)
    z = np.asarray(z)

    interior_mask = (
        (np.abs(x) <= a_mm / 2.0 - wall_margin_mm)
        &
        (np.abs(z) <= ell_mm / 2.0 - wall_margin_mm)
    )

    r = np.sqrt(x**2 + z**2)

    reference_mask = interior_mask.copy()

    if reference_r_min_mm is not None:
        reference_mask &= r >= float(reference_r_min_mm)

    if reference_r_max_mm is not None:
        reference_mask &= r <= float(reference_r_max_mm)

    return reference_mask


def compute_reference_scale(
    xf_data,
    E_reference,
    a_mm,
    ell_mm,
    wall_margin_mm=5.0,
    reference_r_min_mm=None,
    reference_r_max_mm=None,
):
    """
    Determine ONE global field-amplitude scale from a reference region.

    The scale alpha minimizes

        sum_reference (E_XF - alpha E_reference)^2

    so

        alpha = sum(E_XF E_reference*) / sum(|E_reference|^2).

    For the loaded-vs-empty comparison, E_reference should normally be the
    theoretical EMPTY-cavity TE_10n field.  The resulting alpha can then be
    applied to BOTH the empty-cavity field and the loaded eigenmode field.

    This is important: it prevents each model from independently fitting away
    the field enhancement inside the dielectric.
    """

    x = xf_data["x_mm"].to_numpy()
    z = xf_data["z_mm"].to_numpy()
    E_xf = xf_data["E_xf"].to_numpy()

    E_ref = np.asarray(E_reference(x, z))
    E_ref = np.real_if_close(E_ref)

    reference_mask = make_reference_region_mask(
        x=x,
        z=z,
        a_mm=a_mm,
        ell_mm=ell_mm,
        wall_margin_mm=wall_margin_mm,
        reference_r_min_mm=reference_r_min_mm,
        reference_r_max_mm=reference_r_max_mm,
    )

    finite = (
        np.isfinite(E_xf)
        & np.isfinite(E_ref)
    )

    reference_mask &= finite

    if np.sum(reference_mask) < 3:
        raise ValueError(
            "Too few points remain in the reference normalization region. "
            "Adjust reference_r_min_mm, reference_r_max_mm, or wall_margin_mm."
        )

    xf_fit = E_xf[reference_mask]
    ref_fit = E_ref[reference_mask]

    denominator = np.sum(np.abs(ref_fit)**2)

    if denominator == 0:
        raise ValueError(
            "The theoretical reference field is zero throughout the reference region."
        )

    scale = (
        np.sum(xf_fit * np.conj(ref_fit))
        / denominator
    )

    scale = np.real_if_close(scale)

    print()
    print("Reference-region normalization")
    print("------------------------------")
    print("Reference points =", np.sum(reference_mask))
    print("Reference r_min =", reference_r_min_mm, "mm")
    print("Reference r_max =", reference_r_max_mm, "mm")
    print("Wall exclusion =", wall_margin_mm, "mm")
    print("Shared amplitude scale =", scale)

    return scale, reference_mask



# ============================================================
# Neighboring-antinode normalization helpers
# ============================================================

def make_neighbor_antinode_mask(
    x,
    z,
    a_mm,
    ell_mm,
    target_n,
    wall_margin_mm=5.0,
    antinode_x_halfwidth_mm=None,
    antinode_z_halfwidth_mm=None,
    neighbor_order=1,
):
    """
    Build a normalization mask around the two longitudinal antinodes
    neighboring the sample at z=0.

    For the local empty-cavity field

        E_0(x,z) = cos(pi x/a) cos(n pi z/ell),

    the central sample is at the antinode z=0 and the neighboring
    longitudinal antinodes are

        z = +/- neighbor_order * ell / n.

    The mask includes rectangular windows centered on those two antinodes,
    while also excluding points too close to the conducting cavity walls.

    Parameters
    ----------
    antinode_x_halfwidth_mm : float or None
        Half-width of each reference window in x.  If None, use a/4.

    antinode_z_halfwidth_mm : float or None
        Half-width of each reference window in z.  If None, use 15% of the
        antinode spacing ell/n.

    neighbor_order : int
        1 uses the nearest antinodes at +/- ell/n.
        2 would use the next pair at +/- 2 ell/n, etc.
    """

    x = np.asarray(x)
    z = np.asarray(z)

    n = int(target_n)
    if n <= 0:
        raise ValueError("target_n must be positive.")

    neighbor_order = int(neighbor_order)
    if neighbor_order < 1:
        raise ValueError("neighbor_order must be >= 1.")

    antinode_spacing = float(ell_mm) / n
    z_antinode = neighbor_order * antinode_spacing

    if antinode_x_halfwidth_mm is None:
        antinode_x_halfwidth_mm = float(a_mm) / 4.0

    if antinode_z_halfwidth_mm is None:
        antinode_z_halfwidth_mm = 0.15 * antinode_spacing

    antinode_x_halfwidth_mm = float(antinode_x_halfwidth_mm)
    antinode_z_halfwidth_mm = float(antinode_z_halfwidth_mm)

    if antinode_x_halfwidth_mm <= 0 or antinode_z_halfwidth_mm <= 0:
        raise ValueError("Antinode window half-widths must be positive.")

    interior_mask = (
        (np.abs(x) <= float(a_mm) / 2.0 - wall_margin_mm)
        &
        (np.abs(z) <= float(ell_mm) / 2.0 - wall_margin_mm)
    )

    x_window = np.abs(x) <= antinode_x_halfwidth_mm

    plus_window = (
        x_window
        & (np.abs(z - z_antinode) <= antinode_z_halfwidth_mm)
        & interior_mask
    )

    minus_window = (
        x_window
        & (np.abs(z + z_antinode) <= antinode_z_halfwidth_mm)
        & interior_mask
    )

    combined_mask = plus_window | minus_window

    return combined_mask, plus_window, minus_window, z_antinode


def compute_neighbor_antinode_scale(
    xf_data,
    E_reference,
    a_mm,
    ell_mm,
    target_n,
    wall_margin_mm=5.0,
    antinode_x_halfwidth_mm=None,
    antinode_z_halfwidth_mm=None,
    neighbor_order=1,
):
    """
    Determine ONE global amplitude scale from the two neighboring antinodes.

    The reference field should normally be the theoretical EMPTY-cavity
    TE_10n field.  We fit one scalar alpha that minimizes

        sum_ref |E_XF - alpha E_reference|^2

    over windows around the two neighboring antinodes.  The theoretical
    field retains its sign, so the opposite sign of adjacent TE-mode lobes is
    handled automatically.

    The same alpha should then be applied unchanged to BOTH the empty-cavity
    field and the loaded eigenmode field.
    """

    x = xf_data["x_mm"].to_numpy()
    z = xf_data["z_mm"].to_numpy()
    E_xf = xf_data["E_xf"].to_numpy()

    E_ref = np.asarray(E_reference(x, z))
    E_ref = np.real_if_close(E_ref)

    (
        reference_mask,
        plus_mask,
        minus_mask,
        z_antinode,
    ) = make_neighbor_antinode_mask(
        x=x,
        z=z,
        a_mm=a_mm,
        ell_mm=ell_mm,
        target_n=target_n,
        wall_margin_mm=wall_margin_mm,
        antinode_x_halfwidth_mm=antinode_x_halfwidth_mm,
        antinode_z_halfwidth_mm=antinode_z_halfwidth_mm,
        neighbor_order=neighbor_order,
    )

    finite = np.isfinite(E_xf) & np.isfinite(E_ref)
    reference_mask &= finite
    plus_mask &= finite
    minus_mask &= finite

    if np.sum(reference_mask) < 3:
        raise ValueError(
            "Too few XF points fall inside the neighboring-antinode "
            "normalization windows. Increase the antinode window sizes."
        )

    xf_fit = E_xf[reference_mask]
    ref_fit = E_ref[reference_mask]

    denominator = np.sum(np.abs(ref_fit) ** 2)
    if denominator == 0:
        raise ValueError(
            "The theoretical empty-cavity field is zero throughout the "
            "neighboring-antinode reference windows."
        )

    scale = np.sum(xf_fit * np.conj(ref_fit)) / denominator
    scale = np.real_if_close(scale)

    def side_scale(mask):
        if np.sum(mask) < 1:
            return np.nan
        e_xf_side = E_xf[mask]
        e_ref_side = E_ref[mask]
        den = np.sum(np.abs(e_ref_side) ** 2)
        if den == 0:
            return np.nan
        return np.real_if_close(
            np.sum(e_xf_side * np.conj(e_ref_side)) / den
        )

    scale_plus = side_scale(plus_mask)
    scale_minus = side_scale(minus_mask)

    # Reconstruct the actual defaults used, for printing/metadata.
    if antinode_x_halfwidth_mm is None:
        antinode_x_halfwidth_mm = float(a_mm) / 4.0
    if antinode_z_halfwidth_mm is None:
        antinode_z_halfwidth_mm = 0.15 * (float(ell_mm) / int(target_n))

    print()
    print("Neighboring-antinode normalization")
    print("----------------------------------")
    print("Neighbor antinode z positions =", -z_antinode, "and", z_antinode, "mm")
    print("x half-width =", antinode_x_halfwidth_mm, "mm")
    print("z half-width =", antinode_z_halfwidth_mm, "mm")
    print("Wall exclusion =", wall_margin_mm, "mm")
    print("-z antinode points =", np.sum(minus_mask))
    print("+z antinode points =", np.sum(plus_mask))
    print("Combined reference points =", np.sum(reference_mask))
    print("Scale from -z antinode =", scale_minus)
    print("Scale from +z antinode =", scale_plus)
    print("Shared amplitude scale =", scale)

    if np.isfinite(scale_plus) and np.isfinite(scale_minus):
        denom_mean = 0.5 * (abs(scale_plus) + abs(scale_minus))
        if denom_mean > 0:
            side_difference = 100.0 * abs(scale_plus - scale_minus) / denom_mean
            print("Antinode-to-antinode scale difference =", side_difference, "%")

    metadata = {
        "neighbor_order": int(neighbor_order),
        "z_antinode_mm": float(z_antinode),
        "antinode_x_halfwidth_mm": float(antinode_x_halfwidth_mm),
        "antinode_z_halfwidth_mm": float(antinode_z_halfwidth_mm),
        "scale_plus": scale_plus,
        "scale_minus": scale_minus,
    }

    return scale, reference_mask, plus_mask, minus_mask, metadata


# ============================================================
# Compare XF field directly to a theoretical model
# ============================================================

def compare_xfdtd_to_model(
    xf_data,
    E_model,
    a_mm,
    ell_mm,
    normalization="least_squares",
    wall_margin_mm=5.0,
    reference_r_min_mm=None,
    reference_r_max_mm=None,
    scale_override=None,
    residual_mode="magnitude",
    energy_weighted=False,
    percent_cutoff=0.01,
    sample_radius_mm=None,
    eps_sample=1.0,
    vial_outer_radius_mm=None,
    eps_vial=1.0,
    model_label="Eigenmode model",
):
    """
    Compare XFdtd field to a theoretical model at the exact XF sample points.

    normalization options
    ---------------------
    "none"
        No fitted amplitude scale.

    "peak"
        Match the model peak to the XF peak using the interior cavity region.

    "least_squares"
        Fit one amplitude scale over the entire interior cavity region.
        This is useful for pure shape comparison, but it can partially fit away
        a real dielectric field enhancement.

    "reference_region"
        Fit the amplitude only in the user-defined annular reference region.
        This is the preferred option when asking whether the dielectric changes
        the field amplitude near the sample.

    scale_override
        If supplied, use this amplitude scale directly and do not refit it.
        This lets the loaded and empty models use exactly the SAME scale.

    residual_mode options
    ---------------------
    "magnitude"  [recommended]

        residual = |E_XF| - |E_model|

        Positive residual means XF has the larger FIELD MAGNITUDE.
        Negative residual means the theoretical model has the larger magnitude.

    "signed_field"

        residual = E_XF - E_model

        This preserves the instantaneous field sign, but a sign change between
        TE-mode lobes can make "larger amplitude" difficult to interpret.

    energy_weighted
    ---------------
    If True, the plotted percent residual is weighted by

        w = eps_r |E_XF|^2 / max_interior(eps_r |E_XF|^2),

    clipped to 0 <= w <= 1.

    Therefore electric-field nodes naturally contribute almost zero weight.
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
    # Interior region
    # ========================================================

    interior_mask = (
        (np.abs(x) <= a_mm / 2.0 - wall_margin_mm)
        &
        (np.abs(z) <= ell_mm / 2.0 - wall_margin_mm)
    )

    finite = (
        np.isfinite(E_xf)
        & np.isfinite(E_theory_raw)
    )

    interior_mask &= finite

    if np.sum(interior_mask) == 0:
        raise ValueError(
            "wall_margin_mm is too large; no interior points remain."
        )

    # Reference annulus, if requested/provided.
    reference_mask = make_reference_region_mask(
        x=x,
        z=z,
        a_mm=a_mm,
        ell_mm=ell_mm,
        wall_margin_mm=wall_margin_mm,
        reference_r_min_mm=reference_r_min_mm,
        reference_r_max_mm=reference_r_max_mm,
    )
    reference_mask &= finite

    # ========================================================
    # Model amplitude normalization
    # ========================================================

    if scale_override is not None:

        scale = np.real_if_close(scale_override)
        normalization_used = "shared_reference_scale"

    elif normalization == "none":

        scale = 1.0
        normalization_used = "none"

    elif normalization == "peak":

        xf_peak_interior = np.max(np.abs(E_xf[interior_mask]))
        model_peak_interior = np.max(np.abs(E_theory_raw[interior_mask]))

        if model_peak_interior == 0:
            raise ValueError(
                "The theoretical field is zero throughout the interior region."
            )

        scale = xf_peak_interior / model_peak_interior
        normalization_used = "peak"

    elif normalization == "least_squares":

        xf_fit = E_xf[interior_mask]
        model_fit = E_theory_raw[interior_mask]

        denominator = np.sum(np.abs(model_fit)**2)

        if denominator == 0:
            raise ValueError(
                "The theoretical field is zero throughout the fitting region."
            )

        scale = (
            np.sum(xf_fit * np.conj(model_fit))
            / denominator
        )
        scale = np.real_if_close(scale)
        normalization_used = "least_squares"

    elif normalization == "reference_region":

        if reference_r_min_mm is None or reference_r_max_mm is None:
            raise ValueError(
                "normalization='reference_region' requires both "
                "reference_r_min_mm and reference_r_max_mm."
            )

        if np.sum(reference_mask) < 3:
            raise ValueError(
                "Too few points remain in the reference normalization region."
            )

        xf_fit = E_xf[reference_mask]
        model_fit = E_theory_raw[reference_mask]

        denominator = np.sum(np.abs(model_fit)**2)

        if denominator == 0:
            raise ValueError(
                "The theoretical field is zero throughout the reference region."
            )

        scale = (
            np.sum(xf_fit * np.conj(model_fit))
            / denominator
        )
        scale = np.real_if_close(scale)
        normalization_used = "reference_region"

    else:
        raise ValueError(
            "normalization must be 'none', 'peak', 'least_squares', "
            "or 'reference_region'."
        )

    E_theory = scale * E_theory_raw

    # ========================================================
    # Residual definition
    # ========================================================

    if residual_mode == "magnitude":

        # This is the physically useful comparison for the question:
        # does the dielectric increase the field AMPLITUDE?
        residual = np.abs(E_xf) - np.abs(E_theory)

    elif residual_mode == "signed_field":

        residual = E_xf - E_theory

    else:
        raise ValueError(
            "residual_mode must be 'magnitude' or 'signed_field'."
        )

    # ========================================================
    # Local relative permittivity for electric-energy weighting
    # ========================================================

    eps_local = np.ones_like(E_xf, dtype=float)
    r = np.sqrt(x**2 + z**2)

    if vial_outer_radius_mm is not None:
        vial_mask = r <= vial_outer_radius_mm
        eps_local[vial_mask] = eps_vial

    if sample_radius_mm is not None:
        sample_mask = r <= sample_radius_mm
        eps_local[sample_mask] = eps_sample

    # ========================================================
    # Electric-energy-density weight
    # ========================================================

    energy_density = eps_local * np.abs(E_xf)**2

    energy_reference = np.max(
        energy_density[interior_mask]
    )

    if energy_reference <= 0:
        raise ValueError(
            "The interior XF electric-energy density is zero."
        )

    energy_weight = (
        energy_density
        / energy_reference
    )

    # Do not let a strong feed field outside the interior fitting region
    # acquire a weight greater than the cavity reference.
    energy_weight = np.clip(
        energy_weight,
        0.0,
        1.0,
    )

    # ========================================================
    # Ordinary point-by-point percent residual
    # ========================================================

    xf_peak_interior = np.max(
        np.abs(E_xf[interior_mask])
    )

    field_threshold = (
        percent_cutoff
        * xf_peak_interior
    )

    percent_residual = np.full(
        E_xf.shape,
        np.nan,
        dtype=float,
    )

    valid_percent = (
        np.abs(E_xf) >= field_threshold
    )

    percent_residual[valid_percent] = (
        100.0
        * residual[valid_percent]
        / np.abs(E_xf[valid_percent])
    )

    # ========================================================
    # Energy-weighted percent residual
    # ========================================================

    safe_fractional_residual = np.zeros_like(
        E_xf,
        dtype=float,
    )

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

    if energy_weighted:
        residual_for_plot = energy_weighted_percent
    else:
        residual_for_plot = percent_residual

    # ========================================================
    # Store results and metadata
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
    data["reference_mask"] = reference_mask

    data.attrs["energy_weighted"] = bool(energy_weighted)
    data.attrs["percent_cutoff"] = float(percent_cutoff)
    data.attrs["wall_margin_mm"] = float(wall_margin_mm)
    data.attrs["model_label"] = str(model_label)
    data.attrs["residual_mode"] = str(residual_mode)
    data.attrs["normalization"] = str(normalization_used)
    data.attrs["amplitude_scale"] = scale
    data.attrs["reference_r_min_mm"] = reference_r_min_mm
    data.attrs["reference_r_max_mm"] = reference_r_max_mm

    # ========================================================
    # Useful statistics
    # ========================================================

    print()
    print(model_label)
    print("-" * len(model_label))
    print("Amplitude scale =", scale)
    print("Normalization used =", normalization_used)
    print("Residual mode =", residual_mode)
    print("Interior XF peak |E| =", xf_peak_interior, "V/m")
    print("Wall exclusion distance =", wall_margin_mm, "mm")
    print("Maximum interior relative energy density =", energy_reference)

    if reference_r_min_mm is not None or reference_r_max_mm is not None:
        print("Reference-region points =", np.sum(reference_mask))

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
# Theoretical empty-cavity TE_10n field
# ============================================================

def make_empty_cavity_field(a_mm, ell_mm, target_n):
    """
    Return the theoretical EMPTY rectangular-cavity TE_10n field shape
    in local coordinates, with X=Z=0 at the cavity/sample center.

        E_empty(X,Z) = cos(pi X/a) cos(n pi Z/ell)

    For odd n this has unit magnitude at the center.
    """

    a = float(a_mm)
    ell = float(ell_mm)
    n = int(target_n)

    if n % 2 == 0:
        raise ValueError(
            "This centered empty-cavity helper is intended for odd TE_10n modes."
        )

    def E_empty(X, Z):
        X = np.asarray(X)
        Z = np.asarray(Z)

        return (
            np.cos(np.pi * X / a)
            * np.cos(n * np.pi * Z / ell)
        )

    return E_empty


# ============================================================
# Compare XF directly to theoretical EMPTY cavity field
# ============================================================

def compare_xfdtd_to_empty_cavity(
    xf_data,
    a_mm,
    ell_mm,
    target_n,
    normalization="reference_region",
    wall_margin_mm=5.0,
    reference_r_min_mm=None,
    reference_r_max_mm=None,
    scale_override=None,
    residual_mode="magnitude",
    energy_weighted=False,
    percent_cutoff=0.01,
    sample_radius_mm=None,
    eps_sample=1.0,
    vial_outer_radius_mm=None,
    eps_vial=1.0,
):
    """
    Compare loaded or empty XF data to the theoretical EMPTY-cavity TE_10n field.

    For a LOADED XF simulation, pass the actual sample/vial geometry so the
    energy weighting reflects the real loaded material distribution.
    """

    E_empty = make_empty_cavity_field(
        a_mm=a_mm,
        ell_mm=ell_mm,
        target_n=target_n,
    )

    return compare_xfdtd_to_model(
        xf_data=xf_data,
        E_model=E_empty,
        a_mm=a_mm,
        ell_mm=ell_mm,
        normalization=normalization,
        wall_margin_mm=wall_margin_mm,
        reference_r_min_mm=reference_r_min_mm,
        reference_r_max_mm=reference_r_max_mm,
        scale_override=scale_override,
        residual_mode=residual_mode,
        energy_weighted=energy_weighted,
        percent_cutoff=percent_cutoff,
        sample_radius_mm=sample_radius_mm,
        eps_sample=eps_sample,
        vial_outer_radius_mm=vial_outer_radius_mm,
        eps_vial=eps_vial,
        model_label=f"Theoretical empty-cavity TE_10{target_n}",
    )


# ============================================================
# Preferred loaded-vs-empty comparison using ONE shared scale
# ============================================================

def compare_loaded_and_empty_shared_reference(
    xf_data,
    E_loaded_model,
    a_mm,
    ell_mm,
    target_n,
    reference_r_min_mm,
    reference_r_max_mm,
    wall_margin_mm=5.0,
    residual_mode="magnitude",
    energy_weighted=True,
    percent_cutoff=0.01,
    sample_radius_mm=None,
    eps_sample=1.0,
    vial_outer_radius_mm=None,
    eps_vial=1.0,
):
    """
    Compare the SAME loaded XF field to both:

        1. the loaded-cavity eigenmode model, and
        2. the theoretical empty-cavity TE_10n field,

    using ONE COMMON amplitude scale.

    The common scale is fitted between XF and the EMPTY-cavity field only in
    the specified annular reference region.  That same scale is then applied
    unchanged to the loaded eigenmode model.

    This preserves the physical question of interest: after establishing the
    common cavity amplitude away from the dielectric, is the field magnitude
    inside the sample larger than the empty-cavity prediction?

    IMPORTANT:
    The loaded eigenmode should be generated with

        normalization="target_coefficient"

    so its TE_10n coefficient is on the same amplitude convention as the
    empty-cavity field.
    """

    E_empty = make_empty_cavity_field(
        a_mm=a_mm,
        ell_mm=ell_mm,
        target_n=target_n,
    )

    shared_scale, reference_mask = compute_reference_scale(
        xf_data=xf_data,
        E_reference=E_empty,
        a_mm=a_mm,
        ell_mm=ell_mm,
        wall_margin_mm=wall_margin_mm,
        reference_r_min_mm=reference_r_min_mm,
        reference_r_max_mm=reference_r_max_mm,
    )

    common_kwargs = dict(
        xf_data=xf_data,
        a_mm=a_mm,
        ell_mm=ell_mm,
        normalization="none",
        wall_margin_mm=wall_margin_mm,
        reference_r_min_mm=reference_r_min_mm,
        reference_r_max_mm=reference_r_max_mm,
        scale_override=shared_scale,
        residual_mode=residual_mode,
        energy_weighted=energy_weighted,
        percent_cutoff=percent_cutoff,
        sample_radius_mm=sample_radius_mm,
        eps_sample=eps_sample,
        vial_outer_radius_mm=vial_outer_radius_mm,
        eps_vial=eps_vial,
    )

    comparison_loaded = compare_xfdtd_to_model(
        E_model=E_loaded_model,
        model_label=f"Loaded-cavity TE_10{target_n} eigenmode",
        **common_kwargs,
    )

    comparison_empty = compare_xfdtd_to_model(
        E_model=E_empty,
        model_label=f"Theoretical empty-cavity TE_10{target_n}",
        **common_kwargs,
    )

    comparison_loaded.attrs["shared_reference_scale"] = shared_scale
    comparison_empty.attrs["shared_reference_scale"] = shared_scale

    print()
    print("Shared-reference comparison complete")
    print("------------------------------------")
    print("Both models use amplitude scale =", shared_scale)
    print(
        "Residual sign convention:",
        "positive => |XF| > |model|" if residual_mode == "magnitude"
        else "positive => XF > model",
    )

    return comparison_loaded, comparison_empty



# ============================================================
# Preferred comparison using neighboring-antinode normalization
# ============================================================

def compare_loaded_and_empty_neighbor_antinode(
    xf_data,
    E_loaded_model,
    a_mm,
    ell_mm,
    target_n,
    wall_margin_mm=5.0,
    antinode_x_halfwidth_mm=None,
    antinode_z_halfwidth_mm=None,
    neighbor_order=1,
    residual_mode="magnitude",
    energy_weighted=True,
    percent_cutoff=0.01,
    sample_radius_mm=None,
    eps_sample=1.0,
    vial_outer_radius_mm=None,
    eps_vial=1.0,
):
    """
    Compare the same loaded XF field to both the loaded-cavity eigenmode and
    the theoretical empty-cavity TE_10n field using ONE shared amplitude scale.

    The shared scale is determined only from windows around the two neighboring
    longitudinal antinodes of the EMPTY-cavity TE_10n mode.  This avoids using
    the sample region itself to set the normalization and exploits the fact
    that the dielectric perturbation should be much weaker farther from the
    centered sample.

    IMPORTANT:
    E_loaded_model should come from solve_loaded_cavity(...,
    normalization="target_coefficient") so its TE_10n coefficient uses the
    same amplitude convention as the empty-cavity field.
    """

    E_empty = make_empty_cavity_field(
        a_mm=a_mm,
        ell_mm=ell_mm,
        target_n=target_n,
    )

    (
        shared_scale,
        antinode_mask,
        plus_mask,
        minus_mask,
        antinode_meta,
    ) = compute_neighbor_antinode_scale(
        xf_data=xf_data,
        E_reference=E_empty,
        a_mm=a_mm,
        ell_mm=ell_mm,
        target_n=target_n,
        wall_margin_mm=wall_margin_mm,
        antinode_x_halfwidth_mm=antinode_x_halfwidth_mm,
        antinode_z_halfwidth_mm=antinode_z_halfwidth_mm,
        neighbor_order=neighbor_order,
    )

    common_kwargs = dict(
        xf_data=xf_data,
        a_mm=a_mm,
        ell_mm=ell_mm,
        normalization="none",
        wall_margin_mm=wall_margin_mm,
        reference_r_min_mm=None,
        reference_r_max_mm=None,
        scale_override=shared_scale,
        residual_mode=residual_mode,
        energy_weighted=energy_weighted,
        percent_cutoff=percent_cutoff,
        sample_radius_mm=sample_radius_mm,
        eps_sample=eps_sample,
        vial_outer_radius_mm=vial_outer_radius_mm,
        eps_vial=eps_vial,
    )

    comparison_loaded = compare_xfdtd_to_model(
        E_model=E_loaded_model,
        model_label=f"Loaded-cavity TE_10{target_n} eigenmode",
        **common_kwargs,
    )

    comparison_empty = compare_xfdtd_to_model(
        E_model=E_empty,
        model_label=f"Theoretical empty-cavity TE_10{target_n}",
        **common_kwargs,
    )

    # Save the actual antinode masks so they can be inspected/plotted later.
    comparison_loaded["antinode_reference_mask"] = antinode_mask
    comparison_loaded["antinode_plus_mask"] = plus_mask
    comparison_loaded["antinode_minus_mask"] = minus_mask

    comparison_empty["antinode_reference_mask"] = antinode_mask
    comparison_empty["antinode_plus_mask"] = plus_mask
    comparison_empty["antinode_minus_mask"] = minus_mask

    for comp in (comparison_loaded, comparison_empty):
        comp.attrs["normalization"] = "shared_neighbor_antinode_scale"
        comp.attrs["shared_neighbor_antinode_scale"] = shared_scale
        comp.attrs["neighbor_order"] = antinode_meta["neighbor_order"]
        comp.attrs["neighbor_antinode_z_mm"] = antinode_meta["z_antinode_mm"]
        comp.attrs["antinode_x_halfwidth_mm"] = antinode_meta["antinode_x_halfwidth_mm"]
        comp.attrs["antinode_z_halfwidth_mm"] = antinode_meta["antinode_z_halfwidth_mm"]
        comp.attrs["neighbor_scale_plus"] = antinode_meta["scale_plus"]
        comp.attrs["neighbor_scale_minus"] = antinode_meta["scale_minus"]

    print()
    print("Neighbor-antinode shared-scale comparison complete")
    print("--------------------------------------------------")
    print("Both models use amplitude scale =", shared_scale)
    print(
        "Residual sign convention:",
        "positive => |XF| > |model|"
        if residual_mode == "magnitude"
        else "positive => XF > model",
    )

    return comparison_loaded, comparison_empty


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
        2. Theoretical/model electric field
        3. Residual map

    The residual label is selected automatically from the comparison metadata.
    """

    x = comparison["x_mm"].to_numpy()
    z = comparison["z_mm"].to_numpy()

    E_xf = comparison["E_xf"].to_numpy()
    E_model = comparison["E_model"].to_numpy()
    residual_for_plot = comparison["residual_for_plot"].to_numpy()

    energy_weighted = comparison.attrs.get(
        "energy_weighted",
        False,
    )

    model_label = comparison.attrs.get(
        "model_label",
        "Eigenmode model",
    )

    residual_mode = comparison.attrs.get(
        "residual_mode",
        "magnitude",
    )

    # ========================================================
    # Common field scale for XF and model
    # ========================================================

    common_peak = max(
        np.max(np.abs(E_xf)),
        np.max(np.abs(E_model)),
    )

    # ========================================================
    # Common interpolation grid
    # ========================================================

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

    # ========================================================
    # Residual interpolation
    # ========================================================

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

    # ========================================================
    # 1. XF field
    # ========================================================

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

    # ========================================================
    # 2. Model field
    # ========================================================

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

    ax2.set_title(f"{model_label} electric field")
    ax2.set_xlabel("x from sample center [mm]")
    ax2.set_ylabel("z from sample center [mm]")
    ax2.set_aspect("equal")
    fig2.tight_layout()

    # ========================================================
    # 3. Residual
    # ========================================================

    if percent_limit is None:
        residual_peak = np.nanmax(np.abs(residual_grid))
    else:
        residual_peak = float(percent_limit)

    if residual_peak == 0:
        residual_peak = 1.0

    fig3, ax3 = plt.subplots(figsize=(7, 10))

    p3 = ax3.pcolormesh(
        XX,
        ZZ,
        residual_grid,
        shading="auto",
        cmap="RdBu_r",
        vmin=-residual_peak,
        vmax=residual_peak,
    )

    ax3.set_xlim(-xlim, xlim)
    ax3.set_ylim(-ylim, ylim)
    ax3.set_aspect("equal")

    if residual_mode == "magnitude":
        if energy_weighted:
            residual_title = "Energy-weighted field-magnitude residual"
            residual_label = (
                "Energy-weighted (|XF| - |model|) / |XF| [%]"
            )
        else:
            residual_title = "Field-magnitude percent residual"
            residual_label = "(|XF| - |model|) / |XF| [%]"
    else:
        if energy_weighted:
            residual_title = "Energy-weighted signed-field residual"
            residual_label = "Energy-weighted (XF - model) / |XF| [%]"
        else:
            residual_title = "Signed-field percent residual"
            residual_label = "(XF - model) / |XF| [%]"

    ax3.set_title(residual_title)
    ax3.set_xlabel("x from sample center [mm]")
    ax3.set_ylabel("z from sample center [mm]")

    # Colorbar height follows the actual equal-aspect residual axes.
    divider = make_axes_locatable(ax3)

    cax = divider.append_axes(
        "right",
        size="5%",
        pad=0.12,
    )

    cbar = fig3.colorbar(
        p3,
        cax=cax,
    )

    cbar.set_label(residual_label)

    fig3.tight_layout()
    plt.show()


# ============================================================
# Histogram of residuals inside the sample
# ============================================================

def plot_sample_residual_histograms(
    comparison_loaded,
    comparison_empty,
    sample_radius_mm,
    bins=50,
    residual_limit=None,
    density=False,
):
    """
    Plot the residual distribution INSIDE the dielectric sample for both
    the loaded eigenmode model and the empty-cavity model.

    The values come from residual_for_plot, so the histogram automatically
    uses energy-weighted residuals when energy_weighted=True.

    The histogram counts are NOT energy-weighted a second time.  Each XF cell
    contributes one entry whose residual value has already been energy weighted.
    """

    x_loaded = comparison_loaded["x_mm"].to_numpy()
    z_loaded = comparison_loaded["z_mm"].to_numpy()
    residual_loaded = comparison_loaded["residual_for_plot"].to_numpy()

    r_loaded = np.sqrt(x_loaded**2 + z_loaded**2)

    sample_mask_loaded = (
        (r_loaded <= sample_radius_mm)
        & np.isfinite(residual_loaded)
    )

    loaded_values = residual_loaded[sample_mask_loaded]

    x_empty = comparison_empty["x_mm"].to_numpy()
    z_empty = comparison_empty["z_mm"].to_numpy()
    residual_empty = comparison_empty["residual_for_plot"].to_numpy()

    r_empty = np.sqrt(x_empty**2 + z_empty**2)

    sample_mask_empty = (
        (r_empty <= sample_radius_mm)
        & np.isfinite(residual_empty)
    )

    empty_values = residual_empty[sample_mask_empty]

    if len(loaded_values) == 0:
        raise ValueError(
            "No valid loaded-model residual points were found inside the sample."
        )

    if len(empty_values) == 0:
        raise ValueError(
            "No valid empty-cavity residual points were found inside the sample."
        )

    if residual_limit is None:
        max_abs = max(
            np.max(np.abs(loaded_values)),
            np.max(np.abs(empty_values)),
        )
    else:
        max_abs = float(residual_limit)

    if max_abs == 0:
        max_abs = 1.0

    bin_edges = np.linspace(
        -max_abs,
        max_abs,
        bins + 1,
    )

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.hist(
        empty_values,
        bins=bin_edges,
        histtype="step",
        linewidth=2,
        density=density,
        label="Empty-cavity model",
    )

    ax.hist(
        loaded_values,
        bins=bin_edges,
        histtype="step",
        linewidth=2,
        density=density,
        label="Loaded-cavity eigenmode model",
    )

    ax.axvline(
        0.0,
        linestyle="--",
        linewidth=1,
    )

    residual_mode = comparison_loaded.attrs.get(
        "residual_mode",
        "magnitude",
    )

    energy_weighted = comparison_loaded.attrs.get(
        "energy_weighted",
        False,
    )

    if residual_mode == "magnitude":
        if energy_weighted:
            xlabel = "Energy-weighted (|XF| - |model|) / |XF| [%]"
        else:
            xlabel = "(|XF| - |model|) / |XF| [%]"
    else:
        if energy_weighted:
            xlabel = "Energy-weighted (XF - model) / |XF| [%]"
        else:
            xlabel = "(XF - model) / |XF| [%]"

    ax.set_xlabel(xlabel)

    if density:
        ax.set_ylabel("Probability density")
    else:
        ax.set_ylabel("Number of XF cells")

    ax.set_title("Residual distribution inside sample")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()

    # ========================================================
    # Print statistics
    # ========================================================

    print()
    print("Inside-sample residual statistics")
    print("---------------------------------")
    print("Number of XF cells inside sample:", len(loaded_values))

    def print_stats(label, values):
        print()
        print(label + ":")
        print("  Mean =", np.mean(values), "%")
        print("  Median =", np.median(values), "%")
        print("  RMS =", np.sqrt(np.mean(values**2)), "%")
        print("  Std =", np.std(values), "%")
        print("  Min =", np.min(values), "%")
        print("  Max =", np.max(values), "%")

    print_stats(
        "Loaded-cavity eigenmode model",
        loaded_values,
    )

    print_stats(
        "Empty-cavity model",
        empty_values,
    )

    if residual_mode == "magnitude":
        print()
        print("Sign convention:")
        print("  positive -> |XF| > |model|")
        print("  negative -> |model| > |XF|")

    plt.show()

    return loaded_values, empty_values


# ============================================================
# USER SETTINGS / EXAMPLE
# ============================================================

a = 58.16
ell = 304.8
n = 9
Rs = 9 / 2
Rh = 9.1 / 2
eps_sample = 3.8
eps_vial = 1.0


result = solve_loaded_cavity(
    a_mm=a,
    ell_mm=ell,
    target_n=n,
    sample_radius_mm=Rs,
    eps_sample=eps_sample,
    vial_outer_radius_mm=Rh,
    eps_vial=eps_vial,
    m_max=41,
    p_max=61,
    # IMPORTANT for the shared empty/loaded amplitude convention:
    normalization="target_coefficient",
)


xf = read_xfdtd_plane(
    "../../Downloads/WR229_Quartz_TE109_xz_centerfield.csv",

    # In the WR229 XF model the long axis is along z, width is x.
    # In the UHF1to3 model the long axis is along y, width is x.
    x_column="vertex_X (mm)",
    z_column="vertex_Z (mm)",

    field_column="Total E y (V/m)",

    x_center_mm=0.0,
    z_center_mm=210.0,

    # Rexolite_WR229_109 158.615
    # PET_WR229_109       158.72
    # PVDF_WR229_109      160.17
    # KelF_WR229_109      158.93
    # Regolith_WR229_109  160.664
    # KelF_WR229_109      160.706
    # Acetal_WR229_109    161.578
    # Quartz_WR229_109    161.315
    #
    # PET_UHF_107         349.108
    # Regolith_UHF_107    340.422
    # Rexolite_UHF_107    344.471
    # KelF_UHF_107        193.747
    # KelF_UHF_107        195.356   raised planar sensor by 2 mm
    # PVDF_UHF_107        190.057
    # Quartz_UHF_107      192.536
    time_ns=161.315,
)


# ============================================================
# Neighboring-antinode normalization
# ============================================================
#
# The centered sample is at z=0.  For
#
#     E_empty = cos(pi x/a) cos(n pi z/ell),
#
# the nearest longitudinal antinodes are at
#
#     z = +/- ell/n.
#
# We use windows around BOTH of those antinodes to determine the underlying
# TE_10n standing-wave amplitude.  The same amplitude scale is then applied
# unchanged to BOTH the empty-cavity model and the loaded eigenmode model.
#
# This is preferable to fitting the sample region itself, because the local
# dielectric enhancement should be much weaker at the neighboring antinodes.
# ============================================================

# Nearest pair of antinodes: +/- ell/n.
neighbor_order = 1

# Width of the normalization windows.
# These are deliberately restricted around x=0 to avoid the side walls/feed.
antinode_x_halfwidth_mm = a / 4.0
antinode_z_halfwidth_mm = 0.15 * ell / n

print()
print("Chosen neighboring-antinode windows:")
print("  z centers = +/-", neighbor_order * ell / n, "mm")
print("  x half-width =", antinode_x_halfwidth_mm, "mm")
print("  z half-width =", antinode_z_halfwidth_mm, "mm")


# ============================================================
# Preferred comparison:
# one neighboring-antinode scale for BOTH theoretical models
# ============================================================

comparison, comparison_empty = compare_loaded_and_empty_neighbor_antinode(
    xf_data=xf,
    E_loaded_model=result["E"],
    a_mm=a,
    ell_mm=ell,
    target_n=n,
    wall_margin_mm=5.0,

    antinode_x_halfwidth_mm=antinode_x_halfwidth_mm,
    antinode_z_halfwidth_mm=antinode_z_halfwidth_mm,
    neighbor_order=neighbor_order,

    # Compare FIELD MAGNITUDES, not signed instantaneous values.
    # Positive residual => |XF| > |model|.
    residual_mode="magnitude",

    # Suppress nodes according to their electric-energy contribution.
    energy_weighted=True,
    percent_cutoff=0.01,

    # Use the real loaded material distribution for energy weighting.
    sample_radius_mm=Rs,
    eps_sample=eps_sample,
    vial_outer_radius_mm=Rh,
    eps_vial=eps_vial,
)


# Loaded eigenmode residual map
plot_comparison(
    comparison,
    xlim=a / 2,
    ylim=ell / n / 2,
    percent_limit=10,
)


# Empty-cavity residual map
plot_comparison(
    comparison_empty,
    xlim=a / 2,
    ylim=ell / n / 2,
    percent_limit=10,
)


# Inside-sample residual histograms for BOTH models
plot_sample_residual_histograms(
    comparison,
    comparison_empty,
    sample_radius_mm=Rs,
    bins=80,
    residual_limit=20,
)
