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

'''def solve_loaded_cavity(
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
    }'''




def solve_loaded_cavity(
    cavity,
    sample,
    n,
    m_max=21,
    p_max=61,
):
    """
    Solve the 2-D loaded rectangular-cavity eigenproblem

        K c = k^2 M c

    for the loaded mode connected to the empty TE_10n mode.

    NORMALIZATION
    -------------
    TE_101:
        There are no neighboring longitudinal antinodes inside
        the cavity, so the mode is normalized using the
        target-coefficient convention.

        The TE_101 component of the loaded eigenmode is scaled
        so that its contribution at the cavity center has
        amplitude +1.

    TE_10n for odd n >= 3:
        The loaded field is normalized using the two neighboring
        longitudinal antinodes

            x = 0
            z = +/- ell/n

        such that

            mean(
                |E(0,+ell/n)|,
                |E(0,-ell/n)|
            ) = 1.

    In all cases, the overall sign is chosen so that

        E(0,0) > 0.

    The function signature and returned CavitySolution are kept
    identical to the previous implementation.
    """

    # ========================================================
    # Basic checks
    # ========================================================

    if n % 2 == 0:
        raise ValueError(
            "This implementation assumes an odd TE_10n target mode."
        )

    if sample.radius_mm <= 0:
        raise ValueError(
            "sample.radius_mm must be positive."
        )

    if sample.eps_r <= 0:
        raise ValueError(
            "sample.eps_r must be positive."
        )

    if sample.vial_outer_radius_mm is not None:

        if sample.vial_outer_radius_mm <= sample.radius_mm:
            raise ValueError(
                "vial_outer_radius_mm must exceed "
                "sample.radius_mm."
            )

        if sample.eps_vial <= 0:
            raise ValueError(
                "sample.eps_vial must be positive."
            )

    # ========================================================
    # Build empty-cavity basis
    # ========================================================

    modes = make_basis(
        cavity,
        m_max=m_max,
        p_max=p_max,
    )

    N = len(modes)

    k2_empty = np.array(
        [
            mode[2]
            for mode in modes
        ],
        dtype=float,
    )

    K = np.diag(
        k2_empty
    )

    # ========================================================
    # Build dielectric/permittivity matrix
    #
    #       M_ij = integral epsilon_r phi_i phi_j dA
    # ========================================================

    I = np.eye(N)

    S_sample = build_disk_overlap_matrix(
        modes,
        cavity,
        sample.radius_mm,
    )

    # --------------------------------------------------------
    # No vial
    # --------------------------------------------------------

    if sample.vial_outer_radius_mm is None:

        M = (
            I
            +
            (sample.eps_r - 1.0)
            * S_sample
        )

    # --------------------------------------------------------
    # Sample inside concentric vial
    # --------------------------------------------------------

    else:

        S_vial = build_disk_overlap_matrix(
            modes,
            cavity,
            sample.vial_outer_radius_mm,
        )

        M = (
            I
            +
            (sample.eps_vial - 1.0)
            * S_vial
            +
            (sample.eps_r - sample.eps_vial)
            * S_sample
        )

    # ========================================================
    # Solve generalized eigenvalue problem
    #
    #       K c = k^2 M c
    # ========================================================

    eigvals, eigvecs = eigh(
        K,
        M,
    )

    # ========================================================
    # Locate the empty TE_10n basis function
    # ========================================================

    target_index = None

    for i, (m, p, _) in enumerate(modes):

        if (
            m == 1
            and
            p == n
        ):

            target_index = i
            break

    if target_index is None:

        raise ValueError(
            "Target TE_10n is outside the selected basis. "
            "Increase p_max."
        )

    # ========================================================
    # Identify the loaded eigenmode most strongly associated
    # with the empty TE_10n mode
    # ========================================================

    target_projection = np.abs(
        eigvecs[
            target_index,
            :
        ]
    )

    loaded_index = int(
        np.argmax(
            target_projection
        )
    )

    k2_loaded = float(
        eigvals[
            loaded_index
        ]
    )

    if k2_loaded <= 0:

        raise RuntimeError(
            "Loaded eigenvalue is non-positive."
        )

    k_loaded = np.sqrt(
        k2_loaded
    )

    coeff_raw = eigvecs[
        :,
        loaded_index
    ].copy()

    # ========================================================
    # NORMALIZATION
    # ========================================================

    # --------------------------------------------------------
    # TE_101
    #
    # There are no neighboring longitudinal antinodes inside
    # the cavity.
    #
    # Therefore use target-coefficient normalization.
    # --------------------------------------------------------

    if n == 1:

        target_coeff = coeff_raw[
            target_index
        ]

        if np.isclose(
            target_coeff,
            0.0
        ):

            raise RuntimeError(
                "Target TE_101 coefficient is zero; "
                "cannot normalize."
            )

        # ----------------------------------------------------
        # First make the coefficient of the TE_101 basis
        # function equal to one.
        # ----------------------------------------------------

        coeff = (
            coeff_raw
            /
            target_coeff
        )

        # ----------------------------------------------------
        # The normalized rectangular basis function has
        # center amplitude
        #
        #       2 / sqrt(a*ell)
        #
        # multiplied by sin(n*pi/2).
        #
        # Scale so that the TE_101 contribution at the center
        # has amplitude +1.
        # ----------------------------------------------------

        a = float(
            cavity.a_mm
        )

        ell = float(
            cavity.ell_mm
        )

        center_sign = np.sin(
            n * np.pi / 2.0
        )

        coeff *= (
            center_sign
            * np.sqrt(
                a * ell
            )
            / 2.0
        )

        # ----------------------------------------------------
        # These fields existed in the previous return object.
        #
        # For TE_101 there are no physical neighboring
        # antinodes, so store NaN instead of evaluating
        # outside the cavity.
        # ----------------------------------------------------

        E_plus_raw = np.nan
        E_minus_raw = np.nan
        scale = np.nan

    # --------------------------------------------------------
    # TE_103, TE_105, TE_107, ...
    #
    # Normalize using neighboring antinodes.
    # --------------------------------------------------------

    else:

        z_neighbor = (
            cavity.ell_mm
            /
            n
        )

        # ----------------------------------------------------
        # Raw field at + neighboring antinode
        # ----------------------------------------------------

        E_plus_raw = float(
            np.asarray(
                _evaluate_coefficients(
                    coeff_raw,
                    modes,
                    cavity,
                    0.0,
                    +z_neighbor,
                )
            )
        )

        # ----------------------------------------------------
        # Raw field at - neighboring antinode
        # ----------------------------------------------------

        E_minus_raw = float(
            np.asarray(
                _evaluate_coefficients(
                    coeff_raw,
                    modes,
                    cavity,
                    0.0,
                    -z_neighbor,
                )
            )
        )

        # ----------------------------------------------------
        # Average magnitude of the two neighboring antinodes
        # ----------------------------------------------------

        neighbor_mean_magnitude = (
            0.5
            *
            (
                abs(E_plus_raw)
                +
                abs(E_minus_raw)
            )
        )

        if neighbor_mean_magnitude <= 0:

            raise RuntimeError(
                "Neighboring-antinode amplitude is zero; "
                "cannot normalize the field."
            )

        # ----------------------------------------------------
        # Scale whole eigenvector so neighboring antinodes
        # have mean |E| = 1
        # ----------------------------------------------------

        scale = (
            1.0
            /
            neighbor_mean_magnitude
        )

        coeff = (
            coeff_raw
            *
            scale
        )

    # ========================================================
    # Fix arbitrary eigenvector sign
    #
    # We always choose
    #
    #       E(0,0) > 0
    #
    # so the loaded and empty fields use the same center sign.
    # ========================================================

    E_center = float(
        np.asarray(
            _evaluate_coefficients(
                coeff,
                modes,
                cavity,
                0.0,
                0.0,
            )
        )
    )

    if E_center < 0:

        coeff *= -1.0

    # ========================================================
    # Empty-cavity quantities
    # ========================================================

    k_empty = empty_te10n_wavenumber(
        cavity,
        n,
    )

    f_empty = empty_te10n_frequency(
        cavity,
        n,
    )

    # ========================================================
    # Loaded resonance frequency
    # ========================================================

    c0 = 299_792_458.0

    f_loaded = (
        c0
        *
        (
            1000.0
            * k_loaded
        )
        /
        (
            2.0
            * np.pi
        )
    )

    # ========================================================
    # Return same CavitySolution object as before
    # ========================================================

    return CavitySolution(

        cavity=cavity,

        sample=sample,

        n=n,

        modes=modes,

        coefficients=coeff,

        k_loaded_mm_inv=
            k_loaded,

        k_empty_mm_inv=
            k_empty,

        f_loaded_Hz=
            f_loaded,

        f_empty_Hz=
            f_empty,

        loaded_mode_index=
            loaded_index,

        target_basis_index=
            target_index,

        neighbor_antinode_raw_plus=
            E_plus_raw,

        neighbor_antinode_raw_minus=
            E_minus_raw,

        neighbor_antinode_scale=
            scale,
    )





import numpy as np


def find_linear_region(
    x,
    y,
    min_points=5,
    r2_threshold=0.995,
    start_sample=0
):
    """
    Find the longest contiguous region of a dataset that is well described
    by a straight line, while optionally ignoring the beginning of the data.

    The function searches through every possible contiguous region starting
    at or after `start_sample`.

    For each candidate region, a linear fit

        y = slope*x + intercept

    is performed.

    A region is considered sufficiently linear if

        R^2 >= r2_threshold

    Among all regions satisfying this requirement, the function chooses the
    region containing the largest number of points.

    If two acceptable regions contain the same number of points, the one
    with the larger R^2 is selected.


    Parameters
    ----------
    x : array-like
        Independent-variable data.

    y : array-like
        Dependent-variable data.

    min_points : int, optional
        Minimum number of points that must be included in a candidate
        linear region.

        Default:
            min_points = 5


    r2_threshold : float, optional
        Minimum R^2 required for a region to be considered linear.

        Example values:
            0.99
            0.995
            0.999

        Default:
            r2_threshold = 0.995


    start_sample : int, optional
        Index of the first data point that is allowed to be considered
        when searching for a linear region.

        All points before this index are completely ignored during the
        search.

        For example:

            start_sample = 0

        means search the entire dataset.

            start_sample = 20

        means indices 0 through 19 will be ignored, and the search will
        begin at index 20.

        This is useful if the beginning of the data corresponds to a
        region where the sample has not yet entered the cavity.

        Default:
            start_sample = 0


    Returns
    -------
    best_region : dict

        Dictionary containing:

        "x"
            x-values in the selected linear region.

        "y"
            y-values in the selected linear region.

        "indices"
            Tuple containing the starting and ending indices in the
            ORIGINAL x and y arrays.

        "start_index"
            Starting index of the selected region.

        "end_index"
            Ending index of the selected region.

        "slope"
            Slope of the linear fit.

        "intercept"
            Intercept of the linear fit.

        "r2"
            R^2 value of the selected fit.

        "n_points"
            Number of points in the selected region.

        "y_fit"
            Fitted y-values evaluated at the selected x-values.
    """


    # ===============================================================
    # CONVERT INPUTS TO NUMPY ARRAYS
    # ===============================================================

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)


    # ===============================================================
    # BASIC INPUT CHECKS
    # ===============================================================

    # x and y must contain the same number of measurements.
    if len(x) != len(y):
        raise ValueError(
            "x and y must contain the same number of data points."
        )


    # start_sample must be an integer index.
    if not isinstance(start_sample, (int, np.integer)):
        raise TypeError(
            "start_sample must be an integer."
        )


    # Negative indices are not allowed here because we want
    # start_sample to represent a clear physical starting point
    # in the dataset.
    if start_sample < 0:
        raise ValueError(
            "start_sample must be greater than or equal to 0."
        )


    # Make sure start_sample actually lies inside the dataset.
    if start_sample >= len(x):
        raise ValueError(
            f"start_sample={start_sample} is outside the dataset. "
            f"The dataset contains {len(x)} points."
        )


    # Make sure enough points remain after start_sample to construct
    # at least one candidate region.
    if len(x) - start_sample < min_points:
        raise ValueError(
            f"Only {len(x) - start_sample} points remain after "
            f"start_sample={start_sample}, but min_points={min_points}."
        )


    # ===============================================================
    # TOTAL NUMBER OF DATA POINTS
    # ===============================================================

    N = len(x)


    # ===============================================================
    # STORAGE FOR THE BEST REGION FOUND
    # ===============================================================

    # We have not yet found an acceptable region, so begin with None.
    best_region = None


    # ===============================================================
    # SEARCH THROUGH ALL POSSIBLE CONTIGUOUS REGIONS
    # ===============================================================
    #
    # The important change here is:
    #
    #     range(start_sample, N)
    #
    # instead of:
    #
    #     range(0, N)
    #
    # Therefore NO candidate region can begin before start_sample.
    #
    # Example:
    #
    #     start_sample = 20
    #
    # Candidate regions might be:
    #
    #     20 -> 24
    #     20 -> 25
    #     20 -> 26
    #     ...
    #
    #     21 -> 25
    #     21 -> 26
    #     ...
    #
    # But something like:
    #
    #     5 -> 30
    #
    # is NEVER considered.
    #
    # Thus the initial portion of the data has absolutely no effect
    # on the fitted linear region.
    # ===============================================================

    for start_index in range(start_sample, N):


        # -----------------------------------------------------------
        # The smallest candidate region must contain min_points.
        #
        # Because Python slicing excludes the final value, the first
        # allowed end point is:
        #
        #     start_index + min_points
        # -----------------------------------------------------------

        for end_index_exclusive in range(
            start_index + min_points,
            N + 1
        ):


            # =======================================================
            # EXTRACT THIS CANDIDATE REGION
            # =======================================================

            x_region = x[
                start_index:end_index_exclusive
            ]

            y_region = y[
                start_index:end_index_exclusive
            ]


            # =======================================================
            # PERFORM LINEAR FIT
            # =======================================================
            #
            # np.polyfit(..., 1) fits
            #
            #     y = slope*x + intercept
            # =======================================================

            slope, intercept = np.polyfit(
                x_region,
                y_region,
                1
            )


            # =======================================================
            # EVALUATE THE FITTED LINE
            # =======================================================

            y_fit = (
                slope * x_region
                + intercept
            )


            # =======================================================
            # CALCULATE RESIDUAL SUM OF SQUARES
            # =======================================================
            #
            # This tells us how far the measurements lie from the
            # fitted straight line.
            # =======================================================

            ss_residual = np.sum(
                (y_region - y_fit)**2
            )


            # =======================================================
            # CALCULATE TOTAL SUM OF SQUARES
            # =======================================================
            #
            # This tells us how much the y-data vary around their mean.
            # =======================================================

            ss_total = np.sum(
                (
                    y_region
                    - np.mean(y_region)
                )**2
            )


            # =======================================================
            # CALCULATE R^2
            # =======================================================

            if ss_total == 0:

                # If every y-value is identical and the fit is exact,
                # treat it as perfectly linear.
                if ss_residual == 0:
                    r2 = 1.0

                else:
                    r2 = 0.0

            else:

                r2 = (
                    1
                    - ss_residual / ss_total
                )


            # =======================================================
            # NUMBER OF POINTS IN THIS CANDIDATE REGION
            # =======================================================

            n_points = len(x_region)


            # =======================================================
            # CHECK IF THE REGION IS LINEAR ENOUGH
            # =======================================================

            if r2 >= r2_threshold:


                # ---------------------------------------------------
                # If this is the first acceptable region, keep it.
                # ---------------------------------------------------

                if best_region is None:

                    choose_this_region = True


                # ---------------------------------------------------
                # Prefer a candidate containing MORE points.
                #
                # This makes the algorithm favor the longest region
                # that still satisfies the required R^2.
                # ---------------------------------------------------

                elif n_points > best_region["n_points"]:

                    choose_this_region = True


                # ---------------------------------------------------
                # If both regions contain the same number of points,
                # prefer the one with the better R^2.
                # ---------------------------------------------------

                elif (
                    n_points == best_region["n_points"]
                    and
                    r2 > best_region["r2"]
                ):

                    choose_this_region = True


                else:

                    choose_this_region = False


                # ===================================================
                # SAVE THIS REGION IF IT IS CURRENTLY THE BEST
                # ===================================================

                if choose_this_region:


                    # Python slicing excludes the final endpoint,
                    # so subtract one to get the actual ending index.
                    final_index = (
                        end_index_exclusive - 1
                    )


                    best_region = {

                        "x": x_region.copy(),

                        "y": y_region.copy(),

                        "indices": (
                            start_index,
                            final_index
                        ),

                        "start_index": start_index,

                        "end_index": final_index,

                        "slope": slope,

                        "intercept": intercept,

                        "r2": r2,

                        "n_points": n_points,

                        "y_fit": y_fit.copy(),

                        # Save this too, just so the output records
                        # what restriction was used during the search.
                        "search_start_sample": start_sample,
                    }


    # ===============================================================
    # CHECK WHETHER ANY ACCEPTABLE REGION WAS FOUND
    # ===============================================================

    if best_region is None:

        raise ValueError(
            "No sufficiently linear region was found.\n"
            f"Search began at sample {start_sample}.\n"
            f"No region containing at least {min_points} points "
            f"had R^2 >= {r2_threshold}."
        )


    # ===============================================================
    # RETURN THE SELECTED REGION
    # ===============================================================

    return best_region