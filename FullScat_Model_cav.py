import numpy as np
import matplotlib.pyplot as plt

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
# Example
# ============================================================

if __name__ == "__main__":

    # --------------------------------------------------------
    # Cavity dimensions [mm]
    # --------------------------------------------------------

    a = 6.0 * 25.4
    ell = 18.0 * 25.4

    a = 58.16
    ell = 304.8

    # TE_107
    n = 9

    # --------------------------------------------------------
    # Sample
    # --------------------------------------------------------

    Rs = 9/2
    eps_sample = 3.0

    # --------------------------------------------------------
    # Vial
    # --------------------------------------------------------

    Rh = 9.11/2
    eps_vial = 1.0

    # --------------------------------------------------------
    # Solve loaded cavity
    # --------------------------------------------------------

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

    #print()
    #print(
    #    "Empty frequency:",
    #    result["f_empty_Hz"] / 1e9,
    #    "GHz"
    #)

    #print(
    #    "Loaded frequency:",
    #    result["f_loaded_Hz"] / 1e9,
    #    "GHz"
    #)

    #print(
    #    "Frequency shift:",
    #    result["frequency_shift_Hz"] / 1e6,
    #    "MHz"
    #)

    #print(
    #    "Fractional shift:",
    #    result["fractional_shift"]
    #)

    # --------------------------------------------------------
    # Field at center
    # --------------------------------------------------------

    E_center = result["E"](
        0.0,
        0.0,
    )

    #print()
    #print(
    #    "Loaded field at sample center:",
    #    E_center
    #)

    #print(
    #    "Empty field at sample center:",
    #    result["E_empty"](0.0, 0.0)
    #)

    # --------------------------------------------------------
    # Make field map
    # --------------------------------------------------------

    Nx = 350
    Nz = 700

    X = np.linspace(
        -a/2,
        a/2,
        Nx,
    )

    Z = np.linspace(
        -ell/2,
        ell/2,
        Nz,
    )

    XX, ZZ = np.meshgrid(
        X,
        Z,
    )

    E = result["E"](
        XX,
        ZZ,
    )

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(6, 12)
    )

    pcm = ax.pcolormesh(
        X,
        Z,
        E,
        shading="auto",
        cmap="RdBu_r",
    )

    plt.colorbar(
        pcm,
        ax=ax,
        label="Normalized $E_y$",
    )

    # Sample circle
    sample_circle = plt.Circle(
        (0, 0),
        Rs,
        fill=False,
        color="black",
        linewidth=1.5,
    )

    ax.add_patch(
        sample_circle
    )

    # Vial outer circle
    vial_circle = plt.Circle(
        (0, 0),
        Rh,
        fill=False,
        color="black",
        linestyle="--",
        linewidth=1.5,
    )

    ax.add_patch(
        vial_circle
    )

    ax.set_xlabel(
        "x from cavity center [mm]"
    )

    ax.set_ylabel(
        "z from cavity center [mm]"
    )

    ax.set_title(
        f"Loaded TE_10{n} cavity field"
    )

    ax.set_aspect(
        "equal"
    )

    plt.tight_layout()
    plt.show()



    N = int(1e6)




    for n in [3, 5, 7, 9]:
        
        a10n_vs = []
        a10n_s = []

        for epsr in [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:



            xsamples = 2*Rs*np.random.random(int(N)) - Rs
            zsamples = 2*Rs*np.random.random(int(N)) - Rs

            mask = (xsamples**2 + zsamples**2 <= Rs**2)

            #print(len(xsamples), len(zsamples))
            I = solve_loaded_cavity(
            a_mm=a,
            ell_mm=ell,

            target_n=n,

            sample_radius_mm=Rs,
            eps_sample=epsr,

            vial_outer_radius_mm=Rh,
            eps_vial=eps_vial,

            m_max=21,
            p_max=61,

            normalization="target_coefficient",
            )["E"](xsamples[mask], zsamples[mask])*np.cos(np.pi*xsamples[mask]/a)*np.cos(n*np.pi*zsamples[mask]/ell)

            Iold = np.cos(np.pi*xsamples[mask]/a)*np.cos(n*np.pi*zsamples[mask]/ell)*np.cos(np.pi*xsamples[mask]/a)*np.cos(n*np.pi*zsamples[mask]/ell)

            Aeff_old = np.mean(Iold)*np.pi*Rs**2

            print("Sample Effective Area is: ", np.mean(I)*np.pi*Rs**2)
            a10n_vs.append(np.mean(I))
            a10n_s.append(np.mean(I)/np.mean(Iold))

        p1, V1 = np.polyfit([1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0], a10n_s, 1, cov = True)

        m1 = p1[0]
        b1 = p1[1]
        print("Sample Correction function is: ", str(p1[0]), "* epsr0 + ", str(p1[1]))

        plt.errorbar([1, 1.5, 2, 2.5, 3, 3.5, 4], a10n_s, yerr = 0, marker = 'o', label = "n = " + str(n))
        plt.plot([1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0], p1[0]*np.array([1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]) + p1[1])
    plt.show()


    '''xsamples = a*np.random.random(int(1e7)) - a/2
    zsamples = ell*np.random.random(int(1e7)) - ell/2
    n = 7
    #mask = (xsamples**2 + zsamples**2 <= Rs**2)
    #print(len(xsamples), len(zsamples))
    I = solve_loaded_cavity(
    a_mm=a,
    ell_mm=ell,

    target_n=n,

    sample_radius_mm=4.5,
    eps_sample=2.54,

    vial_outer_radius_mm=Rh,
    eps_vial=eps_vial,

    m_max=21,
    p_max=61,

    normalization="target_coefficient",
    )["E"](xsamples, zsamples)*np.cos(np.pi*xsamples/a)*np.cos(n*np.pi*zsamples/ell)
    Iold = np.cos(np.pi*xsamples/a)*np.cos(n*np.pi*zsamples/ell)*np.cos(np.pi*xsamples/a)*np.cos(n*np.pi*zsamples/ell)
    print("Total volume difference of ", (np.mean(I) - np.mean(Iold))/np.mean(Iold)*100, "Percent")'''

