"""
cavity_sample_field.py

Wall-aware 2-D eigenmode model for a centered cylindrical dielectric sample
inside a rectangular cavity operating in an odd TE_10n mode.

NORMALIZATION
-------------
The loaded field is normalized using the TWO NEIGHBORING longitudinal
antinodes of the target TE_10n mode:

    x = 0
    z = +/- ell / n

The scale factor is chosen so that the MEAN magnitude at those two
neighboring antinodes is exactly 1:

    mean( |E(0,+ell/n)| , |E(0,-ell/n)| ) = 1

The overall sign is then chosen so that E(0,0) > 0.

All lengths are in mm.
"""

from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import eigh
from scipy.special import j1


@dataclass(frozen=True)
class Cavity:
    a_mm: float
    ell_mm: float


@dataclass(frozen=True)
class Sample:
    radius_mm: float
    eps_r: float
    vial_outer_radius_mm: float | None = None
    eps_vial: float = 1.0


@dataclass
class CavitySolution:
    cavity: Cavity
    sample: Sample
    n: int
    modes: list
    coefficients: np.ndarray
    k_loaded_mm_inv: float
    k_empty_mm_inv: float
    f_loaded_Hz: float
    f_empty_Hz: float
    loaded_mode_index: int
    target_basis_index: int
    neighbor_antinode_raw_plus: float
    neighbor_antinode_raw_minus: float
    neighbor_antinode_scale: float


def empty_te10n_field(x_mm, z_mm, cavity, n):
    if n % 2 == 0:
        raise ValueError("This centered formulation is intended for odd n.")

    a = float(cavity.a_mm)
    ell = float(cavity.ell_mm)

    x = np.asarray(x_mm, dtype=float)
    z = np.asarray(z_mm, dtype=float)

    return (
        np.cos(np.pi * x / a)
        * np.cos(n * np.pi * z / ell)
    )


def empty_te10n_wavenumber(cavity, n):
    a = float(cavity.a_mm)
    ell = float(cavity.ell_mm)

    return np.sqrt(
        (np.pi / a)**2
        + (n * np.pi / ell)**2
    )


def empty_te10n_frequency(cavity, n):
    c0 = 299_792_458.0
    k_mm = empty_te10n_wavenumber(cavity, n)
    return c0 * (1000.0 * k_mm) / (2.0 * np.pi)


def cavity_basis_function(m, p, x_mm, z_mm, cavity):
    a = float(cavity.a_mm)
    ell = float(cavity.ell_mm)

    x = np.asarray(x_mm, dtype=float)
    z = np.asarray(z_mm, dtype=float)

    x_global = x + a / 2.0
    z_global = z + ell / 2.0

    return (
        2.0 / np.sqrt(a * ell)
        * np.sin(m * np.pi * x_global / a)
        * np.sin(p * np.pi * z_global / ell)
    )


def make_basis(cavity, m_max=21, p_max=61):
    a = float(cavity.a_mm)
    ell = float(cavity.ell_mm)

    modes = []

    for m in range(1, int(m_max) + 1, 2):
        for p in range(1, int(p_max) + 1, 2):
            k2 = (
                (m * np.pi / a)**2
                + (p * np.pi / ell)**2
            )
            modes.append((m, p, k2))

    return modes


def disk_fourier_integral(Q_mm_inv, radius_mm):
    Q = float(Q_mm_inv)
    R = float(radius_mm)

    if np.isclose(Q, 0.0):
        return np.pi * R**2

    return 2.0 * np.pi * R * j1(Q * R) / Q


def disk_mode_overlap(m1, p1, m2, p2, cavity, radius_mm):
    a = float(cavity.a_mm)
    ell = float(cavity.ell_mm)
    R = float(radius_mm)

    s_m1 = np.sin(m1 * np.pi / 2.0)
    s_m2 = np.sin(m2 * np.pi / 2.0)
    s_p1 = np.sin(p1 * np.pi / 2.0)
    s_p2 = np.sin(p2 * np.pi / 2.0)

    sign = s_m1 * s_m2 * s_p1 * s_p2

    total = 0.0

    for sx in (-1, 1):
        qx = (m1 + sx * m2) * np.pi / a

        for sz in (-1, 1):
            qz = (p1 + sz * p2) * np.pi / ell
            Q = np.sqrt(qx**2 + qz**2)
            total += disk_fourier_integral(Q, R)

    return sign * total / (a * ell)


def build_disk_overlap_matrix(modes, cavity, radius_mm):
    N = len(modes)
    S = np.zeros((N, N), dtype=float)

    for i in range(N):
        m1, p1, _ = modes[i]

        for j in range(i, N):
            m2, p2, _ = modes[j]

            value = disk_mode_overlap(
                m1, p1, m2, p2,
                cavity, radius_mm
            )

            S[i, j] = value
            S[j, i] = value

    return S


def _evaluate_coefficients(coefficients, modes, cavity, x_mm, z_mm):
    x = np.asarray(x_mm, dtype=float)
    z = np.asarray(z_mm, dtype=float)

    shape = np.broadcast(x, z).shape
    E = np.zeros(shape, dtype=float)

    for c, (m, p, _) in zip(coefficients, modes):
        E += c * cavity_basis_function(
            m, p, x, z, cavity
        )

    return E


def solve_loaded_cavity(
    cavity,
    sample,
    n,
    m_max=21,
    p_max=61,
):
    if n % 2 == 0:
        raise ValueError(
            "This implementation assumes an odd TE_10n target mode."
        )

    if sample.radius_mm <= 0:
        raise ValueError("sample.radius_mm must be positive.")

    if sample.eps_r <= 0:
        raise ValueError("sample.eps_r must be positive.")

    if sample.vial_outer_radius_mm is not None:
        if sample.vial_outer_radius_mm <= sample.radius_mm:
            raise ValueError(
                "vial_outer_radius_mm must exceed sample.radius_mm."
            )

    modes = make_basis(
        cavity,
        m_max=m_max,
        p_max=p_max,
    )

    N = len(modes)

    k2_empty = np.array(
        [mode[2] for mode in modes],
        dtype=float,
    )

    K = np.diag(k2_empty)
    I = np.eye(N)

    S_sample = build_disk_overlap_matrix(
        modes,
        cavity,
        sample.radius_mm,
    )

    if sample.vial_outer_radius_mm is None:
        M = (
            I
            + (sample.eps_r - 1.0)
            * S_sample
        )
    else:
        S_vial = build_disk_overlap_matrix(
            modes,
            cavity,
            sample.vial_outer_radius_mm,
        )

        M = (
            I
            + (sample.eps_vial - 1.0) * S_vial
            + (sample.eps_r - sample.eps_vial) * S_sample
        )

    eigvals, eigvecs = eigh(K, M)

    target_index = None

    for i, (m, p, _) in enumerate(modes):
        if m == 1 and p == n:
            target_index = i
            break

    if target_index is None:
        raise ValueError(
            "Target TE_10n is outside the selected basis. Increase p_max."
        )

    loaded_index = int(
        np.argmax(
            np.abs(
                eigvecs[target_index, :]
            )
        )
    )

    k2_loaded = float(eigvals[loaded_index])

    if k2_loaded <= 0:
        raise RuntimeError("Loaded eigenvalue is non-positive.")

    k_loaded = np.sqrt(k2_loaded)

    coeff_raw = eigvecs[:, loaded_index].copy()

    # Neighboring antinodes of central antinode: z = +/- ell/n.
    z_neighbor = cavity.ell_mm / n

    E_plus_raw = float(
        np.asarray(
            _evaluate_coefficients(
                coeff_raw, modes, cavity,
                0.0, +z_neighbor
            )
        )
    )

    E_minus_raw = float(
        np.asarray(
            _evaluate_coefficients(
                coeff_raw, modes, cavity,
                0.0, -z_neighbor
            )
        )
    )

    neighbor_mean_magnitude = 0.5 * (
        abs(E_plus_raw)
        + abs(E_minus_raw)
    )

    if neighbor_mean_magnitude <= 0:
        raise RuntimeError(
            "Neighboring-antinode amplitude is zero; cannot normalize."
        )

    scale = 1.0 / neighbor_mean_magnitude
    coeff = coeff_raw * scale

    # Fix arbitrary sign so the central antinode is positive.
    E_center = float(
        np.asarray(
            _evaluate_coefficients(
                coeff, modes, cavity,
                0.0, 0.0
            )
        )
    )

    if E_center < 0:
        coeff *= -1.0

    c0 = 299_792_458.0

    f_loaded = (
        c0
        * (1000.0 * k_loaded)
        / (2.0 * np.pi)
    )

    k_empty = empty_te10n_wavenumber(
        cavity,
        n,
    )

    f_empty = empty_te10n_frequency(
        cavity,
        n,
    )

    return CavitySolution(
        cavity=cavity,
        sample=sample,
        n=n,
        modes=modes,
        coefficients=coeff,
        k_loaded_mm_inv=k_loaded,
        k_empty_mm_inv=k_empty,
        f_loaded_Hz=f_loaded,
        f_empty_Hz=f_empty,
        loaded_mode_index=loaded_index,
        target_basis_index=target_index,
        neighbor_antinode_raw_plus=E_plus_raw,
        neighbor_antinode_raw_minus=E_minus_raw,
        neighbor_antinode_scale=scale,
    )


def loaded_field(solution, x_mm, z_mm):
    return _evaluate_coefficients(
        solution.coefficients,
        solution.modes,
        solution.cavity,
        x_mm,
        z_mm,
    )


def loaded_field_magnitude(solution, x_mm, z_mm):
    return np.abs(
        loaded_field(
            solution,
            x_mm,
            z_mm,
        )
    )


def sample_contains(solution, x_mm, z_mm):
    x = np.asarray(x_mm, dtype=float)
    z = np.asarray(z_mm, dtype=float)

    return (
        x**2 + z**2
        <= solution.sample.radius_mm**2
    )


def field_inside_sample(
    solution,
    x_mm,
    z_mm,
    outside_value=np.nan,
):
    E = np.asarray(
        loaded_field(
            solution,
            x_mm,
            z_mm,
        )
    )

    inside = np.asarray(
        sample_contains(
            solution,
            x_mm,
            z_mm,
        )
    )

    result = np.where(
        inside,
        E,
        outside_value,
    )

    if result.ndim == 0:
        return float(result)

    return result


def center_field(solution):
    return float(
        np.asarray(
            loaded_field(
                solution,
                0.0,
                0.0,
            )
        )
    )


def neighboring_antinode_fields(solution):
    z0 = solution.cavity.ell_mm / solution.n

    E_plus = float(
        np.asarray(
            loaded_field(
                solution,
                0.0,
                +z0,
            )
        )
    )

    E_minus = float(
        np.asarray(
            loaded_field(
                solution,
                0.0,
                -z0,
            )
        )
    )

    return {
        "z_mm": z0,
        "E_plus": E_plus,
        "E_minus": E_minus,
        "abs_plus": abs(E_plus),
        "abs_minus": abs(E_minus),
        "mean_abs": 0.5 * (
            abs(E_plus)
            + abs(E_minus)
        ),
    }


def field_enhancement_over_empty(
    solution,
    x_mm,
    z_mm,
):
    E_loaded = np.asarray(
        loaded_field_magnitude(
            solution,
            x_mm,
            z_mm,
        ),
        dtype=float,
    )

    E_empty = np.asarray(
        np.abs(
            empty_te10n_field(
                x_mm,
                z_mm,
                solution.cavity,
                solution.n,
            )
        ),
        dtype=float,
    )

    E_loaded, E_empty = np.broadcast_arrays(
        E_loaded,
        E_empty,
    )

    ratio = np.full(
        E_loaded.shape,
        np.nan,
        dtype=float,
    )

    valid = E_empty > 1e-12

    ratio[valid] = (
        E_loaded[valid]
        / E_empty[valid]
    )

    if ratio.ndim == 0:
        return float(ratio)

    return ratio


def make_sample_grid(solution, points=201):
    R = float(solution.sample.radius_mm)

    coords = np.linspace(
        -R,
        R,
        int(points),
    )

    X, Z = np.meshgrid(
        coords,
        coords,
    )

    return X, Z


def sample_field_map(solution, points=201):
    X, Z = make_sample_grid(
        solution,
        points=points,
    )

    E = field_inside_sample(
        solution,
        X,
        Z,
    )

    return X, Z, E


def plot_sample_field(
    solution,
    points=201,
    magnitude=False,
):
    X, Z, E = sample_field_map(
        solution,
        points=points,
    )

    if magnitude:
        data = np.abs(E)
        label = "|E|"
        title = "Loaded field magnitude inside sample"
        cmap = "viridis"
        vmin = None
        vmax = None
    else:
        data = E
        label = "Normalized E"
        title = "Loaded electric field inside sample"
        cmap = "RdBu_r"

        peak = np.nanmax(
            np.abs(data)
        )

        vmin = -peak
        vmax = +peak

    fig, ax = plt.subplots(
        figsize=(7, 6)
    )

    pcm = ax.pcolormesh(
        X,
        Z,
        data,
        shading="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )

    fig.colorbar(
        pcm,
        ax=ax,
        label=label,
    )

    circle = plt.Circle(
        (0.0, 0.0),
        solution.sample.radius_mm,
        fill=False,
        linewidth=1.5,
    )

    ax.add_patch(circle)

    ax.set_xlabel(
        "x from sample center [mm]"
    )

    ax.set_ylabel(
        "z from sample center [mm]"
    )

    ax.set_title(title)
    ax.set_aspect("equal")

    fig.tight_layout()

    return fig, ax


def plot_center_line(
    solution,
    axis="z",
    points=501,
    include_empty=True,
):
    R = float(solution.sample.radius_mm)

    s = np.linspace(
        -R,
        R,
        int(points),
    )

    if axis.lower() == "z":
        x = np.zeros_like(s)
        z = s
        xlabel = "z from sample center [mm]"

    elif axis.lower() == "x":
        x = s
        z = np.zeros_like(s)
        xlabel = "x from sample center [mm]"

    else:
        raise ValueError(
            "axis must be 'x' or 'z'."
        )

    E_loaded = loaded_field(
        solution,
        x,
        z,
    )

    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    ax.plot(
        s,
        E_loaded,
        label="Loaded eigenmode",
    )

    if include_empty:
        E_empty = empty_te10n_field(
            x,
            z,
            solution.cavity,
            solution.n,
        )

        ax.plot(
            s,
            E_empty,
            linestyle="--",
            label="Empty cavity",
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Normalized E")
    ax.set_title("Field through sample center")
    ax.legend()
    ax.grid(alpha=0.25)

    fig.tight_layout()

    return fig, ax


if __name__ == "__main__":

    cavity = Cavity(
        a_mm=6.0 * 25.4,
        ell_mm=18.0 * 25.4,
    )

    sample = Sample(
        radius_mm=4.5,
        eps_r=3.08,
        vial_outer_radius_mm=4.55,
        eps_vial=1.0,
    )

    n = 7

    solution = solve_loaded_cavity(
        cavity=cavity,
        sample=sample,
        n=n,
        m_max=31,
        p_max=81,
    )

    print()
    print("Empty frequency =", solution.f_empty_Hz / 1e9, "GHz")
    print("Loaded frequency =", solution.f_loaded_Hz / 1e9, "GHz")
    print(
        "Frequency shift =",
        (solution.f_empty_Hz - solution.f_loaded_Hz) / 1e6,
        "MHz",
    )

    neighbors = neighboring_antinode_fields(solution)

    print()
    print("Neighboring-antinode normalization")
    print("z_neighbor =", neighbors["z_mm"], "mm")
    print("E(+z_neighbor) =", neighbors["E_plus"])
    print("E(-z_neighbor) =", neighbors["E_minus"])
    print("Mean neighboring |E| =", neighbors["mean_abs"])

    print()
    print("Loaded center field =", center_field(solution))
    print(
        "Empty center field =",
        empty_te10n_field(0.0, 0.0, cavity, n),
    )
    print(
        "Center |E_loaded|/|E_empty| =",
        field_enhancement_over_empty(
            solution,
            0.0,
            0.0,
        ),
    )

    x_test = 1.0
    z_test = 2.0

    print()
    print(
        f"E({x_test} mm, {z_test} mm) =",
        field_inside_sample(
            solution,
            x_test,
            z_test,
        ),
    )

    plot_sample_field(
        solution,
        magnitude=False,
    )

    plot_center_line(
        solution,
        axis="z",
        include_empty=True,
    )

    plt.show()
