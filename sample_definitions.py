# sample_definitions.py
#
# Centralized definitions for cavity resonator samples.
#
# Each dictionary entry contains the sample-specific parameters used by
# the main analysis script. The Sample object itself should be created
# in the analysis script using these stored values.


SAMPLES = {

    # ==============================================================
    # Teflon1
    # ==============================================================

    "Teflon1": {
        "sample_range": [1.9, 2.1],
        "radius_mm": 6.79 / 2,
        "eps_r": 2.04,
        "vial_outer_radius_mm": 9.08 / 2,
        "eps_vial": 2.54,
    },


    # ==============================================================
    # KelF1
    # ==============================================================

    "KelF1": {
        "sample_range": [2.2, 2.4],
        "radius_mm": 6.85 / 2,
        "eps_r": 2.35,
        "vial_outer_radius_mm": 9.08 / 2,
        "eps_vial": 2.54,
    },


    # ==============================================================
    # Rexolite
    # ==============================================================

    "Rexolite": {
        "sample_range": [2.4, 2.6],
        "radius_mm": 6.87 / 2,
        "eps_r": 2.54,
        "vial_outer_radius_mm": 9.08 / 2,
        "eps_vial": 2.54,

        # Previously noted possible imaginary-permittivity range:
        # "sample_i_range": [1e-4, 2e-2],
    },


    # ==============================================================
    # Acrylic
    # ==============================================================

    "Acrylic": {
        "sample_range": [2.5, 2.8],
        "radius_mm": 6.79 / 2,
        "eps_r": 2.56,
        "vial_outer_radius_mm": 9.08 / 2,
        "eps_vial": 2.54,

        # Previously noted possible imaginary-permittivity range:
        # "sample_i_range": [1e-3, 2e-1],
    },


    # ==============================================================
    # Acetal1
    # ==============================================================

    "Acetal1": {
        "sample_range": [2.8, 3.1],
        "radius_mm": 6.85 / 2,
        "eps_r": 2.90,
        "vial_outer_radius_mm": 9.08 / 2,
        "eps_vial": 2.54,

        # Previously noted possible imaginary-permittivity range:
        # "sample_i_range": [1e-3, 2e-1],
    },


    # ==============================================================
    # PET
    # ==============================================================

    "PET": {
        "sample_range": [2.9, 3.2],
        "radius_mm": 6.8 / 2,
        "eps_r": 3.08,
        "vial_outer_radius_mm": 9.08 / 2,
        "eps_vial": 2.54,

        # Original note:
        # LunarSample1 Rs = 7.06/2, Rv = 9.09/2
    },


    # ==============================================================
    # PVDF1
    # ==============================================================

    "PVDF1": {
        "sample_range": [2.7, 2.9],
        "radius_mm": 6.87 / 2,
        "eps_r": 2.83,
        "vial_outer_radius_mm": 9.08 / 2,
        "eps_vial": 2.54,
    },

}
