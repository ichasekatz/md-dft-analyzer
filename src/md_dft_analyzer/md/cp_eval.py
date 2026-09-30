"""Constant-pressure heat capacity (Cp) from LAMMPS NPT output files.

Reads a series of LAMMPS ``thermo_style custom`` output files sampled at
different temperatures, computes average enthalpy at each temperature, and
extracts Cp as the slope of the enthalpy–temperature curve via polynomial fit.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from md_dft_analyzer.numerics import poly_eval, poly_fit, r_squared

# LAMMPS thermo column names for NPT enthalpy runs
_lammps_columns = ["step", "ke", "pe", "te", "enthalpy", "vol", "pres", "pxx", "pyy", "pzz", "pxy", "vol2", "lx"]


def load_enthalpy_series(
    file_pattern: str,
    temperatures: list[int],
    skiprows: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load average enthalpy at each temperature from LAMMPS sample files.

    Expects one LAMMPS log file per temperature matching ``file_pattern``
    with ``{T}`` as a placeholder for the temperature value.

    Args:
        file_pattern (str): File path template with ``{T}`` placeholder,
            e.g. ``"outputs/samp_{T}.out"``.
        temperatures (list[int]): List of temperatures in K.
        skiprows (int): Header rows to skip in each file. Defaults to 1.

    Returns:
        Tuple ``(temperatures, avg_enthalpy, std_dev)`` as float arrays,
        where enthalpy is in eV/atom.

    Raises:
        FileNotFoundError: If any expected output file is missing.
    """
    avg_h: list[float] = []
    std_h: list[float] = []

    for temp in temperatures:
        path = Path(file_pattern.format(T=temp))
        if not path.exists():
            raise FileNotFoundError(path)
        df = pd.read_csv(path, skiprows=skiprows, sep=r"\s+", names=_lammps_columns)
        h = df["enthalpy"].to_numpy()
        avg_h.append(float(np.mean(h)))
        std_h.append(float(np.std(h)))

    return np.array(temperatures, dtype=float), np.array(avg_h), np.array(std_h)


def compute_cp(
    temperatures: np.ndarray,
    enthalpy: np.ndarray,
    degree: int = 1,
) -> dict[str, object]:
    """Fit a polynomial to H(T) and extract Cp as dH/dT.

    Args:
        temperatures (np.ndarray): Temperatures in K.
        enthalpy (np.ndarray): Average enthalpy in eV/atom.
        degree (int): Polynomial degree for the H(T) fit. Degree 1 gives a
            constant Cp; degree 2 gives a temperature-dependent Cp.
            Defaults to 1.

    Returns:
        Dictionary with keys:
            - ``"coeffs"`` (np.ndarray): Polynomial coefficients ``[c0, c1, ...]``.
            - ``"cp_ev_per_k"`` (float): Linear Cp coefficient in eV/(atom·K).
            - ``"r2"`` (float): Coefficient of determination.
            - ``"y_fit"`` (np.ndarray): Fitted enthalpy values at the input temperatures.
    """
    coeffs = poly_fit(temperatures, enthalpy, degree)
    y_fit = poly_eval(coeffs, temperatures)
    r2 = r_squared(enthalpy, y_fit)
    return {
        "coeffs": coeffs,
        "cp_ev_per_k": float(coeffs[1]),
        "r2": r2,
        "y_fit": y_fit,
    }
