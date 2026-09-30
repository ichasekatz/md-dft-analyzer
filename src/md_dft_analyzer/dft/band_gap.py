"""Band gap and band structure extraction from VASP EIGENVAL files.

Parses the VASP EIGENVAL file to extract k-points, eigenvalues, and
the indirect/direct band gap.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def parse_eigenval(eigenval_path: Path | str) -> dict[str, object]:
    """Parse a VASP EIGENVAL file and return k-points and eigenvalues.

    Args:
        eigenval_path (Path | str): Path to the EIGENVAL file.

    Returns:
        Dictionary with keys:
            - ``"n_kpoints"`` (int): Number of k-points.
            - ``"n_bands"`` (int): Number of bands.
            - ``"kpoints"`` (list[list[float]]): k-point fractional coordinates.
            - ``"eigenvalues"`` (np.ndarray): Shape ``(n_kpoints, n_bands)`` in eV.

    Raises:
        FileNotFoundError: If the EIGENVAL file does not exist.
    """
    eigenval_path = Path(eigenval_path)
    if not eigenval_path.exists():
        raise FileNotFoundError(eigenval_path)

    with eigenval_path.open() as fh:
        data = fh.readlines()

    n_kpts = int(data[5].split()[1])
    n_bands = int(data[5].split()[2])

    kpoints: list[list[float]] = []
    eigenvalues: list[list[float]] = []

    for ki in range(n_kpts):
        line_idx = 7 + (n_bands + 2) * ki
        kpoints.append([float(x) for x in data[line_idx].split()[:3]])
        evals = [float(data[line_idx + bi + 1].split()[1]) for bi in range(n_bands)]
        eigenvalues.append(evals)

    return {
        "n_kpoints": n_kpts,
        "n_bands": n_bands,
        "kpoints": kpoints,
        "eigenvalues": np.array(eigenvalues),
    }


def compute_band_gap(eigenval_path: Path | str, n_valence_bands: int) -> dict[str, float]:
    """Compute the band gap from a VASP EIGENVAL file.

    Assumes spin-unpolarized calculation. The valence band maximum (VBM) and
    conduction band minimum (CBM) are found across all k-points.

    Args:
        eigenval_path (Path | str): Path to the EIGENVAL file.
        n_valence_bands (int): Number of occupied (valence) bands. Typically
            half the total number of electrons.

    Returns:
        Dictionary with keys:
            - ``"vbm_ev"`` (float): Valence band maximum in eV.
            - ``"cbm_ev"`` (float): Conduction band minimum in eV.
            - ``"gap_ev"`` (float): Band gap in eV (0.0 if metallic).
            - ``"is_metal"`` (bool): True if VBM > CBM (band overlap).

    Raises:
        FileNotFoundError: If the EIGENVAL file does not exist.
    """
    parsed = parse_eigenval(eigenval_path)
    eigenvalues: np.ndarray = parsed["eigenvalues"]  # type: ignore[assignment]

    vbm = float(np.max(eigenvalues[:, n_valence_bands - 1]))
    cbm = float(np.min(eigenvalues[:, n_valence_bands]))
    gap = max(0.0, cbm - vbm)

    return {"vbm_ev": vbm, "cbm_ev": cbm, "gap_ev": gap, "is_metal": cbm < vbm}
