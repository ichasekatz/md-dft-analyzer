"""ENCUT convergence analysis from VASP OUTCAR files.

Reads total energy and pressure from a set of OUTCAR files produced at
different ENCUT values and plots convergence curves.
"""

from __future__ import annotations

import re
from pathlib import Path


_energy_pattern = re.compile(r"free  energy   TOTEN\s+=\s+([-+]?\d+\.\d+)")
_pressure_pattern = re.compile(r"external pressure\s*=\s*([-+]?\d+\.\d+)")


def parse_outcar(outcar_path: Path | str) -> dict[str, float | None]:
    """Extract total energy and pressure from a VASP OUTCAR file.

    Reads until both values are found in the OUTCAR, stopping early
    for efficiency.

    Args:
        outcar_path (Path | str): Path to the OUTCAR file.

    Returns:
        Dictionary with keys ``"total_energy_ev"`` (eV) and
        ``"pressure_kbar"`` (kBar). Values are ``None`` if not found.

    Raises:
        FileNotFoundError: If the OUTCAR file does not exist.
    """
    outcar_path = Path(outcar_path)
    if not outcar_path.exists():
        raise FileNotFoundError(outcar_path)

    total_energy: float | None = None
    pressure: float | None = None

    with outcar_path.open() as fh:
        for line in fh:
            if total_energy is None:
                m = _energy_pattern.search(line)
                if m:
                    total_energy = float(m.group(1))
            if pressure is None:
                m = _pressure_pattern.search(line)
                if m:
                    pressure = float(m.group(1))
            if total_energy is not None and pressure is not None:
                break

    return {"total_energy_ev": total_energy, "pressure_kbar": pressure}


def encut_convergence(
    outcar_dir: Path | str,
    encut_values: list[int],
    outcar_subdir: str = "{encut}",
) -> list[dict[str, object]]:
    """Collect ENCUT convergence data from a directory of OUTCAR files.

    Expects one OUTCAR per ENCUT value at ``outcar_dir/{encut}/OUTCAR``
    (configurable via ``outcar_subdir``).

    Args:
        outcar_dir (Path | str): Root directory containing per-ENCUT subdirs.
        encut_values (list[int]): ENCUT values in eV to scan.
        outcar_subdir (str): Subdirectory template with ``{encut}`` placeholder.
            Defaults to ``"{encut}"`` (e.g. ``outcar_dir/300/OUTCAR``).

    Returns:
        List of dicts with keys ``"encut_ev"``, ``"total_energy_ev"``,
        ``"pressure_kbar"``. Missing OUTCARs are skipped with a warning.
    """
    outcar_dir = Path(outcar_dir)
    results: list[dict[str, object]] = []

    for encut in encut_values:
        path = outcar_dir / outcar_subdir.format(encut=encut) / "OUTCAR"
        if not path.exists():
            print(f"Warning: {path} not found — skipping ENCUT={encut}")
            continue
        parsed = parse_outcar(path)
        results.append({"encut_ev": encut, **parsed})

    return results
