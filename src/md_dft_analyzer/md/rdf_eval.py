"""Radial distribution function (RDF) parsing and visualization from LAMMPS output.

Reads a LAMMPS ``compute rdf`` output file (3-column format: N, r, g(r))
and returns per-frame DataFrames for downstream analysis or plotting.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_rdf(
    rdf_path: Path | str,
    n_bins: int = 200,
    skiprows: int = 3,
) -> list[pd.DataFrame]:
    """Parse a LAMMPS multi-frame RDF output file.

    Each frame in the file contains ``n_bins`` rows with columns
    ``[frame_index, r (Å), g(r)]``. Rows are split into per-frame DataFrames.

    Args:
        rdf_path (Path | str): Path to the LAMMPS ``.rdf`` file.
        n_bins (int): Number of histogram bins per frame. Defaults to 200.
        skiprows (int): Header lines to skip at the top of the file.
            Defaults to 3 (LAMMPS ``fix ave/time`` header).

    Returns:
        List of DataFrames, one per RDF frame, each with columns
        ``["frame", "r", "g_r"]``.

    Raises:
        FileNotFoundError: If ``rdf_path`` does not exist.
    """
    rdf_path = Path(rdf_path)
    if not rdf_path.exists():
        raise FileNotFoundError(rdf_path)

    df_raw = pd.read_csv(
        rdf_path,
        skiprows=skiprows,
        sep=r"\s+",
        names=["frame", "r", "g_r", "coord"],
        on_bad_lines="skip",
    )
    df_raw = df_raw.apply(pd.to_numeric, errors="coerce").dropna()
    df_valid = df_raw[(df_raw["frame"] >= 1) & (df_raw["frame"] <= n_bins)]

    frames: list[pd.DataFrame] = []
    for _, group in df_valid.groupby(df_valid.index // n_bins):
        frames.append(group[["frame", "r", "g_r"]].reset_index(drop=True))
    return frames
