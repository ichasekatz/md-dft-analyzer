"""DFT ENCUT convergence analysis from VASP OUTCAR files.

Reads OUTCAR files from a directory of per-ENCUT subdirectories, extracts
total energy and pressure, and prints a convergence table.

Set ``outcar_dir`` to the root of your DFT outputs.
"""

from __future__ import annotations

from pathlib import Path

from md_dft_analyzer.dft.encut_eval import encut_convergence

# Root directory containing ENCUT subdirectories: outcar_dir/300/OUTCAR etc.
outcar_dir = Path("dft_outputs/Pt")
encut_values = [300, 400, 500, 600, 700, 800, 900]

if __name__ == "__main__":
    results = encut_convergence(outcar_dir, encut_values)

    print(f"{'ENCUT (eV)':<12} {'Energy (eV)':<18} {'Pressure (kBar)':<15}")
    for r in results:
        e = r["total_energy_ev"]
        p = r["pressure_kbar"]
        e_str = f"{e:.6f}" if e is not None else "N/A"
        p_str = f"{p:.2f}" if p is not None else "N/A"
        print(f"{r['encut_ev']:<12} {e_str:<18} {p_str:<15}")
