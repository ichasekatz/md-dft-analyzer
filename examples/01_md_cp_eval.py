"""MD heat capacity (Cp) evaluation from LAMMPS output.

Reads LAMMPS NPT enthalpy output files at a series of temperatures,
fits a polynomial H(T), and extracts Cp = dH/dT.

Set ``file_pattern`` and ``temperatures`` to match your simulation layout.
"""

from __future__ import annotations

from md_dft_analyzer.md.cp_eval import compute_cp, load_enthalpy_series

# Update this path to point to your LAMMPS sample files
# Expected: one file per temperature, named samp_{T}.out
file_pattern = "outputs/samp_{T}.out"
temperatures = list(range(200, 1300, 100))  # 200, 300, ..., 1200 K

if __name__ == "__main__":
    temps, avg_h, std_h = load_enthalpy_series(file_pattern, temperatures)
    result = compute_cp(temps, avg_h, degree=1)

    print(f"Cp (linear fit):  {result['cp_ev_per_k']:.6f} eV/(atom·K)")
    print(f"R²:               {result['r2']:.6f}")
