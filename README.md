<div align="center">

# md-dft-analyzer

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)

**Analysis tools for MD (LAMMPS) and DFT (VASP) simulation outputs.**

</div>

## Overview

`md-dft-analyzer` extracts thermodynamic and electronic properties from LAMMPS and VASP output files.

**MD analysis (`md_dft_analyzer.md`)**:
- `cp_eval` — constant-pressure heat capacity (Cp) from NPT enthalpy runs
- `rdf_eval` — radial distribution function (RDF) from LAMMPS `compute rdf` output

**DFT analysis (`md_dft_analyzer.dft`)**:
- `encut_eval` — ENCUT convergence: total energy and pressure from OUTCAR
- `band_gap` — band gap and band structure from EIGENVAL

**Shared utilities (`md_dft_analyzer.numerics`)**:
- `poly_fit`, `poly_eval`, `gauss_pivot`, `r_squared` — polynomial least-squares fitting

Consolidates the `MD-Toy-Codes` `Evaluation_Codes/` directory and the DFT `Evaluation_Codes/` from `DFT-Toy-Codes`.

## Installation

```bash
git clone https://github.com/ichasekatz/md-dft-analyzer.git
cd md-dft-analyzer
uv sync
```

## Quick Start

```python
# MD: heat capacity from LAMMPS output
from md_dft_analyzer.md.cp_eval import load_enthalpy_series, compute_cp

temps, avg_h, std_h = load_enthalpy_series("outputs/samp_{T}.out", [200, 400, 600, 800, 1000])
result = compute_cp(temps, avg_h)
print(f"Cp = {result['cp_ev_per_k']:.4f} eV/(atom·K)")

# DFT: ENCUT convergence from OUTCAR
from md_dft_analyzer.dft.encut_eval import encut_convergence

results = encut_convergence("dft_outputs/Pt", encut_values=[400, 500, 600, 700])
```

## Module Reference

| Module | Key functions |
|--------|--------------|
| `md_dft_analyzer.md.cp_eval` | `load_enthalpy_series`, `compute_cp` |
| `md_dft_analyzer.md.rdf_eval` | `load_rdf` |
| `md_dft_analyzer.dft.encut_eval` | `parse_outcar`, `encut_convergence` |
| `md_dft_analyzer.dft.band_gap` | `parse_eigenval`, `compute_band_gap` |
| `md_dft_analyzer.numerics` | `poly_fit`, `poly_eval`, `r_squared` |

## Examples

| Script | Description |
|--------|-------------|
| `examples/01_md_cp_eval.py` | Cp from LAMMPS NPT runs |
| `examples/02_dft_encut_eval.py` | ENCUT convergence table from OUTCAR files |

```bash
uv run python examples/02_dft_encut_eval.py
```

## Running Tests

```bash
uv run pytest -v
```

## License

GPL-3.0-or-later — see [LICENSE](LICENSE).

## Citation

```bibtex
@software{katz_md_dft_analyzer_2026,
  author = {Katz, Chase},
  title  = {md-dft-analyzer: MD and DFT simulation output analysis},
  year   = {2026},
  url    = {https://github.com/ichasekatz/md-dft-analyzer},
}
```
