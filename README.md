# NH3-GCMC — Ammonia Capture Screening via GCMC Simulations

[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![RASPA2](https://img.shields.io/badge/engine-RASPA2-green.svg)](https://github.com/iRASPA/RASPA2)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A computational pipeline for screening Metal–Organic Frameworks (MOFs) for **ammonia (NH₃) capture** using Grand Canonical Monte Carlo (GCMC) simulations. Built on top of [RASPA2](https://github.com/iRASPA/RASPA2) and adapted from the [MOFDiff](https://github.com/microsoft/MOFDiff) GCMC workflow for CO₂.

> **Reference Paper:** *"Computational Screening of MOFs for Ammonia Capture from H₂/N₂/NH₃ Mixtures"* — [ACS Omega 2022](https://doi.org/10.1021/acsomega.2c04517)

---

## ✨ Features

- **Pure NH₃ adsorption** — Single-component GCMC at 1 bar / 298 K for benchmarking
- **3-component air removal** — NH₃/H₂O/N₂ mixture with Vacuum Swing Adsorption (VSA) working capacity
- **Automated charge assignment** — PACMOF2 (ML-predicted DDEC6) or MEPO-Qeq (eGULP) partial charges
- **Helium void fraction** — Widom insertion for accurate uptake normalization
- **Batch screening** — Screen directories of CIF files with optional parallel processing
- **Built-in validation** — Compare simulated uptakes against the `raw_nh3_core` reference dataset

---

## 📁 Project Structure

```
GCMC/
├── gcmc_wrapper_nh3.py         # Core RASPA2 wrapper (charge calc, void fraction, GCMC)
├── nh3_simulation.py           # High-level simulation scenarios (pure / air removal)
├── gcmc_screen_nh3.py          # Batch screening entry point
├── validate_nh3_uptake.py      # Validation against raw_nh3_core dataset
├── force_field/                # Local force field files for RASPA2
│   ├── force_field.def         # Force field master definition
│   ├── force_field_mixing_rules.def  # Lorentz-Berthelot mixing rules
│   └── pseudo_atoms.def        # Framework atom LJ parameters
├── molecule_definitions/       # Sorbate molecule models
│   └── TraPPE/
│       ├── Ammonia.def         # 4-site TraPPE-EH ammonia model
│       └── H2O.def             # Water model for mixture simulations
├── raw_nh3_core/               # Validation dataset (CIFs + expected uptakes)
├── .env.example                # Environment variable template
└── NH3_gcmc.txt                # Technical notes (CO₂ → NH₃ adaptation)
```

---

## 🚀 Quick Start

### 1. Prerequisites

| Dependency | Purpose | Install |
|---|---|---|
| [RASPA2](https://github.com/iRASPA/RASPA2) | GCMC simulation engine | `conda install -c conda-forge raspa2` |
| [PACMOF2](https://github.com/snurr-group/PACMOF2) | ML-predicted DDEC6 charges | `pip install -e .` (from source) |
| [eGULP](https://github.com/uowoolab/eGULP) *(optional)* | MEPO-Qeq charges (fallback) | Build from source |
| [p_tqdm](https://github.com/swansonk14/p_tqdm) *(optional)* | Parallel screening | `pip install p_tqdm` |
| NumPy | Numerical utilities | `pip install numpy` |

### 2. Environment Setup

```bash
# Clone the repository
git clone https://github.com/<your-username>/NH3-GCMC.git
cd NH3-GCMC

# Copy and configure environment variables
cp .env.example .env
```

Edit `.env` with the paths for your system:

```env
# RASPA2 — find path with: python -c "import RASPA2; print(RASPA2.__file__)"
RASPA_PATH=/path/to/raspa2
RASPA_SIM_PATH=simulate

# PACMOF2 (required for default charge method)
PACMOF2_PATH=/path/to/pacmof2

# eGULP (optional — only if using mepo_qeq charge method)
EGULP_PATH=/path/to/egulp
EGULP_PARAMETER_PATH=/path/to/egulp/parameters
```

### 3. Run a Single MOF

```python
from nh3_simulation import nh3_uptake_pure

result = nh3_uptake_pure(
    "path/to/structure.cif",
    calc_charges=True,
    charge_method="pacmof2",   # or "mepo_qeq"
    temperature=298,            # K
    pressure=100000,            # 1 bar in Pa
)

print(f"NH3 uptake: {result['NH3_uptake_mmol_g']:.4f} mmol/g")
print(f"Qst: {result['NH3_heat_of_adsorption_kJ_mol']:.2f} kJ/mol")
```

### 4. Batch Screening

```bash
# Pure NH₃ adsorption (validation mode)
python gcmc_screen_nh3.py --input ./raw_nh3_core/cifs --mode pure --max_mofs 5

# 3-component air removal with vacuum swing
python gcmc_screen_nh3.py --input ./path/to/cifs --mode air_removal --max_mofs 10

# Parallel screening with MEPO-Qeq charges
python gcmc_screen_nh3.py --input ./path/to/cifs --charge_method mepo_qeq --parallel --ncpu 8
```

### 5. Validate Against Reference Data

```bash
python validate_nh3_uptake.py --n_mofs 5 --charge_method pacmof2
```

This runs GCMC on a diverse subset (low, mid, high uptake) from `raw_nh3_core` and reports per-MOF error percentages.

---

## 🔬 Simulation Details

### Simulation Modes

| Mode | Conditions | Use Case |
|---|---|---|
| **`pure`** | 100% NH₃, 1 bar, 298 K | Validation & benchmarking |
| **`air_removal`** | NH₃:H₂O:N₂ = 0.1%:3%:96.9% @ 298 K, 1 bar → desorption @ 400 K, 0.01 bar | Realistic air capture / VSA |

### Key Parameters

| Parameter | Value | Rationale |
|---|---|---|
| Force field | TraPPE-EH | Industry standard for polar molecules like NH₃ |
| Mixing rules | Lorentz-Berthelot | Standard combining rules |
| Cutoff | 12 Å | Long-range interactions for polar NH₃ |
| Eq / Prod cycles | 20,000 / 20,000 | Paper-validated convergence |
| Swap probability | 2.0 | Higher insertion rate for polar molecules in tight pores |
| Charge method | PACMOF2 (default) | ML-predicted DDEC6 — DFT-quality at classical speed |

### Pipeline Workflow

```
CIF File
  │
  ├── 1. Charge Assignment (PACMOF2 / MEPO-Qeq)
  │
  ├── 2. Helium Void Fraction (Widom insertion)
  │
  ├── 3. Unit Cell Replication (minimum image convention)
  │
  └── 4. GCMC Simulation (RASPA2)
         │
         └── Output: uptake (mmol/g), Qst (kJ/mol), selectivity
```

---

## 📊 Output Format

Results are saved as JSON. Example for `pure` mode:

```json
{
  "NH3_uptake_mol_kg": 5.123,
  "NH3_uptake_mmol_g": 5.123,
  "NH3_heat_of_adsorption_kJ_mol": -42.5
}
```

For `air_removal` mode:

```json
{
  "adsorption_NH3_mmol_g": 3.21,
  "desorption_NH3_mmol_g": 0.15,
  "working_capacity_nh3_mmol_g": 3.06,
  "NH3_N2_selectivity": 1250.3
}
```

---

## 🧪 Validation Dataset

The `raw_nh3_core/` directory contains ~1000 MOF CIF files with experimentally-derived NH₃ uptake values at 298 K / 1 bar. Use `validate_nh3_uptake.py` to benchmark simulation accuracy against this dataset.

---

## 📖 Citation

If you use this pipeline, please cite the reference paper:

```bibtex
@article{nh3_mof_screening_2022,
  title={Computational Screening of Metal--Organic Frameworks for Ammonia Capture from H2/N2/NH3 Mixtures},
  journal={ACS Omega},
  year={2022},
  doi={10.1021/acsomega.2c04517}
}
```

---

## 🤝 Acknowledgements

- **[RASPA2](https://github.com/iRASPA/RASPA2)** — GCMC simulation engine  
- **[MOFDiff](https://github.com/microsoft/MOFDiff)** — Original CO₂ GCMC pipeline (adapted here for NH₃)  
- **[PACMOF2](https://github.com/snurr-group/PACMOF2)** — ML-predicted DDEC6 partial charges  
- **[TraPPE Force Field](http://trappe.oit.umn.edu/)** — Transferable potentials for NH₃ and H₂O  

---

## 📝 License

This project is licensed under the [MIT License](LICENSE).
