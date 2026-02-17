---
description: Setup and run NH3 GCMC simulations for MOF screening
---

# NH3 GCMC Simulation Setup

## 1. Install RASPA2 (GCMC Engine)

// turbo
```bash
conda install -c conda-forge raspa2
```

Verify:
// turbo
```bash
simulate -h
```

## 2. Install PACMOF2 (ML DDEC6 Charges)

```bash
git clone https://github.com/snurr-group/PACMOF2.git
cd PACMOF2
pip install -e .
```

## 3. Install Python Dependencies

// turbo
```bash
pip install pymatgen ase numpy openbabel-wheel python-dotenv p_tqdm
```

## 4. Configure Environment

```bash
cd d:\FINAL_YEAR\GCMC
copy .env.example .env
```

Edit `.env` and set the paths:
- `RASPA_PATH`: path to the RASPA2 installation (e.g. `C:\Users\SIBI\miniconda3\envs\gcmc\share\raspa`)
- `RASPA_SIM_PATH`: path to `simulate` binary (usually just `simulate` if in PATH)

Find RASPA path with:
// turbo
```bash
python -c "import RASPA2; print(RASPA2.__file__)"
```

## 5. Run Validation (Pure NH3 against reference dataset)

```bash
cd d:\FINAL_YEAR\GCMC
python validate_nh3_uptake.py --n_mofs 3 --charge_method pacmof2
```

## 6. Run NH3 Air Removal Screening

```bash
cd d:\FINAL_YEAR\GCMC
python gcmc_screen_nh3.py --input ./path/to/cifs --mode air_removal --max_mofs 5
```

## 7. Run Pure NH3 Screening

```bash
cd d:\FINAL_YEAR\GCMC
python gcmc_screen_nh3.py --input ./raw_nh3_core/cifs --mode pure --max_mofs 5
```
