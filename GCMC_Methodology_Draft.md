# Grand Canonical Monte Carlo (GCMC) Simulation Methodology
## Full Technical Specification — NH₃ Capture Screening

Compiled from actual implementation files: `gcmc_wrapper_nh3.py`, `nh3_simulation.py`,
`force_field_mixing_rules.def`, `pseudo_atoms.def`, `Ammonia.def`, `H2O.def`

---

## 1. Goal and Target Metrics

### Task
Two simulation scenarios are evaluated:

| Scenario | Description |
|---|---|
| **Pure NH₃** | Benchmark/validation against CoRE MOF reference data |
| **VSA Air Removal** | 3-component NH₃/H₂O/N₂ humid air capture with vacuum-swing regeneration |

### Exact Operating Points

| State | T (K) | P (Pa / bar) | Composition |
|---|---|---|---|
| Adsorption (pure) | 298 | 100,000 Pa (1.0 bar) | y(NH₃) = 1.0 |
| Adsorption (VSA) | 298 | 100,000 Pa (1.0 bar) | y(NH₃)=0.001, y(H₂O)=0.03, y(N₂)=0.969 |
| Desorption (VSA) | 400 | 1,000 Pa (0.01 bar) | y(NH₃) = 1.0 |

> **Note on NH₃ concentration:** 1,000 ppm NH₃ (y = 0.001) represents a realistic indoor/agricultural air scenario. The humid-air VSA scenario (NH₃/H₂O/N₂) is an original contribution of this work; for a comparable humid-air NH₃ capture study see, e.g., Brandt *et al.* and related DAC/air-purification literature. The ACS Omega 2022 paper (DOI: 10.1021/acsomega.2c04517) addresses H₂/N₂/NH₃ *dry* mixtures (NH₃ = 1.0, 0.35, 0.2 mol%) with no water; it is cited only for TraPPE NH₃/N₂ force-field practice and general RASPA screening methodology, **not** for the operating conditions of our VSA scenario.

### Primary Labels and Units

| Label | Formula | Units |
|---|---|---|
| **Absolute NH₃ uptake** | Direct RASPA output | **mmol/g** (= mol/kg) |
| **VSA Working Capacity** | WC = N_ads − N_des | **mmol/g** |
| **NH₃/N₂ Selectivity** | S = (x_NH₃/y_NH₃)/(x_N₂/y_N₂) | dimensionless |
| **NH₃/H₂O Selectivity** | S = (x_NH₃/y_NH₃)/(x_H₂O/y_H₂O) | dimensionless |
| **Heat of Adsorption** | Q_st = −H_ads | **kJ/mol** |

**All uptake quantities are reported as absolute loading in mol/kg (equivalently mmol/g).** Since the molar mass of 1 kg of framework is the denominator, mol/kg = mmol/g numerically.

### Gas Phase Compositions (VSA only — 1 case)

| Component | Mole Fraction |
|---|---|
| NH₃ | 0.001 (1,000 ppm) |
| H₂O | 0.030 (3%) |
| N₂ | 0.969 (balance) |

### H₂O Treatment
H₂O is explicitly included in the gas-phase mixture during adsorption GCMC (co-adsorption, not sequential). For desorption, pure NH₃ at 400 K / 0.01 bar is used (H₂O omitted, representing complete water desorption during regeneration). This is a deliberate simplification that must be acknowledged as a limitation.

### KPI for Mixture
- **Primary:** VSA working capacity (mmol/g NH₃)
- **Secondary:** NH₃/N₂ and NH₃/H₂O adsorption selectivities from the mixture simulation
- **Selectivity formula** (adsorption-phase mole-fraction ratio):

  S(NH₃/X) = [x_NH₃ / y_NH₃] / [x_X / y_X]

  where x = adsorbed loading (mol/kg) and y = gas-phase mole fraction. Note: this is the *loading-based separation factor*, not IAST selectivity.

---

## 2. MOF Dataset and Preprocessing

### Source
**CoRE MOF 2019** dataset — CIF files from `raw_nh3_core/` directory (1,030 structures confirmed present).

### Filters Applied
The following filters are applied **before** simulation (consistent with standard CoRE MOF screening workflows):

1. **Disorder removal:** Only ordered structures included (CIF files must resolve disorder).
2. **Missing hydrogen:** Structures with missing H atoms on organic linkers are excluded.
3. **Charge neutrality:** After partial charges are assigned by PACMOF2 or MEPO-Qeq, the net charge of the unit cell is verified explicitly using `verify_charge_neutrality()` in `gcmc_wrapper_nh3.py` (tolerance: |q_net| ≤ 0.05 e). Structures failing this check are excluded from screening and logged.
4. **Geometric viability:** Pore-limiting diameter (PLD) > kinetic diameter of NH₃ (3.62 Å). Structures with PLD < 3.62 Å are excluded as NH₃ cannot diffuse through them.

### Open Metal Sites
CIF files are used **as-is** from the CoRE MOF 2019 release (already solvent-removed / activated). Open metal sites are not explicitly blocked or specially treated — their GCMC interaction is determined by framework atom LJ parameters and PACMOF2 charges. This is a known limitation for OMS-containing MOFs.

### Geometric Pre-screening
- PLD and ASA are evaluated using **Zeo++ with a helium-probe radius of 1.2 Å** (standard for CoRE MOF geometric characterization). The PLD exclusion threshold is the NH₃ kinetic diameter (3.62 Å); this is intentionally a **coarse pre-filter** — a He-probe PLD is a geometric lower bound that will pass some structures inaccessible to NH₃. The actual molecular accessibility is determined by GCMC acceptance during simulation: if no NH₃ molecules are inserted (loading = 0), the structure is effectively inaccessible.
- For a rigorous NH₃-size filter, Zeo++ should be re-run with an NH₃ probe radius of 1.81 Å (r = d_kin/2 = 3.62/2). This is noted as a recommended improvement.

### Duplicate Handling
Duplicates are not explicitly handled in the current pipeline; the CoRE MOF 2019 dataset itself curates for crystallographic uniqueness by refcode.

### Framework Relaxation
**No geometry relaxation.** CIF files are used as-is. This is standard practice in high-throughput GCMC screening; DFT relaxation is reserved for top candidates only.

---

## 3. Simulation Engine and Reproducibility

### Code
**RASPA2** — `simulate` binary invoked via `subprocess` in `gcmc_wrapper_nh3.py`.
- The RASPA2 path is configured via the environment variable `RASPA_PATH` and `RASPA_SIM_PATH` (see `.env.example`).
- Target: RASPA2 v2.0.x (exact commit hash depends on cloud instance build; should be documented upon final runs).

### Input File Publication
All simulation inputs are reproduced programmatically by `gcmc_wrapper_nh3.py`. Each run writes:
- `simulation.input` — full RASPA2 input
- `force_field.def`, `force_field_mixing_rules.def`, `pseudo_atoms.def` — local force field
- `molecules/TraPPE/Ammonia.def`, `H2O.def` — molecule definitions

These files will be published with the final dataset.

### Random Seed
Fixed seeds: `random.seed(4)` and `np.random.seed(4)` (set in `nh3_simulation.py`). A single seed is used per MOF per condition.

---

## 4. Force Field Choices

### MOF Framework
Framework atoms interact via Lennard-Jones parameters taken from the **Universal Force Field (UFF)**. The `force_field.def` file specifies 0 custom rules and 0 overwritten rules, meaning RASPA2 looks up all framework atom types from its internal UFF library at runtime.

### NH₃ — TraPPE-EH (Transferable Potentials for Phase Equilibria — Explicit Hydrogen)
A **5-site rigid** model (from `Ammonia.def`):

| Site | ε/k_B (K) | σ (Å) | Charge (e) | Notes |
|---|---|---|---|---|
| N_nh3 | 185.0 | 3.42 | 0.0 | from `force_field_mixing_rules.def` |
| H_nh3 (×3) | 0.0 | 0.0 | +0.41 | explicit H; LJ off |
| com_nh3 | 0.0 | 0.0 | −1.23 | center-of-mass charge site |

Geometry from `Ammonia.def`: N at origin, three explicit H at positions (±0.4714, ±0.8165, 0.3351) Å (pyramidal C₃ᵥ symmetry). The molecule is **rigid**.

> Note: The `pseudo_atoms.def` lists charges on H_nh3 = +0.41 e and com_nh3 = −1.23 e, giving a net dipole moment that accurately captures NH₃ polarity. This is the standard TraPPE-EH model (Eckl et al. / Lorentz-Berthelot).

### H₂O — TIP4P/2005
A **4-site rigid** model (from `H2O.def`):

| Site | ε/k_B (K) | σ (Å) | Charge (e) |
|---|---|---|---|
| O_w | 93.2 | 3.1589 | 0.0 |
| H_w (×2) | 0.0 | 0.0 | +0.5564 |
| M_w | 0.0 | 0.0 | −1.1128 |

### N₂ — TraPPE (3-site, Potoff & Siepmann 2001)
Defined in `molecule_definitions/TraPPE/N2.def` (local file, not RASPA internal library):

| Site | ε/k_B (K) | σ (Å) | Charge (e) | Position (z, Å) |
|---|---|---|---|---|
| N_n2 (×2) | 36.0 | 3.31 | −0.482 | ±0.55 |
| com_n2 | 0.0 | 0.0 | +0.964 | 0.0 (COM) |

Bond length 1.10 Å. The COM site carries a neutralizing charge. Net molecular charge = 2×(−0.482) + 0.964 = **0.0 e** ✓
The `com_n2` entry is registered in `pseudo_atoms.def` (8 total sites) and `force_field_mixing_rules.def` (8 entries).

### Combining Rules
**Lorentz–Berthelot**: ε_ij = √(ε_i · ε_j), σ_ij = (σ_i + σ_j)/2

Explicitly set as `Lorentz-Berthelot` in `force_field_mixing_rules.def`.

### Cross-interaction Tuning
**None.** No special NH₃–metal-site cross terms are added. All cross-interactions are derived from LB mixing rules applied to the TraPPE-EH NH₃ LJ parameters and UFF framework parameters.

---

## 5. Potential Truncation

### Shifting vs. Truncation
From `force_field_mixing_rules.def` line 2: **`shifted`**. The LJ potential is shifted to zero at the cutoff. This eliminates the discontinuity at the cutoff but removes tail corrections.

### Tail Corrections
`force_field_mixing_rules.def` line 4: **`no`**. Long-range tail corrections are **not applied**.

> **Important note for reviewers:** Disabling tail corrections with a shifted potential is internally consistent (tail corrections assume an unshifted truncated potential). However, this is known to underestimate vaporization enthalpies by ~2–5%. For comparative high-throughput screening this is acceptable, but should be noted explicitly.

### Cutoffs (from `simulation.input` generated by `gcmc_wrapper_nh3.py`)
- `CutOffVDW`: **12.0 Å**
- `CutOffChargeCharge`: **12.0 Å**
- `CutOffChargeBondDipole`: **12.0 Å**
- `CutOffBondDipoleBondDipole`: **12.0 Å**

---

## 6. Electrostatics and Partial Charges

### Framework Charge Assignment
- **Primary method:** **PACMOF2** — machine-learning model trained to reproduce DDEC6 charges. Predicts near DFT-quality charges at negligible cost. Implemented via `_calculate_pacmof2_charges()` in `gcmc_wrapper_nh3.py` using `pymatgen` + `pacmof` package.
- **Fallback / comparison method:** **MEPO-Qeq** via eGULP, implemented in `_calculate_mepo_qeq_charges()`.

### PACMOF2 Validation
PACMOF2 was validated by the Snurr group against DFT-DDEC6 charges on a diverse MOF set. For NH₃ specifically, the sensitivity of uptake to charge method is quantified using `validate_charges.py`:

```bash
python validate_charges.py --input ./raw_nh3_core/cifs --max_mofs 10 --threshold 0.5
```

This script runs pure NH₃ GCMC at 298 K / 1 bar for a 10-MOF subset with both PACMOF2 and MEPO-Qeq, prints a side-by-side comparison table, and flags MOFs where |ΔUPTAKE| > 0.5 mmol/g. Results are saved to `gcmc_runs/charge_validation_results.json`.

### Adsorbate Charges
Fixed — from the TraPPE-EH and TIP4P/2005 model definitions. Not scaled.

### Electrostatics Method
**Ewald summation** — `ChargeMethod Ewald`, `EwaldPrecision 1e-6` (set in `simulation.input`).

### Charge Neutrality
RASPA2 handles charge neutrality internally when `UseChargesFromCIFFile yes` is set. The CIF charges from PACMOF2 are used directly; PACMOF2 enforces charge neutrality during prediction.

---

## 7. Framework Treatment

**Rigid framework.** All MOFs are treated as rigid during GCMC sampling. This is the standard approximation for high-throughput screening.

**Justification for NH₃:** While NH₃ is known to induce framework breathing in some flexible MOFs (e.g., MIL-53 family), the current screening covers broad structural diversity. Flexible GCMC would require per-structure force field parameterization and is computationally prohibitive at scale. Rigidity is a recognized limitation for known flexible MOFs.

---

## 8. Simulation Cell and Boundary Conditions

### Supercell Construction
Computed automatically in `calculate_unit_cells()` (`gcmc_wrapper_nh3.py`, lines 115–128):

```
unit_cells[i] = ceil(2 × cutoff / perpendicular_length[i])
```

This enforces the **minimum image convention**: the simulation box width in each direction ≥ 2 × 12 Å = **24 Å**.

### Boundary Conditions
3D periodic boundary conditions enforced by RASPA2.

### Handling Low-Symmetry Cells
The supercell is chosen so that the **shortest perpendicular cell dimension** of the replicated box exceeds 2 × cutoff (24 Å). For a general triclinic cell with lattice vectors a, b, c and angles α, β, γ, the perpendicular widths are:

```
d_a = a · √det / (sin β · sin γ)
d_b = b · √det / (sin α · sin γ)        where det = 1 − cos²α − cos²β − cos²γ + 2·cosα·cosβ·cosγ
d_c = c · √det / (sin α · sin β)
```

This is the exact formula derived from the volume V = a·b·c·√det divided by the face areas. The code (`calculate_unit_cells()` in `gcmc_wrapper_nh3.py`) implements this formula, replacing the earlier approximation `dim[i] × |cos(angle[i] − 90°)|` which was only correct for pseudo-orthogonal cells.

---

## 9. GCMC Move Set and Sampling

Move probabilities (set per sorbate component, from `simulation.input`):

| Move Type | Relative Probability |
|---|---|
| Translation | 0.5 |
| Rotation | 0.5 |
| Reinsertion | 0.5 |
| Insertion/Deletion (Swap) | **2.0** |
| Identity Change (multi-component) | 1.0 |
| CBMC | Not used (rigid molecules) |

**Swap probability is intentionally elevated to 2.0** to improve acceptance of NH₃ insertion into tight pores at low gas-phase mole fractions.

### CBMC
Not used. All adsorbate molecules (NH₃, H₂O, N₂) are rigid, so configurational-bias MC is not required for bond/torsion sampling. Standard random insertion is used.

### Biased Insertion
No site-biased insertion. Uniform random insertion is used.

---

## 10. Equilibration, Production, and Convergence

| Parameter | Value |
|---|---|
| Initialization cycles | 0 |
| Equilibration cycles | **20,000** |
| Production cycles | **20,000** |
| Print every | 1,000 cycles |

### Definition of a Cycle
In RASPA2, one cycle = N_molecules attempted MC moves, where N_molecules is the current number of molecules in the system (minimum 20). For GCMC this increases as loading grows.

### Convergence Assessment
Convergence is monitored via RASPA's printed running averages (output every 1,000 cycles). A flat running average in the second half of production is taken as convergence. Block SE is now computed from per-checkpoint history using `block_standard_error()` in `nh3_simulation.py` (Flyvbjerg & Petersen method, 5 blocks over the production half).

### Independent Replicas
**1 replica per MOF per condition** (fixed seed = 4). Block SE from a single run is a lower bound on true SE — for top candidates, ≥3 independent seeds are recommended.

### Acceptance at Low y(NH₃)
At y(NH₃) = 0.001, NH₃ insertion acceptance can be very low, especially in dense frameworks. The elevated swap probability (2.0) partially mitigates this. For MOFs with near-zero acceptance at this concentration, the result should be flagged and optionally re-simulated with more cycles or a pre-equilibration in pure NH₃ atmosphere.

### Tail Corrections
Not used (consistent with shifted potential — see Section 5).

---

## 11. Outputs and Derived Quantities

The following are parsed from RASPA `.data` output files by `extract_raspa_output_nh3()` in `nh3_simulation.py`:

| Quantity | RASPA Field | Units |
|---|---|---|
| Absolute NH₃ loading | `Average loading absolute [mol/kg framework]` | mol/kg |
| H₂O loading (ads) | same field for H₂O component | mol/kg |
| N₂ loading (ads) | same field for N₂ component | mol/kg |
| Enthalpy of adsorption | `Total enthalpy of adsorption [K]` | K → kJ/mol |
| Helium void fraction | `Average Widom Rosenbluth-weight` | dimensionless |

### Unit Conversions
- `mol/kg = mmol/g` (identical numerically)
- Enthalpy [K] × R (8.314 J mol⁻¹ K⁻¹) / 1000 = enthalpy [kJ/mol]
- Q_st = −H_ads (positive value = exothermic)

### Heat of Adsorption Method
**Fluctuation method** (ensemble-average enthalpy directly from RASPA). Clausius-Clapeyron multi-temperature method is not used.

### Selectivity Formula
```
S(NH₃/X) = [N_NH₃(mol/kg) / y_NH₃] / [N_X(mol/kg) / y_X]
```
Evaluated at the single adsorption composition (y_NH₃=0.001, y_H₂O=0.03, y_N₂=0.969). This is a **loading-based separation factor** (not IAST selectivity).

### Working Capacity
```
WC = N_NH₃(ads, 298K, 1bar, mixture) − N_NH₃(des, 400K, 0.01bar, pure NH₃)
```
Units: mol/kg (mmol/g).

### Uncertainty Reporting
A convergence trend indicator is computed from RASPA cumulative-average checkpoints (printed every `PrintEvery=1000` cycles). The spread (std) of equal-width window means over the production phase provides a rough convergence check: a near-zero spread indicates the average has stabilised. This is **not a statistically valid standard error** because RASPA's cumulative averages are strongly correlated by construction — the Flyvbjerg-Petersen block method requires independent or block-resampled stationary time series. All output dicts include a `*_se_mmol_g` field containing this trend spread as a convergence diagnostic. **For publication, SE must be computed from ≥3 independent seeds as: SE = std(seed_means) / √n_seeds.**

---

## 12. Helium Void Fraction

Computed via **Widom particle insertion** of helium at 298 K with 10,000 MC cycles (`compute_helium_void_fraction()`, `gcmc_wrapper_nh3.py`).

**What is being parsed:** In RASPA2, a Widom insertion simulation outputs `Average Widom Rosenbluth-weight`, which equals ⟨exp(−βu_host)⟩, the Boltzmann-weighted test-particle acceptance. For helium at 298 K, where physisorption is extremely weak (ε/k_B ≈ 10 K for He–framework), this quantity converges to the accessible geometric void fraction θ_He, as used by the Snurr group and described in Dubbeldam et al. RASPA manual (Eq. for Widom insertion). We parse this field directly as the helium void fraction and pass it to RASPA via `HeliumVoidFraction`. This is the standard approach in published GCMC screening studies. The resulting void fraction is used by RASPA to report excess uptake; we report **absolute loading** directly.

---

## 13. Humidity / Competitive Adsorption

H₂O is **explicitly included in the gas-phase mixture** during adsorption GCMC (co-adsorption). This is the thermodynamically rigorous approach — both NH₃ and H₂O compete simultaneously for adsorption sites.

H₂O model: TIP4P/2005 (4-site rigid). Charges: O_w = 0.0, H_w = +0.5564 e, M_w = −1.1128 e.

Water clustering and slow equilibration are a genuine concern: at 3% H₂O mole fraction and 1.0 bar total pressure, the water partial pressure is ~30 mbar, which is **near saturation** (saturation pressure at 298 K ≈ 31.7 mbar). This corresponds to approximately 95% relative humidity — a condition where water clustering and slow equilibration are plausible, particularly in hydrophilic frameworks and frameworks with open metal sites. The 20,000-cycle equilibration may be insufficient in such cases; this is a recognised limitation that should be stated in the paper.

---

## 14. Screening Filters and Post-Processing

Pre-simulation descriptors:
- **PLD** (pore-limiting diameter) — computed with Zeo++, probe radius = 1.2 Å
- **LCD** (largest cavity diameter)
- **ASA** (accessible surface area)
- **Pore volume**

Exclusion rule: PLD < 3.62 Å → excluded (NH₃ kinetic diameter).

No second-stage high-fidelity rerun is currently automated, but the pipeline supports it by re-calling the same functions with higher cycle counts.

---

## 15. Validation and Sanity Checks

### Validation Against Known Data
The pure NH₃ mode (`mode="pure"`) is designed explicitly for validation against the `raw_nh3_core/` reference dataset (CoRE MOF NH₃ uptake benchmark). Discrepancies > 50% on this dataset indicate a force field or charge assignment issue.

### Cross-Check: Charge Methods
The pipeline implements both PACMOF2 and MEPO-Qeq. A subset comparison (5–10 MOFs) using both methods is recommended to bound charge-method sensitivity. This is implemented via `--charge_method mepo_qeq` flag.

### Rigid vs Flexible
Not currently compared; rigid only. Limitation acknowledged.

### Cutoff Sensitivity
12 Å cutoff is standard. Sensitivity (12 vs 14 Å) on a subset has not been performed in the current pipeline. For final publication, a 5-MOF cutoff sensitivity test is recommended.

### OMS Validation
Expected behavior: MOFs with open metal sites should show high Q_st. Density maps (not currently enabled — `Movies no`) would confirm NH₃ localization near metal sites.

---

## 16. Integration with ML Pipeline

### Labels
The 1,023 labels (one per CoRE MOF, from `raw_nh3_core/`) are generated under a **single fixed condition**: 298 K, 1.0 bar, pure NH₃ (validation set). For VSA-mode labels, the condition is 298 K / 1 bar (adsorption) → 400 K / 0.01 bar (desorption).

### Predicted Scalars
The NH₃ head predicts:
- NH₃ absolute uptake (mol/kg) — primary label
- (Optionally) Q_st, working capacity, selectivity

### Label Noise from Non-Convergence
Currently detected by checking if RASPA output parse fails (exception → `info = "failed"`, `result = None`). A more robust approach is to check running averages for trends in the last 25% of production. Non-converged runs should be excluded from training or down-weighted.

### Per-Sample GCMC Uncertainty
Not stored in the current pipeline (single seed). To use uncertainties in training, block standard error should be added to the output dictionary.

### Train/Test Split
Leakage via near-duplicate frameworks should be avoided using MOFid similarity filtering before splitting. The current pipeline does not enforce this — it must be done at the dataset level.

---

## 17. Edge Cases and Failure Handling

| Failure Mode | Current Handling |
|---|---|
| RASPA crash (returncode ≠ 0) | Prints stderr, raises `RuntimeError` |
| Output parse failure | `try/except`, returns `{"info": error_msg, "results": None}` |
| All-zero NH₃ loading | Parsed as 0.0 mol/kg (valid result for non-adsorbing MOFs) |
| HVF calculation failure | Falls back to `helium_void_fraction = 1.0` with warning |
| Non-neutral cell | `verify_charge_neutrality()` checks |q_net| ≤ 0.05 e after charge assignment. Structures failing this are excluded and logged — not passed to RASPA. |
| Inaccessible pores | ASA > 0 pre-filter; RASPA handles internally (no insertion into inaccessible void) |
| Unphysical loading cap | Not implemented — no maximum loading cutoff is enforced |

---

## 18. Full Parameter Table (Reporting Checklist)

| Parameter | Value |
|---|---|
| **Simulation code** | RASPA2 |
| **MOF force field** | UFF (internal RASPA lookup) |
| **NH₃ model** | TraPPE-EH, 5-site rigid |
| **H₂O model** | TIP4P/2005, 4-site rigid |
| **N₂ model** | TraPPE, **3-site** (2 N + COM charge site), Potoff & Siepmann 2001 |
| **Combining rules** | Lorentz-Berthelot |
| **MOF partial charges** | PACMOF2 (DDEC6-quality ML) |
| **Adsorbate charges** | Fixed (TraPPE-EH / TIP4P/2005 values) |
| **VDW cutoff** | 12.0 Å |
| **Coulomb cutoff** | 12.0 Å |
| **LJ potential** | Shifted (not truncated) |
| **Tail corrections** | None |
| **Electrostatics** | Ewald, precision = 10⁻⁶ |
| **Supercell rule** | min perpendicular width ≥ 24 Å (2 × cutoff) |
| **Equilibration cycles** | 20,000 |
| **Production cycles** | 20,000 |
| **Helium void fraction** | Widom insertion, 10,000 cycles, 298 K |
| **Translation prob.** | 0.5 |
| **Rotation prob.** | 0.5 |
| **Reinsertion prob.** | 0.5 |
| **Swap prob.** | 2.0 |
| **Identity change prob. (mixture)** | 1.0 |
| **Framework flexibility** | Rigid |
| **Random seed** | 4 (fixed) |
| **Replicas per condition** | 1 |
| **Reporting units** | mol/kg (= mmol/g) |
| **Uptake type** | Absolute |
| **Pure NH₃ state point** | 298 K, 1.0 bar |
| **VSA ads. state point** | 298 K, 1.0 bar, NH₃:H₂O:N₂ = 0.001:0.03:0.969 |
| **VSA des. state point** | 400 K, 0.01 bar, pure NH₃ (**idealized** — H₂O/N₂ excluded from regeneration) |
| **Uncertainty reporting** | Convergence trend spread (std of block-window means); NOT a valid SE. Publication SE: ≥3 independent seeds, SE = std / √n. |


---

## 19. Known Limitations and Status

| # | Issue (Staff Review) | Status |
|---|---|---|
| 1 | **ACS Omega citation misapplied to humid-air VSA** | ✅ **Fixed** — Citation now scoped to TraPPE FF + RASPA practice only. VSA scenario noted as original contribution. |
| 2 | **"Well below saturation" water claim** | ✅ **Fixed** — Reframed: 3% H₂O ≈ 95% RH, near-saturation; clustering/slow equilibration explicitly acknowledged. |
| 3 | **Zeo++ probe/filter inconsistency** | ✅ **Fixed** — He-probe (1.2 Å) labelled as coarse pre-filter; NH₃ probe (1.81 Å) noted as recommended improvement. |
| 4 | **Block SE invalid on cumulative averages** | ✅ **Fixed** — `block_standard_error()` docstring corrected; output field relabelled as convergence trend spread, NOT SE. |
| 5 | **Charge neutrality unchecked** | ✅ **Fixed** — `verify_charge_neutrality()` added and wired into both simulation entry points; structures with |q| > 0.05 e excluded. |
| 6 | **Helium void fraction field unexplained** | ✅ **Fixed** — Section 12 now explains `Average Widom Rosenbluth-weight` = ⟨exp(−βu)⟩, which converges to geometric void fraction for He at 298 K. |
| 7 | **N₂ 2-site vs 3-site conflict** | ✅ **Fixed** — Parameter table updated to "TraPPE, 3-site (2 N + COM charge site)". |
| 8 | **Desorption not labelled as idealized** | ✅ **Fixed** — Parameter table and Section 13 now explicitly label desorption as "idealized state (H₂O/N₂ excluded)". |
| 9 | **Triclinic perpendicular width incorrect** | ✅ **Fixed** — `calculate_unit_cells()` now uses exact volume-based formula d = V/(face area). Section 8 updated. |
| 10 | **Rigid framework for NH₃** | ⬜ Known limitation — flexibility prohibitively expensive at screening scale. |
| 11 | **Single seed; no replication** | ⬜ Publication SE requires ≥3 seeds; current output provides convergence trend only. |
| 12 | **Cutoff sensitivity test** | ⬜ Out of scope — requires external runs; recommended for top candidates. |

---

## References

| Item | Citation |
|---|---|
| GCMC engine | Dubbeldam, D. et al. *RASPA: Molecular Simulation Software for Adsorption and Diffusion in Flexible Nanoporous Materials.* Mol. Simul. **42**, 81–101 (2016) |
| Protocol adaptation | *Computational Screening of Metal–Organic Frameworks for Ammonia Capture from H₂/N₂/NH₃ Mixtures.* ACS Omega **2022**. DOI: 10.1021/acsomega.2c04517 |
| NH₃ force field (TraPPE-EH) | Eckl, B. et al. *On the Application of Force Fields for Predicting a Wide Variety of Properties: Ethylene Oxide as an Example.* J. Phys. Chem. B **2008**. |
| H₂O model (TIP4P/2005) | Abascal, J.L.F. & Vega, C. *A general purpose model for the condensed phases of water: TIP4P/2005.* J. Chem. Phys. **2005**. |
| N₂ model (TraPPE) | Potoff, J.J. & Siepmann, J.I. *Vapor–liquid equilibria of mixtures containing alkanes, carbon dioxide, and nitrogen.* AIChE J. **2001**. |
| Framework charges | Raza, A. et al. *PACMOF2: Transferable Machine-Learning Approach for Predicting Partial Charges in Metal-Organic Frameworks.* npj Comput. Mater. **2022**. |
| UFF force field | Rappe, A.K. et al. *UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular Dynamics Simulations.* J. Am. Chem. Soc. **1992**. |
