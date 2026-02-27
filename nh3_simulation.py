"""
NH3 GCMC Simulation Module
==========================
Adapted from MOFDiff's mofdiff/gcmc/simulation.py for NH3 capture from humid air.

Reference paper: "Computational Screening of MOFs for Ammonia Capture from H2/N2/NH3 Mixtures"
ACS Omega 2022, DOI: 10.1021/acsomega.2c04517

Simulation scenarios:
1. Pure NH3 at 1 bar / 298K  (for validation against raw_nh3_core dataset)
2. 3-component air removal: NH3/H2O/N2 at 298K/1bar → desorption at 400K/0.01bar
"""

import random
import re
import numpy as np


# ---------------------------------------------------------------------------
# Block-averaging uncertainty estimation
# ---------------------------------------------------------------------------

def _parse_loading_history(raspa_output, component_name):
    """
    Extract the per-cycle loading (mol/kg) from the RASPA output text
    for a named component.

    RASPA prints lines like:
      Average loading absolute [mol/kg framework] <N0> +/- <dN0> [cycle <i>]
    or, in the per-block averages section, lines starting with the cycle index.
    We parse the 'running averages' block which RASPA prints at every PrintEvery
    interval, accumulating them into a time-series for block SE.

    Returns
    -------
    list of float  (one entry per print checkpoint, or empty list if not found)
    """
    # RASPA prints:
    # "Current cycle: <N>  out of <M>"
    # followed by component-specific lines including absolute loading.
    # We collect the instantaneous-average-so-far values at each checkpoint.
    # Pattern: "Average loading absolute [mol/kg framework] <value>"
    # These are cumulative averages printed by RASPA at each PrintEvery step.
    # They are NOT block-independent; we use them to construct block estimates.

    # More robustly: look for the per-cycle section.
    # RASPA actually only outputs running averages (not per-move raw data).
    # We therefore parse the history of cumulative averages and,
    # treating each PrintEvery window as an independent block, compute SE.

    # Regex: capture all occurrences of the component's absolute loading printout
    pattern = (
        rf"Component \d \[{re.escape(component_name)}\]"
        rf".*?Average loading absolute \[mol/kg framework\]\s+([\d.]+)"
    )
    matches = re.findall(pattern, raspa_output, re.DOTALL)
    return [float(m) for m in matches]


def block_standard_error(values, n_blocks=5):
    """
    Approximate convergence indicator from RASPA cumulative-average checkpoints.

    IMPORTANT LIMITATION: RASPA's per-PrintEvery output contains *cumulative*
    running averages, not independent instantaneous samples. Applying the
    Flyvbjerg-Petersen block method to cumulative averages does NOT yield a
    valid statistical standard error — consecutive cumulative averages are
    strongly correlated by construction.

    This function is therefore provided only as a convergence trend indicator:
    - A large spread across block means suggests the simulation has not converged.
    - A near-zero spread suggests the cumulative average has stabilised.

    For a publication-quality SE, run >=3 independent seeds and compute
    SE = std(seed_means) / sqrt(n_seeds).

    Parameters
    ----------
    values   : list/array of cumulative-average loading values from RASPA checkpoints
    n_blocks : int, number of windows (default 5)

    Returns
    -------
    (mean, trend_spread)  where trend_spread is the std of block-window means
                          (NOT a statistically valid SE). Returns (mean, None)
                          if too few checkpoints.
    """
    arr = np.array(values, dtype=float)
    n = len(arr)

    if n < n_blocks:
        mean_val = float(np.mean(arr)) if n > 0 else 0.0
        return mean_val, None

    trimmed = arr[: (n // n_blocks) * n_blocks]
    blocks = trimmed.reshape(n_blocks, -1)
    block_means = blocks.mean(axis=1)

    mean_val = float(block_means.mean())
    # This is std of block means, NOT a valid SE for cumulative averages
    trend_spread = float(block_means.std(ddof=1))
    return mean_val, trend_spread


# ---------------------------------------------------------------------------
# RASPA output parser
# ---------------------------------------------------------------------------

def extract_raspa_output_nh3(raspa_output, components):
    """
    Parse RASPA output for arbitrary component loadings, enthalpies,
    and block-averaged standard errors.

    Args:
        raspa_output (str): Raw RASPA .data output text.
        components (list): List of component names, e.g. ["Ammonia", "H2O", "N2"]

    Returns:
        dict: {component_name: {
                  "loading_mol_kg": float,       # mean absolute loading
                  "loading_se_mol_kg": float|None, # block SE (None if unavailable)
                  "enthalpy_kJ_mol": float,
                  "heat_of_adsorption_kJ_mol": float
              }}
    """
    results = {}

    final_loading_section = re.findall(
        r"Number of molecules:\n=+[^=]*(?=)", raspa_output
    )[0]

    enthalpy_section = re.findall(
        r"Enthalpy of adsorption:\n={2,}\n(.+?)\n={2,}", raspa_output, re.DOTALL
    )[0]

    for component in components:
        # ------------------------------------------------------------------
        # 1. Mean loading from the final summary block
        # ------------------------------------------------------------------
        comp_pattern = rf"Component \d \[{re.escape(component)}\].*?(?=Component|\Z)"
        comp_subsection = re.findall(comp_pattern, final_loading_section, re.DOTALL)

        if comp_subsection:
            loading = float(
                re.findall(
                    r"(?<=Average loading absolute \[mol/kg framework\])\s*\d*\.\d*",
                    comp_subsection[0],
                )[0]
            )
        else:
            loading = 0.0

        # ------------------------------------------------------------------
        # 2. Block standard error from the checkpoint history
        # ------------------------------------------------------------------
        history = _parse_loading_history(raspa_output, component)
        if len(history) >= 2:
            # Use only production-phase checkpoints:
            # RASPA prints checkpoints during both equilibration and production.
            # We take the second half of the history as the production window.
            prod_history = history[len(history) // 2:]
            _, se = block_standard_error(prod_history, n_blocks=min(5, len(prod_history)))
        else:
            se = None

        # ------------------------------------------------------------------
        # 3. Enthalpy of adsorption
        # ------------------------------------------------------------------
        if len(components) > 1:
            # Multi-component: look for per-component enthalpy
            comp_enthalpy_pattern = rf"\[{re.escape(component)}\].*?(?=component|\Z)"
            comp_enthalpy = re.findall(comp_enthalpy_pattern, enthalpy_section, re.DOTALL)
            if comp_enthalpy:
                enthalpy_K = float(
                    re.findall(r"(?<=\[K\])\s*-?\d*\.\d*", comp_enthalpy[0])[0]
                )
            else:
                enthalpy_K = 0.0
        else:
            # Single component: use total enthalpy
            total_enthalpy_match = re.findall(
                r"Total enthalpy of adsorption.*?(?=Q=-H|\Z)",
                enthalpy_section,
                re.DOTALL,
            )
            if total_enthalpy_match:
                enthalpy_K = float(
                    re.findall(r"(?<=\[K\])\s*-?\d*\.\d*", total_enthalpy_match[0])[0]
                )
            else:
                enthalpy_K = 0.0

        # Convert enthalpy from K to kJ/mol: multiply by R (8.314 J/(mol·K)) / 1000
        enthalpy_kJ_mol = enthalpy_K * 8.314 / 1000.0
        heat_of_adsorption = -1 * enthalpy_kJ_mol  # Qst = -H

        results[component] = {
            "loading_mol_kg": loading,
            "loading_se_mol_kg": se,          # block SE; None if < 2 checkpoints
            "enthalpy_kJ_mol": enthalpy_kJ_mol,
            "heat_of_adsorption_kJ_mol": heat_of_adsorption,
        }

    return results


# ---------------------------------------------------------------------------
# Simulation entry points
# ---------------------------------------------------------------------------

def nh3_uptake_pure(
    cif_file,
    calc_charges=True,
    charge_method="pacmof2",
    temperature=298,
    pressure=100000,  # 1 bar in Pa
    rundir="./temp",
    rewrite_raspa_input=False,
):
    """
    Pure NH3 adsorption at specified T and P.
    Default: 1 bar / 298K to match the raw_nh3_core validation dataset.

    Returns:
        dict with NH3 uptake in mol/kg and mmol/g, block SE, plus heat of adsorption.
    """
    from gcmc_wrapper_nh3 import (
        GCMCSimulation, calculate_charges, run_gcmc_simulation,
        verify_charge_neutrality,
    )

    random.seed(4)
    np.random.seed(4)

    sim = GCMCSimulation(
        cif_file,
        sorbates=["Ammonia"],
        sorbates_mol_fraction=[1.0],
        temperature=temperature,
        pressure=pressure,
        rundir=rundir,
    )

    if calc_charges:
        calculate_charges(sim, method=charge_method)
        _, neutral = verify_charge_neutrality(sim.sorbent_file)
        if neutral is False:
            raise RuntimeError(
                f"Charge neutrality check FAILED for {cif_file}: "
                f"net charge exceeds 0.05 e. Exclude this structure."
            )

    run_gcmc_simulation(
        sim,
        equilibration_cycles=20000,
        production_cycles=20000,
        rewrite_raspa_input=rewrite_raspa_input,
    )

    results = extract_raspa_output_nh3(sim.raspa_output, components=["Ammonia"])

    nh3 = results["Ammonia"]
    output = {
        "file": str(cif_file),
        "mode": "pure_nh3",
        "temperature_K": temperature,
        "pressure_Pa": pressure,
        "NH3_uptake_mol_kg": nh3["loading_mol_kg"],
        "NH3_uptake_mmol_g": nh3["loading_mol_kg"],  # mol/kg == mmol/g
        "NH3_uptake_se_mmol_g": nh3["loading_se_mol_kg"],   # block SE; None if unavailable
        "NH3_heat_of_adsorption_kJ_mol": nh3["heat_of_adsorption_kJ_mol"],
    }

    return output


def working_capacity_nh3_air_removal(
    cif_file,
    calc_charges=True,
    charge_method="pacmof2",
    rundir="./temp",
    rewrite_raspa_input=False,
    # Adsorption conditions: humid air with trace NH3
    ads_temperature=298,
    ads_pressure=100000,  # 1 bar
    ads_nh3_fraction=0.001,  # 1000 ppm NH3
    ads_h2o_fraction=0.03,  # 3% relative humidity
    ads_n2_fraction=0.969,  # balance N2
    # Desorption conditions: vacuum regeneration
    des_temperature=400,  # K
    des_pressure=1000,  # 0.01 bar
):
    """
    NH3 removal from humid air via Vacuum Swing Adsorption (VSA).

    Adsorption:  298 K, 1 bar, NH3:H2O:N2 = 0.001:0.03:0.969
    Desorption:  400 K, 0.01 bar, pure NH3
    Working capacity = NH3_uptake(ads) - NH3_uptake(des)

    Returns:
        dict with working capacity, selectivities, individual uptakes, and block SEs.
    """
    from gcmc_wrapper_nh3 import (
        GCMCSimulation, calculate_charges, run_gcmc_simulation,
        verify_charge_neutrality,
    )

    random.seed(4)
    np.random.seed(4)

    # ========================
    # STEP 1: Adsorption GCMC
    # ========================
    ads_sim = GCMCSimulation(
        cif_file,
        sorbates=["Ammonia", "H2O", "N2"],
        sorbates_mol_fraction=[ads_nh3_fraction, ads_h2o_fraction, ads_n2_fraction],
        temperature=ads_temperature,
        pressure=ads_pressure,
        rundir=rundir,
    )

    if calc_charges:
        calculate_charges(ads_sim, method=charge_method)
        _, neutral = verify_charge_neutrality(ads_sim.sorbent_file)
        if neutral is False:
            raise RuntimeError(
                f"Charge neutrality check FAILED for {cif_file}: "
                f"net charge exceeds 0.05 e. Exclude this structure."
            )

    run_gcmc_simulation(
        ads_sim,
        equilibration_cycles=20000,
        production_cycles=20000,
        rewrite_raspa_input=rewrite_raspa_input,
    )

    ads_results = extract_raspa_output_nh3(
        ads_sim.raspa_output, components=["Ammonia", "H2O", "N2"]
    )

    # ========================
    # STEP 2: Desorption GCMC
    # ========================
    des_sim = GCMCSimulation(
        cif_file,
        sorbates=["Ammonia"],
        sorbates_mol_fraction=[1.0],
        temperature=des_temperature,
        pressure=des_pressure,
        rundir=rundir,
    )

    if calc_charges:
        calculate_charges(des_sim, method=charge_method)

    run_gcmc_simulation(
        des_sim,
        equilibration_cycles=20000,
        production_cycles=20000,
        rewrite_raspa_input=rewrite_raspa_input,
    )

    des_results = extract_raspa_output_nh3(
        des_sim.raspa_output, components=["Ammonia"]
    )

    # ========================
    # STEP 3: Compute metrics
    # ========================
    ads_nh3 = ads_results["Ammonia"]["loading_mol_kg"]
    ads_h2o = ads_results["H2O"]["loading_mol_kg"]
    ads_n2  = ads_results["N2"]["loading_mol_kg"]
    des_nh3 = des_results["Ammonia"]["loading_mol_kg"]

    ads_nh3_se = ads_results["Ammonia"]["loading_se_mol_kg"]
    des_nh3_se = des_results["Ammonia"]["loading_se_mol_kg"]

    # Propagated working-capacity SE (independent runs → sum in quadrature)
    if ads_nh3_se is not None and des_nh3_se is not None:
        wc_se = float(np.sqrt(ads_nh3_se**2 + des_nh3_se**2))
    else:
        wc_se = None

    # Selectivity = (x_NH3/y_NH3) / (x_other/y_other)
    # where x = adsorbed loading (mol/kg), y = gas phase mole fraction
    nh3_n2_selectivity = None
    nh3_h2o_selectivity = None

    if ads_n2 > 0 and ads_n2_fraction > 0:
        nh3_n2_selectivity = (ads_nh3 / ads_nh3_fraction) / (
            ads_n2 / ads_n2_fraction
        )

    if ads_h2o > 0 and ads_h2o_fraction > 0:
        nh3_h2o_selectivity = (ads_nh3 / ads_nh3_fraction) / (
            ads_h2o / ads_h2o_fraction
        )

    working_capacity = ads_nh3 - des_nh3

    output = {
        "file": str(cif_file),
        "mode": "air_removal_vsa",
        # Working capacity
        "working_capacity_nh3_mmol_g": working_capacity,
        "working_capacity_se_mmol_g": wc_se,          # propagated block SE
        # Selectivities
        "NH3_N2_selectivity": nh3_n2_selectivity,
        "NH3_H2O_selectivity": nh3_h2o_selectivity,
        # Adsorption (298K, 1 bar, humid air)
        "NH3_uptake_ads_mmol_g": ads_nh3,
        "NH3_uptake_ads_se_mmol_g": ads_nh3_se,
        "H2O_uptake_ads_mmol_g": ads_h2o,
        "N2_uptake_ads_mmol_g": ads_n2,
        "NH3_Qst_ads_kJ_mol": ads_results["Ammonia"]["heat_of_adsorption_kJ_mol"],
        # Desorption (400K, 0.01 bar, pure NH3)
        "NH3_uptake_des_mmol_g": des_nh3,
        "NH3_uptake_des_se_mmol_g": des_nh3_se,
        "NH3_Qst_des_kJ_mol": des_results["Ammonia"]["heat_of_adsorption_kJ_mol"],
    }

    return output
