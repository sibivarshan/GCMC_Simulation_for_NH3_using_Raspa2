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


def extract_raspa_output_nh3(raspa_output, components):
    """
    Parse RASPA output for arbitrary component loadings and enthalpies.

    Args:
        raspa_output (str): Raw RASPA .data output text.
        components (list): List of component names, e.g. ["Ammonia", "H2O", "N2"]

    Returns:
        dict: {component_name: {"loading_mol_kg": float, "enthalpy_kJ_mol": float}}
    """
    results = {}

    final_loading_section = re.findall(
        r"Number of molecules:\n=+[^=]*(?=)", raspa_output
    )[0]

    enthalpy_section = re.findall(
        r"Enthalpy of adsorption:\n={2,}\n(.+?)\n={2,}", raspa_output, re.DOTALL
    )[0]

    for component in components:
        # --- Loading ---
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

        # --- Enthalpy of adsorption ---
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
            "enthalpy_kJ_mol": enthalpy_kJ_mol,
            "heat_of_adsorption_kJ_mol": heat_of_adsorption,
        }

    return results


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
        dict with NH3 uptake in mol/kg and mmol/g, plus heat of adsorption.
    """
    from gcmc_wrapper_nh3 import GCMCSimulation, calculate_charges, run_gcmc_simulation

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
        dict with working capacity, selectivities, and individual uptakes.
    """
    from gcmc_wrapper_nh3 import GCMCSimulation, calculate_charges, run_gcmc_simulation

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
    ads_n2 = ads_results["N2"]["loading_mol_kg"]
    des_nh3 = des_results["Ammonia"]["loading_mol_kg"]

    # Selectivity = (x_NH3/y_NH3) / (x_other/y_other)
    # where x = adsorbed mole fraction, y = gas phase mole fraction
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
        # Selectivities
        "NH3_N2_selectivity": nh3_n2_selectivity,
        "NH3_H2O_selectivity": nh3_h2o_selectivity,
        # Adsorption (298K, 1 bar, humid air)
        "NH3_uptake_ads_mmol_g": ads_nh3,
        "H2O_uptake_ads_mmol_g": ads_h2o,
        "N2_uptake_ads_mmol_g": ads_n2,
        "NH3_Qst_ads_kJ_mol": ads_results["Ammonia"]["heat_of_adsorption_kJ_mol"],
        # Desorption (400K, 0.01 bar, pure NH3)
        "NH3_uptake_des_mmol_g": des_nh3,
        "NH3_Qst_des_kJ_mol": des_results["Ammonia"]["heat_of_adsorption_kJ_mol"],
    }

    return output
