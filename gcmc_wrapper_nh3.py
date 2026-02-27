"""
GCMC Wrapper for NH3 Simulations
=================================
Adapted from MOFDiff's mofdiff/gcmc/gcmc_wrapper.py for NH3 capture.

Key changes from the original CO2 wrapper:
- Helium void fraction calculated via Widom insertion (not assumed 1.0)
- Higher equilibration/production cycles (20,000 each - paper validated)
- Configurable SwapProbability (higher for polar molecules)
- Support for PACMOF2 charges (in addition to MEPO-Qeq)
- Molecule definitions from local directory (not RASPA installation)
"""

import os
import subprocess
import re
import shutil
from pathlib import Path
from time import time
from math import cos, sin, radians
from textwrap import dedent

try:
    from openbabel.pybel import readfile
except ImportError:
    from openbabel import pybel
    readfile = pybel.readfile

# ---------------------
# Environment Variables
# ---------------------
# Load from .env if available
try:
    from dotenv import load_dotenv
    load_dotenv(Path(__file__).parent / ".env")
except ImportError:
    pass

RASPA_PATH = os.getenv("RASPA_PATH", "")
RASPA_SIM_PATH = os.getenv("RASPA_SIM_PATH", "simulate")
EGULP_PATH = os.getenv("EGULP_PATH", "")
EGULP_PARAMETER_PATH = os.getenv("EGULP_PARAMETER_PATH", "")
GCMC_DIR = Path(__file__).parent  # directory where this script lives


class GCMCSimulation:
    """Parent class for aggregating GCMC simulation parameters and results."""

    def __init__(
        self,
        cif_file,
        sorbates=None,
        sorbates_mol_fraction=None,
        temperature=298,
        pressure=101325,
        rundir="./temp",
    ):
        if sorbates is None:
            sorbates = ["Ammonia"]
        if sorbates_mol_fraction is None:
            sorbates_mol_fraction = [1.0]

        # Pull unit cell dimensions from cif_file
        self.sorbent = next(readfile("cif", str(cif_file)))
        self.dim = [0, 0, 0]
        self.angle = [0, 0, 0]

        with open(cif_file, "r") as file:
            content = file.read()
        dim_match = re.findall(r"_cell_length_.\s+\d+\.?\d*", content)
        angle_match = re.findall(r"_cell_angle_\S+\s+\d+\.?\d*", content)

        for i in range(3):
            self.dim[i] = float(re.findall(r"\d+\.?\d*", dim_match[i])[0])
            self.angle[i] = float(re.findall(r"\d+\.?\d*", angle_match[i])[0])

        # Unique identifier for this simulation
        basename = ".".join(str(cif_file).split("/")[-1].split("\\")[-1].split(".")[:-1])
        self.identifier = basename + "_" + str(time()).split(".")[1]

        # Run directory
        self.rundir = Path(rundir)
        self.rundir.mkdir(parents=True, exist_ok=True)
        self.sorbent_file = str(cif_file)

        # Simulation conditions
        self.sorbates = sorbates
        self.sorbates_mol_fraction = [
            i / sum(sorbates_mol_fraction) for i in sorbates_mol_fraction
        ]
        self.temperature = temperature
        self.pressure = pressure

        # Calculated / filled in later
        self.helium_void_fraction = None  # Will be computed via Widom insertion
        self.rosenbluth_weights = [1.0 for _ in range(len(sorbates))]
        self.raspa_config = None
        self.raspa_output = None

    def get_sorbate_radius(self, sorbate):
        """Returns kinetic radius for common sorbates (Angstrom)."""
        kinetic_diameter = {
            "He": 2.551, "Ne": 2.82, "Ar": 3.542, "Kr": 3.655, "Xe": 4.047,
            "H2": 2.8585, "N2": 3.72, "O2": 3.467,
            "CO": 3.69, "CO2": 3.3, "NO": 3.492, "N2O": 3.838, "SO2": 4.112,
            "H2O": 2.641, "CH4": 3.758, "NH3": 3.62, "Ammonia": 3.62,
            "H2S": 3.623,
        }
        try:
            return kinetic_diameter[sorbate] * 0.5
        except KeyError:
            print(f"Unknown sorbate: {sorbate}")
            raise

    def calculate_unit_cells(self, forcefield_cutoff):
        """
        Determine unit cell replications to satisfy the minimum-image convention.

        For a general triclinic cell with lattice parameters (a, b, c, α, β, γ),
        the perpendicular distances between opposite faces are:

            d_a = V / (b * c * sin(α))
            d_b = V / (a * c * sin(β))
            d_c = V / (a * b * sin(γ))

        where V = a*b*c * sqrt(1 - cos²α - cos²β - cos²γ + 2·cosα·cosβ·cosγ)

        We require: unit_cells[i] * d_i >= 2 * cutoff  for all i.

        This replaces the previous approximate formula (dim[i] * |cos(angle[i] - 90°)|)
        which is only correct for cells where each axis is orthogonal to the other two,
        and can undercount replications in highly oblique cells.
        """
        from math import sqrt as _sqrt

        a, b, c = self.dim
        # Angles stored as [alpha, beta, gamma] (from CIF _cell_angle_alpha/beta/gamma)
        alpha, beta, gamma = [radians(ang) for ang in self.angle]

        cos_a, cos_b, cos_g = cos(alpha), cos(beta), cos(gamma)
        sin_a, sin_b, sin_g = abs(sin(alpha)), abs(sin(beta)), abs(sin(gamma))

        # Cell volume divided by abc  → normalised volume factor
        det = 1.0 - cos_a**2 - cos_b**2 - cos_g**2 + 2.0*cos_a*cos_b*cos_g
        # Guard against numerical noise (det should be > 0 for a valid cell)
        if det <= 0:
            print("  Warning: degenerate cell angles; falling back to orthogonal approximation")
            perp = [a * sin_a, b * sin_b, c * sin_g]
        else:
            vol_factor = _sqrt(det)         # V / (a*b*c)
            # Perpendicular widths via V / (face area)
            # face bc area (per abc) = sin(alpha); similarly for others
            d_a = a * vol_factor / (sin_b * sin_g) if (sin_b * sin_g) > 1e-10 else a
            d_b = b * vol_factor / (sin_a * sin_g) if (sin_a * sin_g) > 1e-10 else b
            d_c = c * vol_factor / (sin_a * sin_b) if (sin_a * sin_b) > 1e-10 else c
            perp = [d_a, d_b, d_c]

        unit_cells = [1, 1, 1]
        for i in range(3):
            if perp[i] <= 0:
                perp[i] = self.dim[i]   # safe fallback
            while unit_cells[i] * perp[i] < 2 * forcefield_cutoff:
                unit_cells[i] += 1

        return unit_cells

    def write_out(self, output_path):
        """Write RASPA output to a file."""
        with open(output_path, "w") as log_file:
            log_file.write(self.raspa_output)


def compute_helium_void_fraction(simulation, forcefield_cutoff=12):
    """
    Calculate helium void fraction via Widom insertion in RASPA.
    This is critical for accurate uptake calculations.
    The ACS Omega paper uses this approach.
    """
    # Determine force field and molecule definition paths
    ff_dir = GCMC_DIR / "force_field"
    mol_dir = GCMC_DIR / "molecule_definitions"

    workdir = simulation.rundir / "hvf" / simulation.identifier
    workdir.mkdir(exist_ok=True, parents=True)

    sorbent_file = ".".join(
        simulation.sorbent_file.replace("\\", "/").split("/")[-1].split(".")[:-1]
    )

    # Copy CIF to RASPA structures directory if RASPA_PATH is set
    if RASPA_PATH:
        struct_dir = Path(RASPA_PATH) / "share" / "raspa" / "structures" / "cif"
        struct_dir.mkdir(parents=True, exist_ok=True)
        try:
            shutil.copy(simulation.sorbent_file, str(struct_dir))
        except Exception:
            pass  # Fallback: use workdir copy below

    # Always copy CIF to working directory so RASPA can find it
    shutil.copy(simulation.sorbent_file, str(workdir))

    unit_cells = simulation.calculate_unit_cells(forcefield_cutoff)

    hvf_config = dedent(f"""
        SimulationType                Widom
        NumberOfCycles                10000
        PrintEvery                    1000

        Forcefield                    Local
        UseChargesFromCIFFile         yes
        CutOffVDW                     {forcefield_cutoff}
        CutOffChargeCharge            {forcefield_cutoff}
        ChargeMethod                  Ewald
        EwaldPrecision                1e-6

        Framework                     0
        FrameworkName                 {sorbent_file}
        InputFileType                 cif
        UnitCells                     {unit_cells[0]} {unit_cells[1]} {unit_cells[2]}
        ExternalTemperature           298

        Component 0 MoleculeName            helium
                 MoleculeDefinition            TraPPE
                 WidomProbability              1.0
                 CreateNumberOfMolecules       0
    """).strip()

    # Write input file
    with open(workdir / "simulation.input", "w") as f:
        f.write(hvf_config)

    # Copy force field files
    if ff_dir.exists():
        for ff_file in ff_dir.glob("*"):
            shutil.copy(str(ff_file), str(workdir))

    # Run RASPA
    subprocess.run(
        [RASPA_SIM_PATH, "simulation.input"],
        cwd=str(workdir),
        capture_output=True,
    )

    # Parse output for HVF
    try:
        output_dir = workdir / "Output" / "System_0"
        file_list = os.listdir(str(output_dir))
        raspa_log = [item for item in file_list if re.match(r".*\.data", item)][0]

        with open(str(output_dir / raspa_log), "r") as log:
            output_text = log.read()

        hvf_match = re.findall(
            r"Average Widom Rosenbluth-weight:\s*\d*\.\d*", output_text
        )
        if hvf_match:
            hvf = float(re.findall(r"\d+\.\d+", hvf_match[0])[0])
            simulation.helium_void_fraction = hvf
            print(f"  Helium void fraction: {hvf:.4f}")
        else:
            print("  Warning: Could not parse HVF from RASPA output, using 1.0")
            simulation.helium_void_fraction = 1.0
    except Exception as e:
        print(f"  Warning: HVF calculation failed ({e}), using 1.0")
        simulation.helium_void_fraction = 1.0

def verify_charge_neutrality(cif_file, tolerance=0.05):
    """
    Verify that the unit cell net charge is within `tolerance` elementary charges
    after partial charges have been written into the CIF.

    This must be called explicitly after calculate_charges() to fulfil the
    methodological requirement of reporting charge neutrality in the paper.

    Parameters
    ----------
    cif_file : str or Path
        Path to the CIF file that contains _atom_site_charge (or _atom_site_partial_charge)
        columns written by PACMOF2 or MEPO-Qeq.
    tolerance : float
        Maximum absolute net charge (in e) to pass. Structures outside this range
        are flagged for exclusion. Default 0.05 e is conservative but publication-standard.

    Returns
    -------
    (net_charge, passed) : (float, bool)
        net_charge : summed charge of all atoms in one unit cell
        passed     : True if |net_charge| <= tolerance
    """
    import re
    with open(str(cif_file), "r") as f:
        content = f.read()

    # Try both common CIF charge column names
    # Match numeric values (positive or negative, integer or float)
    charge_pattern = re.compile(
        r"_atom_site(?:_partial)?_charge", re.IGNORECASE
    )
    if not charge_pattern.search(content):
        print(f"  Warning: no charge column found in {cif_file}; cannot verify neutrality.")
        return None, None

    # Extract all numeric values that follow the charge column in the atom loop.
    # Strategy: find the loop_ block containing the charge column and parse the
    # column index.
    lines = content.split("\n")
    in_loop = False
    headers = []
    charge_col_idx = None
    net_charge = 0.0
    n_atoms = 0

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.lower() == "loop_":
            in_loop = True
            headers = []
            charge_col_idx = None
            i += 1
            continue
        if in_loop:
            if line.startswith("_"):
                headers.append(line.lower())
                if "charge" in line.lower():
                    charge_col_idx = len(headers) - 1
            elif line == "" or line.startswith("loop_") or line.startswith("#"):
                in_loop = (line.lower() == "loop_")
                if line.lower() == "loop_":
                    headers = []
                    charge_col_idx = None
            else:
                if charge_col_idx is not None and headers:
                    # This is a data row
                    tokens = line.split()
                    if len(tokens) >= len(headers):
                        try:
                            q = float(tokens[charge_col_idx])
                            net_charge += q
                            n_atoms += 1
                        except (ValueError, IndexError):
                            pass
        i += 1

    if n_atoms == 0:
        print(f"  Warning: no atom charge data parsed from {cif_file}.")
        return None, None

    passed = abs(net_charge) <= tolerance
    status = "PASS" if passed else "FAIL (excluded)"
    print(f"  Charge neutrality check: net = {net_charge:+.4f} e over {n_atoms} atoms — {status}")
    return net_charge, passed


def calculate_charges(simulation, method="pacmof2"):
    """
    Calculate partial charges for the MOF framework.

    After charges are written, callers should invoke verify_charge_neutrality()
    on the resulting CIF to confirm |net charge| < 0.05 e before simulation.

    Args:
        simulation: GCMCSimulation object
        method: "pacmof2" (ML-predicted DDEC6) or "mepo_qeq" (eGULP)
    """
    if method == "pacmof2":
        _calculate_pacmof2_charges(simulation)
    elif method == "mepo_qeq":
        _calculate_mepo_qeq_charges(simulation)
    else:
        raise ValueError(f"Unknown charge method: {method}. Use 'pacmof2' or 'mepo_qeq'.")


def _calculate_pacmof2_charges(simulation):
    """
    Use PACMOF2 to predict DDEC6 charges for the MOF.
    PACMOF2 is an ML model that approximates DFT-quality DDEC6 charges.
    """
    try:
        from pacmof2.pacmof2 import get_charges
    except ImportError:
        raise ImportError(
            "PACMOF2 is required. Install from: https://github.com/snurr-group/PACMOF2"
        )

    print(f"  Calculating PACMOF2 (DDEC6) charges for {simulation.identifier}...")

    rundir = simulation.rundir / "charges" / simulation.identifier
    rundir.mkdir(exist_ok=True, parents=True)

    # PACMOF2 API: get_charges(path_to_cif, output_dir)
    # Writes {basename}_pacmof.cif to output_dir
    get_charges(simulation.sorbent_file, str(rundir))

    # The output file is named {basename}_pacmof.cif
    import os
    basename = os.path.splitext(os.path.basename(simulation.sorbent_file))[0]
    charged_cif = str(rundir / f"{basename}_pacmof.cif")

    # Update sorbent file path to the charged CIF
    simulation.sorbent_file = charged_cif
    print(f"  PACMOF2 charges written to {charged_cif}")


def _calculate_mepo_qeq_charges(simulation, egulp_parameter_set="MEPO"):
    """
    Calculate charges using MEPO-Qeq method via eGULP (fallback method).
    Adapted from MOFDiff's original implementation.
    """
    if not EGULP_PATH:
        raise RuntimeError(
            "EGULP_PATH not set. Configure in .env or set environment variable."
        )

    print(f"  Calculating MEPO-Qeq charges for {simulation.identifier}...")

    simulation.sorbent_file = str(simulation.rundir / f"{simulation.identifier}.cif")
    simulation.sorbent.write("cif", simulation.sorbent_file)

    rundir = simulation.rundir / "charges" / simulation.identifier
    rundir.mkdir(exist_ok=True, parents=True)

    config = dedent("""
        build_grid 0
        build_grid_from_scratch 1 none 0.25 0.25 0.25 1.0 2.0 0 0.3
        save_grid 0 grid.cube
        calculate_pot_diff 0
        calculate_pot 0 repeat.cube
        skip_everything 0
        point_charges_present 0
        include_pceq 0
        imethod 0
    """).strip()

    with open(rundir / "temp_config.input", "w") as file:
        file.write(config)

    subprocess.run(
        [
            EGULP_PATH,
            simulation.sorbent_file,
            os.path.join(EGULP_PARAMETER_PATH, egulp_parameter_set + ".param"),
            "temp_config.input",
        ],
        cwd=str(rundir),
    )

    simulation.sorbent_file = str(rundir / "charges.cif")
    print(f"  MEPO-Qeq charges written to {simulation.sorbent_file}")


def run_gcmc_simulation(
    simulation,
    initialization_cycles=0,
    equilibration_cycles=20000,
    production_cycles=20000,
    forcefield="Local",
    forcefield_cutoff=12,
    molecule_definitions="TraPPE",
    unit_cells=None,
    swap_probability=2.0,
    compute_hvf=True,
    cleanup=False,
    rewrite_raspa_input=False,
):
    """
    Run a GCMC simulation using RASPA2.

    Key differences from MOFDiff's original:
    - forcefield="Local" uses the force field files in the GCMC working directory
    - 20,000 eq/prod cycles (paper-validated)
    - swap_probability=2.0 for polar molecule insertion
    - Widom insertion for helium void fraction
    """
    if unit_cells is None:
        unit_cells = [0, 0, 0]

    # Compute helium void fraction if not already done
    if compute_hvf and simulation.helium_void_fraction is None:
        compute_helium_void_fraction(simulation, forcefield_cutoff)

    if simulation.helium_void_fraction is None:
        simulation.helium_void_fraction = 1.0

    workdir = simulation.rundir / "raspa_output" / simulation.identifier
    workdir.mkdir(exist_ok=True, parents=True)

    sorbent_file = ".".join(
        simulation.sorbent_file.replace("\\", "/").split("/")[-1].split(".")[:-1]
    )

    # Copy CIF to RASPA structures directory if RASPA_PATH is set
    if RASPA_PATH:
        struct_dir = Path(RASPA_PATH) / "share" / "raspa" / "structures" / "cif"
        struct_dir.mkdir(parents=True, exist_ok=True)
        try:
            shutil.copy(simulation.sorbent_file, str(struct_dir))
        except Exception:
            pass  # Fallback: use workdir copy below

    # Always copy CIF to working directory so RASPA can find it
    shutil.copy(simulation.sorbent_file, str(workdir))

    # Calculate number of unit cells if not user-defined
    if sum(unit_cells) == 0:
        unit_cells = simulation.calculate_unit_cells(forcefield_cutoff)

    # Copy force field files to working directory
    ff_dir = GCMC_DIR / "force_field"
    if ff_dir.exists():
        for ff_file in ff_dir.glob("*"):
            shutil.copy(str(ff_file), str(workdir))

    # Copy molecule definition files to working directory
    mol_dir = GCMC_DIR / "molecule_definitions"
    if mol_dir.exists():
        mol_dest = workdir / "molecules" / molecule_definitions
        mol_dest.mkdir(parents=True, exist_ok=True)
        src_dir = mol_dir / molecule_definitions
        if src_dir.exists():
            for mol_file in src_dir.glob("*.def"):
                shutil.copy(str(mol_file), str(mol_dest))

    # Build RASPA config
    simulation.raspa_config = dedent(f"""
        SimulationType                MonteCarlo
        NumberOfCycles                {production_cycles}
        NumberOfInitializationCycles  {initialization_cycles}
        NumberOfEquilibrationCycles   {equilibration_cycles}
        PrintEvery                    1000

        Forcefield                    {forcefield}
        UseChargesFromCIFFile         yes
        CutOffVDW                     {forcefield_cutoff}
        CutOffChargeCharge            {forcefield_cutoff}
        CutOffChargeBondDipole        {forcefield_cutoff}
        CutOffBondDipoleBondDipole    {forcefield_cutoff}
        ChargeMethod                  Ewald
        EwaldPrecision                1e-6

        Framework                     0
        FrameworkName                 {sorbent_file}
        InputFileType                 cif
        UnitCells                     {unit_cells[0]} {unit_cells[1]} {unit_cells[2]}
        HeliumVoidFraction            {simulation.helium_void_fraction}
        ExternalTemperature           {simulation.temperature}
        ExternalPressure              {simulation.pressure}

        Movies                        no
    """).strip()

    # Add sorbate components
    total_sorbates = len(simulation.sorbates)
    identity_change_prob = 1.0 if total_sorbates > 1 else 0.0
    sorbate_list = " ".join(str(n) for n in range(total_sorbates))

    for i in range(total_sorbates):
        sorbate = simulation.sorbates[i]
        sorbate_mol_fraction = simulation.sorbates_mol_fraction[i]
        rosenbluth_weight = simulation.rosenbluth_weights[i]

        simulation.raspa_config += "\n\n"
        simulation.raspa_config += dedent(f"""
            Component {i} MoleculeName            {sorbate}
                     MoleculeDefinition            {molecule_definitions}
                     MolFraction                   {sorbate_mol_fraction}
                     IdealGasRosenbluthWeight      {rosenbluth_weight}
                     IdentityChangeProbability     {identity_change_prob}
                       NumberOfIdentityChanges       {total_sorbates}
                       IdentityChangeList            {sorbate_list}
                     TranslationProbability        0.5
                     RotationProbability           0.5
                     ReinsertionProbability        0.5
                     SwapProbability               {swap_probability}
                     CreateNumberOfMolecules       0
        """).strip()

    # Write RASPA input file
    with open(workdir / "simulation.input", "w") as raspa_input:
        raspa_input.write(simulation.raspa_config)

    # Optionally rewrite input file to avoid RASPA parsing issues
    if rewrite_raspa_input:
        raspa_input_path = workdir / "simulation.input"
        # Read and rewrite to ensure proper line endings
        with open(raspa_input_path, "r") as f:
            content = f.read()
        with open(raspa_input_path, "w", newline="\n") as f:
            f.write(content)

    print(f"  Running RASPA GCMC: {sorbent_file}")
    print(f"    Sorbates: {simulation.sorbates}")
    print(f"    T={simulation.temperature}K, P={simulation.pressure}Pa")
    print(f"    Cycles: {equilibration_cycles} eq + {production_cycles} prod")
    print(f"    Unit cells: {unit_cells}")

    # Run RASPA simulation
    result = subprocess.run(
        [RASPA_SIM_PATH, "simulation.input"],
        cwd=str(workdir),
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(f"  RASPA stderr: {result.stderr[:500]}")

    # Collect RASPA output
    try:
        output_dir = workdir / "Output" / "System_0"
        file_list = os.listdir(str(output_dir))
        raspa_log = [item for item in file_list if re.match(r".*\.data", item)][0]

        with open(str(output_dir / raspa_log), "r") as log:
            simulation.raspa_output = log.read()
    except Exception as e:
        raise RuntimeError(
            f"Failed to read RASPA output for {sorbent_file}: {e}\n"
            f"RASPA stdout: {result.stdout[:500]}\n"
            f"RASPA stderr: {result.stderr[:500]}"
        )

    # Cleanup temp directory
    if cleanup:
        shutil.rmtree(str(workdir))

    return simulation
