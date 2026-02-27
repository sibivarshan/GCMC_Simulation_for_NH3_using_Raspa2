#!/bin/bash
###############################################################################
# NH3 GCMC Simulation — GCP Setup & Run Script
# ==============================================
# This script sets up the environment and runs GCMC simulations for NH3
# capture on MOFs using RASPA2 + PACMOF2.
#
# USAGE:
#   chmod +x setup_and_run.sh
#   ./setup_and_run.sh          # Full setup + test run of 1 MOF
#   ./setup_and_run.sh --skip-install --max_mofs 50   # Skip install, run 50
#
# REQUIREMENTS: Ubuntu (GCP VM), Python 3.8+, internet for pip
###############################################################################
set -euo pipefail

# ─── Configurable Parameters ────────────────────────────────────────────────
GCMC_DIR="$(cd "$(dirname "$0")" && pwd)"
MAX_MOFS="${MAX_MOFS:-1}"          # Override: MAX_MOFS=50 ./setup_and_run.sh
MODE="${MODE:-pure}"               # "pure" or "air_removal"
NCPU="${NCPU:-4}"                  # Parallel workers
PARALLEL="${PARALLEL:-false}"      # "true" to enable parallel
SKIP_INSTALL="${SKIP_INSTALL:-false}"
# ─────────────────────────────────────────────────────────────────────────────

# Parse CLI args
for arg in "$@"; do
  case $arg in
    --skip-install) SKIP_INSTALL=true ;;
    --parallel)     PARALLEL=true ;;
    --max_mofs=*)   MAX_MOFS="${arg#*=}" ;;
    --mode=*)       MODE="${arg#*=}" ;;
    --ncpu=*)       NCPU="${arg#*=}" ;;
  esac
done

echo "============================================================"
echo " NH3 GCMC Simulation Pipeline"
echo "============================================================"
echo "  Project dir:  $GCMC_DIR"
echo "  Mode:         $MODE"
echo "  Max MOFs:     $MAX_MOFS"
echo "  Parallel:     $PARALLEL (ncpu=$NCPU)"
echo "  Skip install: $SKIP_INSTALL"
echo "============================================================"

###############################################################################
# STEP 1: Install Dependencies
###############################################################################
if [ "$SKIP_INSTALL" = false ]; then
  echo ""
  echo ">>> STEP 1: Installing dependencies..."

  sudo apt-get update -qq
  sudo apt-get install -y -qq python3-pip build-essential > /dev/null

  pip3 install --quiet --upgrade pip

  # Core packages
  pip3 install --quiet \
    RASPA2==2.0.4 \
    pacmof2 \
    openbabel-wheel \
    pymatgen \
    ase \
    "numpy<2" \
    p_tqdm \
    python-dotenv

  echo "    All Python packages installed."
else
  echo ""
  echo ">>> STEP 1: Skipped (--skip-install)"
fi

###############################################################################
# STEP 2: Verify Installation
###############################################################################
echo ""
echo ">>> STEP 2: Verifying installation..."

# Check simulate binary
SIM_PATH=$(which simulate 2>/dev/null || true)
if [ -z "$SIM_PATH" ]; then
  echo "    ERROR: 'simulate' binary not found in PATH."
  echo "    Try: export PATH=\$HOME/.local/bin:\$PATH"
  export PATH="$HOME/.local/bin:$PATH"
  SIM_PATH=$(which simulate 2>/dev/null || true)
  if [ -z "$SIM_PATH" ]; then
    echo "    FATAL: Still not found. Install RASPA2: pip3 install RASPA2"
    exit 1
  fi
fi
echo "    simulate binary: $SIM_PATH"

# Check RASPA2 Python module
RASPA_SHARE=$(python3 -c "
import RASPA2, os
print(os.path.join(os.path.dirname(RASPA2.__file__), 'share', 'raspa'))
")
echo "    RASPA share dir: $RASPA_SHARE"

# Check PACMOF2
python3 -c "from pacmof2.pacmof2 import get_charges; print('    PACMOF2: OK')"

###############################################################################
# STEP 3: Install Custom Molecule Definitions into RASPA
###############################################################################
echo ""
echo ">>> STEP 3: Installing Ammonia/H2O molecule definitions into RASPA..."

RASPA_MOL_DIR="$RASPA_SHARE/molecules/TraPPE"
if [ -d "$RASPA_MOL_DIR" ]; then
  cp "$GCMC_DIR/molecule_definitions/TraPPE/Ammonia.def" "$RASPA_MOL_DIR/"
  cp "$GCMC_DIR/molecule_definitions/TraPPE/H2O.def"     "$RASPA_MOL_DIR/"
  echo "    Copied Ammonia.def and H2O.def to $RASPA_MOL_DIR"
else
  echo "    WARNING: $RASPA_MOL_DIR does not exist. Creating it..."
  mkdir -p "$RASPA_MOL_DIR"
  cp "$GCMC_DIR/molecule_definitions/TraPPE/Ammonia.def" "$RASPA_MOL_DIR/"
  cp "$GCMC_DIR/molecule_definitions/TraPPE/H2O.def"     "$RASPA_MOL_DIR/"
  echo "    Created and populated $RASPA_MOL_DIR"
fi

###############################################################################
# STEP 4: Fix Line Endings (Windows CRLF → Unix LF)
###############################################################################
echo ""
echo ">>> STEP 4: Fixing line endings (CRLF → LF)..."

sed -i 's/\r$//' "$GCMC_DIR/force_field/"*.def
sed -i 's/\r$//' "$GCMC_DIR/molecule_definitions/TraPPE/"*.def
echo "    Done."

###############################################################################
# STEP 5: Configure .env
###############################################################################
echo ""
echo ">>> STEP 5: Configuring .env..."

cat > "$GCMC_DIR/.env" << ENVEOF
RASPA_PATH=$RASPA_SHARE
RASPA_SIM_PATH=simulate
EGULP_PATH=
EGULP_PARAMETER_PATH=
ENVEOF
echo "    .env written with RASPA_PATH=$RASPA_SHARE"

###############################################################################
# STEP 6: Quick Smoke Test (1 MOF, 10 cycles)
###############################################################################
echo ""
echo ">>> STEP 6: Smoke test — running RASPA directly on ABEXIQ.cif (10 cycles)..."

SMOKE_DIR="/tmp/raspa_smoke_test_$$"
mkdir -p "$SMOKE_DIR"
cp "$GCMC_DIR/raw_nh3_core/cifs/ABEXIQ.cif" "$SMOKE_DIR/"
cp "$GCMC_DIR/force_field/"*                 "$SMOKE_DIR/"

cat > "$SMOKE_DIR/simulation.input" << 'SIMINPUT'
SimulationType                MonteCarlo
NumberOfCycles                10
NumberOfInitializationCycles  0
NumberOfEquilibrationCycles   10
PrintEvery                    5
Forcefield                    Local
UseChargesFromCIFFile         yes
CutOffVDW                     12
CutOffChargeCharge            12
ChargeMethod                  Ewald
EwaldPrecision                1e-6
Framework                     0
FrameworkName                 ABEXIQ
InputFileType                 cif
UnitCells                     3 2 2
HeliumVoidFraction            1.0
ExternalTemperature           298
ExternalPressure              100000
Movies                        no
Component 0 MoleculeName            Ammonia
         MoleculeDefinition            TraPPE
         TranslationProbability        0.5
         RotationProbability           0.5
         ReinsertionProbability        0.5
         SwapProbability               2.0
         CreateNumberOfMolecules       0
SIMINPUT

cd "$SMOKE_DIR"
simulate simulation.input > /dev/null 2> smoke_stderr.log

if [ -d "Output/System_0" ]; then
  echo "    SMOKE TEST PASSED — RASPA produced output successfully."
else
  echo "    SMOKE TEST FAILED — RASPA did not produce output."
  echo "    stderr:"
  cat smoke_stderr.log
  echo ""
  echo "    Common fixes:"
  echo "      1. Ensure force_field_mixing_rules.def has 'lennard-jones' keyword"
  echo "      2. Ensure line endings are Unix (LF not CRLF)"
  echo "      3. Ensure Ammonia.def is in RASPA's TraPPE molecules dir"
  rm -rf "$SMOKE_DIR"
  exit 1
fi
rm -rf "$SMOKE_DIR"

###############################################################################
# STEP 7: Run GCMC Screening
###############################################################################
echo ""
echo ">>> STEP 7: Running GCMC screening ($MAX_MOFS MOFs, mode=$MODE)..."
echo ""

cd "$GCMC_DIR"

PARALLEL_FLAG=""
if [ "$PARALLEL" = true ]; then
  PARALLEL_FLAG="--parallel --ncpu $NCPU"
fi

python3 gcmc_screen_nh3.py \
  --input ./raw_nh3_core/cifs \
  --mode "$MODE" \
  --max_mofs "$MAX_MOFS" \
  --rewrite_raspa_input \
  $PARALLEL_FLAG

echo ""
echo "============================================================"
echo " DONE — Results saved to:"
echo "   raw_nh3_core/gcmc_runs/nh3_${MODE}_results.json"
echo "============================================================"
