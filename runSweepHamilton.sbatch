#!/bin/bash
#SBATCH --job-name=hyphaEcSweep
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -c 1
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --gres=tmp:100G
#SBATCH -p multi
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@durham.ac.uk
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# ============================================================
# Hypha Ec_h Sweep (Hamilton): 200 cases, log-spaced Ec_h
# Ec_h from 0.01 to 10 (inclusive), only Ec_h changes.
# Each case compiled + run sequentially with MPI.
# ============================================================

set -euo pipefail

# ============================================================
# Configuration
# ============================================================
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
SOURCE_FILE_NAME="dropImpact.c"
EXECUTABLE_NAME="${SOURCE_FILE_NAME%.c}"
MPI_TASKS="${SLURM_NTASKS:-128}"

# Case numbering (200 cases)
CASE_START=1000
NCASES=200
CASE_END=$((CASE_START + NCASES - 1))

# Sweep range (log-spaced, inclusive)
EC_MIN="0.01"
EC_MAX="10"

# Base params file: keep everything fixed here
# Must contain at least CaseNo=... (will be overwritten)
BASE_CONFIG="${SCRIPT_DIR}/sweep_base.params"

# Where your source lives (as in your original script)
SRC_FILE_ORIG="${SCRIPT_DIR}/simulationCases/${SOURCE_FILE_NAME}"

cd "$SCRIPT_DIR"

echo "============================================="
echo "Hypha Ec_h Sweep"
echo "============================================="
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Working dir: $(pwd)"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "MPI tasks: ${MPI_TASKS}"
echo "Cases: ${CASE_START}..${CASE_END} (${NCASES})"
echo "Ec_h: ${EC_MIN} .. ${EC_MAX} (log-spaced, inclusive)"
echo "============================================="
echo ""

# ============================================================
# Load Modules
# ============================================================
echo "Loading modules..."
module purge
module load gcc
module load openmpi
echo "Modules loaded"
echo ""

# ============================================================
# Setup Basilisk Environment (ref-locked)
# ============================================================
echo "Setting up Basilisk environment..."
if ! command -v curl &> /dev/null; then
  echo "ERROR: curl not found" >&2
  exit 1
fi

BASILISK_REF="v2026-01-13"
BASILISK_INSTALL_URL="https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh"
LOCK_FILE="${SCRIPT_DIR}/basilisk/.comphy-lock"
LOCK_REF=""
if [ -f "$LOCK_FILE" ]; then
  LOCK_REF="$(awk -F= '$1=="ref"{print $2; exit}' "$LOCK_FILE" 2>/dev/null || true)"
fi

if [ ! -f "$LOCK_FILE" ] || [ "$LOCK_REF" != "$BASILISK_REF" ]; then
  echo "Lock mismatch or missing (found ref='${LOCK_REF:-none}') - hard install..."
  curl -fsSL "$BASILISK_INSTALL_URL" | bash -s -- --ref="$BASILISK_REF" --hard
else
  curl -fsSL "$BASILISK_INSTALL_URL" | bash -s -- --ref="$BASILISK_REF"
fi

if [ -f "${SCRIPT_DIR}/.project_config" ]; then
  source "${SCRIPT_DIR}/.project_config"
  echo "Basilisk env loaded from .project_config"
  echo "BASILISK: ${BASILISK:-unset}"
else
  echo "ERROR: .project_config not found after Basilisk install" >&2
  exit 1
fi
echo ""

# ============================================================
# Validate inputs
# ============================================================
if [ ! -f "$BASE_CONFIG" ]; then
  echo "ERROR: Base params file not found: $BASE_CONFIG" >&2
  echo "Create it with your fixed parameters (everything except Ec_h)." >&2
  exit 1
fi

if [ ! -f "$SRC_FILE_ORIG" ]; then
  echo "ERROR: Source file not found: $SRC_FILE_ORIG" >&2
  exit 1
fi

# ============================================================
# Helper: compute log-spaced Ec_h list using python
# ============================================================
EC_LIST_FILE="${SCRIPT_DIR}/.Ec_h_values_${SLURM_JOB_ID:-$$}.txt"
trap 'rm -f "$EC_LIST_FILE"' EXIT

python3 - << 'PY' > "$EC_LIST_FILE"
import numpy as np
ec_min = 0.01
ec_max = 10.0
n = 200
vals = np.logspace(np.log10(ec_min), np.log10(ec_max), n)
for v in vals:
    # plenty of sig figs; Basilisk atof can read scientific notation
    print(f"{v:.12g}")
PY

# ============================================================
# Run sweep
# ============================================================
echo "============================================="
echo "Running ${NCASES} simulations sequentially"
echo "============================================="
echo ""

SUCCESSFUL=0
FAILED=0
CASE_NO=$CASE_START
IDX=0

while read -r EC_H; do
  IDX=$((IDX + 1))
  CASE_DIR="${SCRIPT_DIR}/simulationCases/${CASE_NO}"

  echo "========================================="
  echo "Case ${CASE_NO}  (index ${IDX}/${NCASES})  Ec_h=${EC_H}"
  echo "Time: $(date)"
  echo "========================================="

  mkdir -p "$CASE_DIR"

  # Write case.params (base + overrides)
  cp "$BASE_CONFIG" "$CASE_DIR/case.params"

  # Ensure CaseNo is set/updated
  if grep -q "^CaseNo=" "$CASE_DIR/case.params"; then
    sed -i'.bak' "s|^CaseNo=.*|CaseNo=${CASE_NO}|" "$CASE_DIR/case.params"
  else
    echo "CaseNo=${CASE_NO}" >> "$CASE_DIR/case.params"
  fi

  # Set/override Ec_h only
  if grep -q "^Ec_h=" "$CASE_DIR/case.params"; then
    sed -i'.bak' "s|^Ec_h=.*|Ec_h=${EC_H}|" "$CASE_DIR/case.params"
  else
    echo "Ec_h=${EC_H}" >> "$CASE_DIR/case.params"
  fi
  rm -f "$CASE_DIR/case.params.bak"

  # Copy source into case directory (so it’s self-contained)
  cp "$SRC_FILE_ORIG" "$CASE_DIR/$SOURCE_FILE_NAME"

  cd "$CASE_DIR"

  # Compile
  echo "Compiling..."
  if CC99='mpicc -std=c99 -D_GNU_SOURCE=1' qcc -I../../src-local \
      -Wall -O2 -D_MPI=1 -disable-dimensions \
      "$SOURCE_FILE_NAME" -o "$EXECUTABLE_NAME" -lm 2>&1; then
    echo "Compilation OK"
  else
    echo "ERROR: Compilation failed for case ${CASE_NO}" >&2
    FAILED=$((FAILED + 1))
    cd "$SCRIPT_DIR"
    CASE_NO=$((CASE_NO + 1))
    continue
  fi

  if [ -f "restart" ]; then
    echo "Restart file found - will resume"
  fi

  # Run
  echo "Running: mpirun -np ${MPI_TASKS} ./${EXECUTABLE_NAME} case.params"
  if mpirun -np "${MPI_TASKS}" ./"$EXECUTABLE_NAME" case.params; then
    echo "Case ${CASE_NO} finished OK"
    SUCCESSFUL=$((SUCCESSFUL + 1))
  else
    code=$?
    echo "ERROR: Case ${CASE_NO} failed (exit ${code})" >&2
    FAILED=$((FAILED + 1))
  fi

  cd "$SCRIPT_DIR"
  echo ""

  CASE_NO=$((CASE_NO + 1))
done < "$EC_LIST_FILE"

# ============================================================
# Summary
# ============================================================
echo "============================================="
echo "Sweep complete"
echo "============================================="
echo "Finished at: $(date)"
echo "Total: ${NCASES}"
echo "Successful: ${SUCCESSFUL}"
echo "Failed: ${FAILED}"
echo "Outputs in: simulationCases/"
echo "============================================="

if [ "$FAILED" -gt 0 ]; then
  exit 1
fi
exit 0
