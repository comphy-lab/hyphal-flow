#!/bin/bash
#SBATCH --job-name=hyphaEcSweepSerial128
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --gres=tmp:100G
#SBATCH -p multi
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vatsal.sanjay@durham.ac.uk
#SBATCH --output=slurm-serial-128-%j.out
#SBATCH --error=slurm-serial-128-%j.err

# ============================================================
# Hypha Ec_h Sweep (Hamilton, serial executable, concurrent launch)
# Reads sweep values from a config file (default: sweep-128.params).
# Compiles without MPI, then launches all compiled cases in parallel.
# Usage:
#   sbatch runSweepHamilton-serial-128.sbatch [sweep_config] [--exec SOURCE]
# Examples:
#   sbatch runSweepHamilton-serial-128.sbatch
#   sbatch runSweepHamilton-serial-128.sbatch sweep-128.params
#   sbatch runSweepHamilton-serial-128.sbatch custom_sweep.params --exec hypha-capillary.c
# ============================================================

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: sbatch runSweepHamilton-serial-128.sbatch [sweep_config] [--exec SOURCE]

Arguments:
  sweep_config   Sweep config file (default: sweep-128.params)

Options:
  --exec SOURCE  Source file in simulationCases/ (default: hypha.c)
  -h, --help     Show this help message
EOF
}

# ============================================================
# Configuration
# ============================================================
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
SWEEP_CONFIG_ARG="sweep-128.params"
EXEC_SOURCE_NAME="hypha.c"
SWEEP_CONFIG_SET=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --exec)
      if [[ -z "${2:-}" ]]; then
        echo "ERROR: --exec requires a source filename." >&2
        usage
        exit 1
      fi
      EXEC_SOURCE_NAME="$2"
      shift 2
      ;;
    --exec=*)
      EXEC_SOURCE_NAME="${1#*=}"
      shift
      ;;
    --)
      shift
      break
      ;;
    -*)
      echo "ERROR: Unknown option: $1" >&2
      usage
      exit 1
      ;;
    *)
      if [[ $SWEEP_CONFIG_SET -eq 0 ]]; then
        SWEEP_CONFIG_ARG="$1"
        SWEEP_CONFIG_SET=1
        shift
      else
        echo "ERROR: Unexpected argument: $1" >&2
        usage
        exit 1
      fi
      ;;
  esac
done

if [[ $# -gt 0 ]]; then
  echo "ERROR: Unexpected trailing arguments: $*" >&2
  usage
  exit 1
fi

if [[ "$EXEC_SOURCE_NAME" != *.c ]]; then
  EXEC_SOURCE_NAME="${EXEC_SOURCE_NAME}.c"
fi

if [[ "$SWEEP_CONFIG_ARG" = /* ]]; then
  SWEEP_CONFIG="$SWEEP_CONFIG_ARG"
else
  SWEEP_CONFIG="${SCRIPT_DIR}/${SWEEP_CONFIG_ARG}"
fi

if [ ! -f "$SWEEP_CONFIG" ]; then
  echo "ERROR: Sweep config not found: $SWEEP_CONFIG" >&2
  exit 1
fi

CONFIG_DIR="$(cd "$(dirname "$SWEEP_CONFIG")" && pwd)"

# shellcheck disable=SC1090
source "$SWEEP_CONFIG"

SOURCE_FILE_NAME="$EXEC_SOURCE_NAME"
EXECUTABLE_NAME="${SOURCE_FILE_NAME%.c}"
CASE_START="${CASE_START:-1506}"
BASE_CONFIG="${BASE_CONFIG:-default.params}"
SWEEP_EC_RAW="${SWEEP_Ec_h:-}"

if [[ "$BASE_CONFIG" != /* ]]; then
  BASE_CONFIG="${CONFIG_DIR}/${BASE_CONFIG}"
fi

SRC_FILE_ORIG="${SCRIPT_DIR}/simulationCases/${SOURCE_FILE_NAME}"

if [ -z "$SWEEP_EC_RAW" ]; then
  echo "ERROR: SWEEP_Ec_h is missing in $SWEEP_CONFIG" >&2
  exit 1
fi

SWEEP_EC_RAW="$(printf '%s' "$SWEEP_EC_RAW" | tr -d '[:space:]')"
IFS=',' read -r -a _EC_VALUES <<< "$SWEEP_EC_RAW"

EC_VALUES=()
for ec in "${_EC_VALUES[@]}"; do
  [ -n "$ec" ] && EC_VALUES+=("$ec")
done

NCASES=${#EC_VALUES[@]}
if [ "$NCASES" -le 0 ]; then
  echo "ERROR: No Ec_h values parsed from SWEEP_Ec_h in $SWEEP_CONFIG" >&2
  exit 1
fi

if [ -n "${CASE_END:-}" ]; then
  EXPECTED_NCASES=$((CASE_END - CASE_START + 1))
  if [ "$EXPECTED_NCASES" -ne "$NCASES" ]; then
    echo "ERROR: CASE_START/CASE_END imply $EXPECTED_NCASES cases, but SWEEP_Ec_h has $NCASES values" >&2
    exit 1
  fi
else
  CASE_END=$((CASE_START + NCASES - 1))
fi

EC_MIN="${EC_VALUES[0]}"
EC_MAX="${EC_VALUES[$((NCASES - 1))]}"

cd "$SCRIPT_DIR"

echo "============================================="
echo "Hypha Ec_h Sweep (serial executable)"
echo "============================================="
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Working dir: $(pwd)"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Sweep config: ${SWEEP_CONFIG}"
echo "Source file: ${SOURCE_FILE_NAME}"
echo "Cases: ${CASE_START}..${CASE_END} (${NCASES})"
echo "Ec_h: ${EC_MIN} .. ${EC_MAX} (from SWEEP_Ec_h list)"
echo "============================================="
echo ""

# ============================================================
# Load Modules
# ============================================================
echo "Loading modules..."
module purge
module load gcc
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
# Prepare + Compile
# ============================================================
echo "============================================="
echo "Preparing and compiling ${NCASES} cases"
echo "============================================="
echo ""

SUCCESSFUL=0
FAILED=0
IDX=0
COMPILED_CASES=()
COMPILED_DIRS=()

for EC_H in "${EC_VALUES[@]}"; do
  IDX=$((IDX + 1))
  CASE_NO=$((CASE_START + IDX - 1))
  CASE_DIR="${SCRIPT_DIR}/simulationCases/${CASE_NO}"

  echo "---------------------------------------------"
  echo "Preparing case ${CASE_NO} (index ${IDX}/${NCASES}) Ec_h=${EC_H}"
  echo "---------------------------------------------"

  mkdir -p "$CASE_DIR"
  cp "$BASE_CONFIG" "$CASE_DIR/case.params"

  if grep -q "^CaseNo=" "$CASE_DIR/case.params"; then
    sed -i'.bak' "s|^CaseNo=.*|CaseNo=${CASE_NO}|" "$CASE_DIR/case.params"
  else
    echo "CaseNo=${CASE_NO}" >> "$CASE_DIR/case.params"
  fi

  if grep -q "^Ec_h=" "$CASE_DIR/case.params"; then
    sed -i'.bak' "s|^Ec_h=.*|Ec_h=${EC_H}|" "$CASE_DIR/case.params"
  else
    echo "Ec_h=${EC_H}" >> "$CASE_DIR/case.params"
  fi
  rm -f "$CASE_DIR/case.params.bak"

  cp "$SRC_FILE_ORIG" "$CASE_DIR/$SOURCE_FILE_NAME"

  if (
    cd "$CASE_DIR"
    qcc -I../../src-local -Wall -O2 -disable-dimensions \
      "$SOURCE_FILE_NAME" -o "$EXECUTABLE_NAME" -lm
  ); then
    COMPILED_CASES+=("$CASE_NO")
    COMPILED_DIRS+=("$CASE_DIR")
    echo "Compilation OK for case ${CASE_NO}"
  else
    echo "ERROR: Compilation failed for case ${CASE_NO}" >&2
    FAILED=$((FAILED + 1))
  fi
done

COMPILED_COUNT=${#COMPILED_CASES[@]}
if [ "$COMPILED_COUNT" -eq 0 ]; then
  echo "ERROR: No compiled cases available to run." >&2
  exit 1
fi

echo ""
echo "============================================="
echo "Launching ${COMPILED_COUNT} compiled cases in parallel"
echo "============================================="
echo ""

# Keep each case single-threaded to avoid CPU oversubscription
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

PIDS=()
PID_CASES=()

for idx in "${!COMPILED_CASES[@]}"; do
  CASE_NO="${COMPILED_CASES[$idx]}"
  CASE_DIR="${COMPILED_DIRS[$idx]}"

  (
    cd "$CASE_DIR"
    ./"$EXECUTABLE_NAME" case.params > run.out 2> run.err
  ) &

  PID=$!
  PIDS+=("$PID")
  PID_CASES+=("$CASE_NO")
  echo "Launched case ${CASE_NO} (pid ${PID})"
done

for idx in "${!PIDS[@]}"; do
  PID="${PIDS[$idx]}"
  CASE_NO="${PID_CASES[$idx]}"

  if wait "$PID"; then
    SUCCESSFUL=$((SUCCESSFUL + 1))
    echo "Case ${CASE_NO} finished OK"
  else
    FAILED=$((FAILED + 1))
    echo "ERROR: Case ${CASE_NO} failed" >&2
  fi
done

# ============================================================
# Summary
# ============================================================
echo "============================================="
echo "Sweep complete"
echo "============================================="
echo "Finished at: $(date)"
echo "Total requested: ${NCASES}"
echo "Compiled: ${COMPILED_COUNT}"
echo "Successful runs: ${SUCCESSFUL}"
echo "Failed (compile+run): ${FAILED}"
echo "Outputs in: simulationCases/"
echo "============================================="

if [ "$FAILED" -gt 0 ]; then
  exit 1
fi
exit 0
