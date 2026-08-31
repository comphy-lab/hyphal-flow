#!/bin/bash
#SBATCH --job-name=hyphaSingleMPI
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -c 1
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --gres=tmp:100G
#SBATCH -p multi
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@durham.ac.uk
#SBATCH --output=slurm-single-%j.out
#SBATCH --error=slurm-single-%j.err

# ============================================================
# Hypha Single Case (Hamilton, MPI)
# Compiles one case and runs it with MPI.
# Usage:
#   sbatch runSingleCaseHamilton.sbatch [params_file] [--exec SOURCE]
# Examples:
#   sbatch runSingleCaseHamilton.sbatch
#   sbatch runSingleCaseHamilton.sbatch default.params
#   sbatch runSingleCaseHamilton.sbatch default.params --exec hypha-capillary.c
# ============================================================

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: sbatch runSingleCaseHamilton.sbatch [params_file] [--exec SOURCE]

Arguments:
  params_file    Parameter file (default: default.params)

Options:
  --exec SOURCE  Source file in simulationCases/ (default: hypha.c)
  -h, --help     Show this help message
EOF
}

get_param_value() {
  local key="$1"
  local file="$2"
  awk -F '=' -v key="$key" '
    /^[[:space:]]*#/ { next }
    {
      k = $1
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", k)
      if (k == key) {
        v = $2
        sub(/[[:space:]]*#.*/, "", v)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", v)
        print v
        exit
      }
    }
  ' "$file"
}

# ============================================================
# Configuration
# ============================================================
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
PARAM_FILE_ARG="default.params"
EXEC_SOURCE_NAME="hypha.c"
PARAM_FILE_SET=0

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
      if [[ $PARAM_FILE_SET -eq 0 ]]; then
        PARAM_FILE_ARG="$1"
        PARAM_FILE_SET=1
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

if [[ "$PARAM_FILE_ARG" = /* ]]; then
  PARAM_FILE="$PARAM_FILE_ARG"
else
  PARAM_FILE="${SCRIPT_DIR}/${PARAM_FILE_ARG}"
fi

SOURCE_FILE_NAME="$EXEC_SOURCE_NAME"
EXECUTABLE_NAME="${SOURCE_FILE_NAME%.c}"
MPI_TASKS="${SLURM_NTASKS:-128}"
SRC_FILE_ORIG="${SCRIPT_DIR}/simulationCases/${SOURCE_FILE_NAME}"

cd "$SCRIPT_DIR"

echo "============================================="
echo "Hypha Single Case (MPI)"
echo "============================================="
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Working dir: $(pwd)"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Params file: ${PARAM_FILE}"
echo "Source file: ${SOURCE_FILE_NAME}"
echo "MPI tasks: ${MPI_TASKS}"
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
  # shellcheck disable=SC1090
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
if [ ! -f "$PARAM_FILE" ]; then
  echo "ERROR: Parameter file not found: $PARAM_FILE" >&2
  exit 1
fi

if [ ! -f "$SRC_FILE_ORIG" ]; then
  echo "ERROR: Source file not found: $SRC_FILE_ORIG" >&2
  exit 1
fi

CASE_NO="$(get_param_value "CaseNo" "$PARAM_FILE")"
if [[ -z "$CASE_NO" ]]; then
  echo "ERROR: CaseNo not found in parameter file: $PARAM_FILE" >&2
  exit 1
fi

if [[ ! "$CASE_NO" =~ ^[0-9]+$ ]]; then
  echo "ERROR: CaseNo must be numeric, got: $CASE_NO" >&2
  exit 1
fi

CASE_DIR="${SCRIPT_DIR}/simulationCases/${CASE_NO}"

# ============================================================
# Prepare + Compile
# ============================================================
echo "Preparing case directory: ${CASE_DIR}"
mkdir -p "$CASE_DIR"
cp "$PARAM_FILE" "$CASE_DIR/case.params"
cp "$SRC_FILE_ORIG" "$CASE_DIR/$SOURCE_FILE_NAME"

# Ensure CaseNo in case.params matches target folder.
if grep -q "^CaseNo=" "$CASE_DIR/case.params"; then
  sed -i'.bak' "s|^CaseNo=.*|CaseNo=${CASE_NO}|" "$CASE_DIR/case.params"
else
  echo "CaseNo=${CASE_NO}" >> "$CASE_DIR/case.params"
fi
rm -f "$CASE_DIR/case.params.bak"

cd "$CASE_DIR"

echo "Compiling..."
if CC99='mpicc -std=c99 -D_GNU_SOURCE=1' qcc -I../../src-local \
    -Wall -O2 -D_MPI=1 -disable-dimensions \
    "$SOURCE_FILE_NAME" -o "$EXECUTABLE_NAME" -lm 2>&1; then
  echo "Compilation OK"
else
  echo "ERROR: Compilation failed for case ${CASE_NO}" >&2
  exit 1
fi
echo ""

if [ -f "restart" ]; then
  echo "Restart file found - will resume"
fi

# ============================================================
# Run
# ============================================================
echo "Running: mpirun -np ${MPI_TASKS} ./${EXECUTABLE_NAME} case.params"
if mpirun -np "${MPI_TASKS}" ./"$EXECUTABLE_NAME" case.params; then
  echo ""
  echo "Case ${CASE_NO} finished OK"
  echo "Output directory: simulationCases/${CASE_NO}/"
  exit 0
else
  code=$?
  echo ""
  echo "ERROR: Case ${CASE_NO} failed (exit ${code})" >&2
  exit "$code"
fi
