#!/usr/bin/env bash
#SBATCH --job-name=hyphal-flow
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --gres=tmp:100G
#SBATCH --partition=multi
#SBATCH --output=slurm-single-%j.out
#SBATCH --error=slurm-single-%j.err

# Thin Hamilton wrapper. The site environment owns the pinned Basilisk
# installation; this job never downloads or installs software.

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
PARAM_FILE="${1:-default.params}"
: "${RUN_ROOT:?Set RUN_ROOT to the reserved external run directory}"
[[ -f "${SCRIPT_DIR}/runSimulation.sh" &&
   -f "${SCRIPT_DIR}/simulationCases/hyphal-flow.c" ]] || {
  printf 'ERROR: submit this wrapper from the hyphal-flow repository root\n' >&2
  exit 2
}
[[ "$RUN_ROOT" = /* && "$RUN_ROOT" != / && -d "$RUN_ROOT" ]] || {
  printf 'ERROR: RUN_ROOT must be an existing absolute directory other than /\n' >&2
  exit 2
}
SCRIPT_REAL="$(cd "$SCRIPT_DIR" && pwd -P)"
RUN_REAL="$(cd "$RUN_ROOT" && pwd -P)"
case "${RUN_REAL}/" in
  "${SCRIPT_REAL}/"*) printf 'ERROR: RUN_ROOT must be outside the source checkout\n' >&2; exit 2 ;;
esac

module purge
module load gcc
module load openmpi

cd "$SCRIPT_DIR"
bash runSimulation.sh "$PARAM_FILE" \
  --output-root "$RUN_REAL" --mpi --cpus "${SLURM_NTASKS:-128}"
