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

# Hamilton single-case wrapper. MPI dumps are written to allocated node-local
# storage, restore-validated after a successful run, checksummed, and only then
# staged to the durable RUN_ROOT. The EXIT trap preserves partial failure
# evidence without claiming validation.

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
PARAM_FILE="${1:-default.params}"
: "${RUN_ROOT:?Set RUN_ROOT to the reserved external run directory}"
: "${TMPDIR:?Hamilton must provide allocated node-local TMPDIR}"
[[ -f "${SCRIPT_DIR}/runSimulation.sh" &&
   -f "${SCRIPT_DIR}/simulationCases/hyphal-flow.c" &&
   -f "${SCRIPT_DIR}/src-local/parse_params.sh" ]] || {
  printf 'ERROR: submit this wrapper from the hyphal-flow repository root\n' >&2
  exit 2
}
[[ "${PARAM_FILE}" = /* ]] || PARAM_FILE="${SCRIPT_DIR}/${PARAM_FILE}"
[[ -f "${PARAM_FILE}" ]] || {
  printf 'ERROR: parameter file not found: %s\n' "${PARAM_FILE}" >&2
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
[[ "${TMPDIR}" = /* && -d "${TMPDIR}" && -w "${TMPDIR}" ]] || {
  printf 'ERROR: TMPDIR must be an allocated writable absolute directory\n' >&2
  exit 2
}
command -v rsync >/dev/null 2>&1 || {
  printf 'ERROR: rsync is required for verified stage-out\n' >&2
  exit 2
}

# shellcheck source=src-local/parse_params.sh
source "${SCRIPT_DIR}/src-local/parse_params.sh"
CASE_NO="$(get_param_value CaseNo "${PARAM_FILE}")"
[[ "${CASE_NO}" =~ ^[0-9]{4}$ ]] || {
  printf 'ERROR: invalid CaseNo in %s\n' "${PARAM_FILE}" >&2
  exit 2
}

LOCAL_OUTPUT_ROOT="${TMPDIR}/hyphal-flow-${SLURM_JOB_ID}"
LOCAL_CASE_DIR="${LOCAL_OUTPUT_ROOT}/${CASE_NO}"
DURABLE_CASE_DIR="${RUN_REAL}/${CASE_NO}"
[[ ! -e "${DURABLE_CASE_DIR}" ]] || {
  printf 'ERROR: durable case path already exists: %s\n' "${DURABLE_CASE_DIR}" >&2
  exit 2
}
mkdir -p "${LOCAL_OUTPUT_ROOT}"

stage_out () {
  status=$?
  trap - EXIT
  if [[ -d "${LOCAL_CASE_DIR}" ]]; then
    {
      printf 'job_id=%s\n' "${SLURM_JOB_ID}"
      printf 'case=%s\n' "${CASE_NO}"
      printf 'exit_status=%s\n' "${status}"
      printf 'staged_at=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } > "${LOCAL_CASE_DIR}/stage-out-status.txt"
    checksum_file="${TMPDIR}/stage-out-${SLURM_JOB_ID}-${CASE_NO}.sha256"
    (
      cd "${LOCAL_CASE_DIR}"
      find . -type f ! -name stage-out-sha256.txt -print0 |
        sort -z | xargs -0 sha256sum
    ) > "${checksum_file}"
    mv "${checksum_file}" "${LOCAL_CASE_DIR}/stage-out-sha256.txt"
    mkdir "${DURABLE_CASE_DIR}"
    rsync -a --partial --delay-updates -- \
      "${LOCAL_CASE_DIR}/" "${DURABLE_CASE_DIR}/"
  fi
  exit "${status}"
}
trap stage_out EXIT

module purge
module load gcc
module load openmpi

cd "$SCRIPT_DIR"
bash runSimulation.sh "$PARAM_FILE" \
  --output-root "$LOCAL_OUTPUT_ROOT" --mpi --cpus "${SLURM_NTASKS:-128}"

if command -v qcc >/dev/null 2>&1; then
  QCC="$(command -v qcc)"
elif [[ -n "${BASILISK:-}" && -x "${BASILISK}/qcc" ]]; then
  QCC="${BASILISK}/qcc"
else
  printf 'ERROR: qcc is unavailable for dump validation\n' >&2
  exit 2
fi

FACET_CHECK="${TMPDIR}/getFacet-threePhase"
(
  cd "${SCRIPT_DIR}/postProcess"
  "${QCC}" -O2 -Wall -disable-dimensions getFacet-threePhase.c \
    -o "${FACET_CHECK}" -lm
)

(
  cd "${LOCAL_CASE_DIR}"
  while IFS= read -r -d '' dump_file; do
    relative="${dump_file#./}"
    "${FACET_CHECK}" "${relative}" true 2>/dev/null
    "${FACET_CHECK}" "${relative}" false 2>/dev/null
  done < <(
    {
      find intermediate -maxdepth 1 -type f -name 'snapshot-*' -print0
      printf '%s\0' restart final
    } | sort -z -V
  )
  printf 'validated_at=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    > dump-validation.txt
)
