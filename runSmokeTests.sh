#!/usr/bin/env bash

# CoMPhy smoke campaign: compile and advance reduced physical simulations.
# This is an execution/restart contract, not unit testing, verification or
# validation of the constitutive model.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  command cat <<'EOF'
Usage: bash runSmokeTests.sh --output-root DIR [--skip-mpi]

Runs four reduced serial simulations (Newtonian control, solid-only,
liquids-only, and full rheology), a two-case full-rheology sweep, a two-rank
MPI full-rheology run, and an input-matched restart of that MPI case.

The output root must be outside the source checkout for managed compute runs.
EOF
}

require_value() {
  if [[ $# -lt 2 || -z "${2:-}" ]]; then
    printf 'ERROR: %s requires a value\n' "$1" >&2
    exit 2
  fi
}

OUTPUT_ROOT=""
SKIP_MPI=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --output-root) require_value "$1" "${2:-}"; OUTPUT_ROOT="$2"; shift 2 ;;
    --output-root=*) OUTPUT_ROOT="${1#*=}"; shift ;;
    --skip-mpi) SKIP_MPI=1; shift ;;
    *) printf 'ERROR: unknown argument: %s\n' "$1" >&2; usage; exit 2 ;;
  esac
done

[[ -n "$OUTPUT_ROOT" ]] || { printf 'ERROR: --output-root is required\n' >&2; exit 2; }
[[ "$OUTPUT_ROOT" = /* ]] || OUTPUT_ROOT="${PWD}/${OUTPUT_ROOT}"
mkdir -p "$OUTPUT_ROOT"
SCRIPT_REAL="$(cd "$SCRIPT_DIR" && pwd -P)"
OUTPUT_REAL="$(cd "$OUTPUT_ROOT" && pwd -P)"
case "${OUTPUT_REAL}/" in
  "${SCRIPT_REAL}/"*)
    printf 'ERROR: --output-root must be outside the source checkout\n' >&2
    exit 2
    ;;
esac
OUTPUT_ROOT="$OUTPUT_REAL"

for case_name in newtonian solid-only liquids-only full; do
  bash "${SCRIPT_DIR}/runSimulation.sh" \
    "${SCRIPT_DIR}/smoke/${case_name}.params" \
    --output-root "${OUTPUT_ROOT}/serial"
done

bash "${SCRIPT_DIR}/runParameterSweep.sh" \
  "${SCRIPT_DIR}/smoke/sweep.params" \
  --output-root "${OUTPUT_ROOT}/sweep"

if [[ $SKIP_MPI -eq 0 ]]; then
  MPI_CASE="${OUTPUT_ROOT}/mpi/9103"
  bash "${SCRIPT_DIR}/runSimulation.sh" "${SCRIPT_DIR}/smoke/full.params" \
    --output-root "${OUTPUT_ROOT}/mpi" --mpi --cpus 2

  [[ -s "${MPI_CASE}/intermediate/snapshot-0.0100" ]] || {
    printf 'ERROR: restart smoke needs intermediate/snapshot-0.0100\n' >&2
    exit 1
  }
  cp "${MPI_CASE}/log" "${MPI_CASE}/log.initial"
  awk 'NR == 1 || $3 <= 0.0100001' "${MPI_CASE}/log.initial" > "${MPI_CASE}/log"
  mv "${MPI_CASE}/final" "${MPI_CASE}/final.initial"
  cp "${MPI_CASE}/intermediate/snapshot-0.0100" "${MPI_CASE}/restart"
  bash "${SCRIPT_DIR}/runSimulation.sh" "${SCRIPT_DIR}/smoke/full.params" \
    --output-root "${OUTPUT_ROOT}/mpi" --mpi --cpus 2 --resume
fi

FULL_LOG="${OUTPUT_ROOT}/serial/9103/log"
awk '
  NR == 1 { next }
  {
    seen = 1
    if ($3 > 0) positive_time = 1
    if ($6 > max_overlap) max_overlap = $6
    if ($10 > 0) drop_stress = 1
    if ($11 > 0) solid_stress = 1
    if ($12 > 0) liquid_stress = 1
  }
  END {
    if (!seen || !positive_time || !drop_stress || !solid_stress ||
        !liquid_stress || max_overlap > 1e-6)
      exit 1
  }
' "$FULL_LOG" || {
  printf 'ERROR: full-rheology diagnostic gate failed (time/stress/overlap)\n' >&2
  exit 1
}

printf 'CoMPhy smoke campaign passed: %s\n' "$OUTPUT_ROOT"
