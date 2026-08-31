#!/usr/bin/env bash

# Generate a Cartesian product of SWEEP_* values and route every case through
# the canonical single-case runner. Dry runs only parse configuration; they do
# not require Basilisk or a machine-local project configuration.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=src-local/parse_params.sh
source "${SCRIPT_DIR}/src-local/parse_params.sh"

usage() {
  command cat <<'EOF'
Usage: bash runParameterSweep.sh [sweep_file] [OPTIONS]

Options:
  --exec FILE       Simulation source (only hyphal-flow.c is supported)
  --output-root DIR Case-output root passed to runSimulation.sh
  --dry-run         Print generated cases without compiling or executing
  --verbose         Print all expanded parameter values
  --mpi             Compile and execute every case with MPI
  --cpus N          MPI process count (alias: --CPUs; default: 4)
  -h, --help        Show this help
EOF
}

EXEC_CODE="hyphal-flow.c"
SWEEP_FILE="sweep.params"
SWEEP_FILE_SET=0
OUTPUT_ROOT="${SCRIPT_DIR}/simulationCases"
DRY_RUN=0
VERBOSE=0
USE_MPI=0
MPI_CPUS=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --exec) EXEC_CODE="${2:-}"; shift 2 ;;
    --exec=*) EXEC_CODE="${1#*=}"; shift ;;
    --output-root) OUTPUT_ROOT="${2:-}"; shift 2 ;;
    --output-root=*) OUTPUT_ROOT="${1#*=}"; shift ;;
    -n|--dry-run) DRY_RUN=1; shift ;;
    -v|--verbose) VERBOSE=1; shift ;;
    --mpi) USE_MPI=1; shift ;;
    --cpus|--CPUs) MPI_CPUS="${2:-}"; shift 2 ;;
    --cpus=*|--CPUs=*) MPI_CPUS="${1#*=}"; shift ;;
    --) shift; break ;;
    -*) printf 'ERROR: unknown option: %s\n' "$1" >&2; usage; exit 2 ;;
    *)
      if [[ $SWEEP_FILE_SET -ne 0 ]]; then
        printf 'ERROR: unexpected argument: %s\n' "$1" >&2
        exit 2
      fi
      SWEEP_FILE="$1"
      SWEEP_FILE_SET=1
      shift
      ;;
  esac
done

if [[ $# -gt 0 ]]; then
  printf 'ERROR: unexpected trailing arguments: %s\n' "$*" >&2
  exit 2
fi
if [[ "$EXEC_CODE" != "hyphal-flow.c" ]]; then
  printf 'ERROR: only hyphal-flow.c is sweep-compatible; got %s\n' "$EXEC_CODE" >&2
  exit 2
fi
if [[ -z "$OUTPUT_ROOT" ]]; then
  printf 'ERROR: --output-root may not be empty\n' >&2
  exit 2
fi
if [[ ! "$MPI_CPUS" =~ ^[1-9][0-9]*$ ]]; then
  printf 'ERROR: --cpus must be a positive integer; got %s\n' "$MPI_CPUS" >&2
  exit 2
fi

[[ "$SWEEP_FILE" = /* ]] || SWEEP_FILE="${SCRIPT_DIR}/${SWEEP_FILE}"
[[ "$OUTPUT_ROOT" = /* ]] || OUTPUT_ROOT="${SCRIPT_DIR}/${OUTPUT_ROOT}"
[[ -f "$SWEEP_FILE" ]] || { printf 'ERROR: sweep file not found: %s\n' "$SWEEP_FILE" >&2; exit 2; }

CONFIG_DIR="$(cd "$(dirname "$SWEEP_FILE")" && pwd)"
BASE_CONFIG="$(get_param_value BASE_CONFIG "$SWEEP_FILE" default.params)"
CASE_START="$(get_param_value CASE_START "$SWEEP_FILE" 1000)"
CASE_END="$(get_param_value CASE_END "$SWEEP_FILE" '')"
[[ "$BASE_CONFIG" = /* ]] || BASE_CONFIG="${CONFIG_DIR}/${BASE_CONFIG}"
[[ -f "$BASE_CONFIG" ]] || { printf 'ERROR: base config not found: %s\n' "$BASE_CONFIG" >&2; exit 2; }

if [[ ! "$CASE_START" =~ ^[0-9]{4}$ ]] || ((10#$CASE_START < 1000 || 10#$CASE_START > 9999)); then
  printf 'ERROR: CASE_START must be four digits from 1000 to 9999; got %s\n' "$CASE_START" >&2
  exit 2
fi
if [[ -n "$CASE_END" ]] && { [[ ! "$CASE_END" =~ ^[0-9]{4}$ ]] || ((10#$CASE_END < 1000 || 10#$CASE_END > 9999)); }; then
  printf 'ERROR: CASE_END must be four digits from 1000 to 9999; got %s\n' "$CASE_END" >&2
  exit 2
fi

SWEEP_VARS=()
SWEEP_VALUES=()
while IFS= read -r line || [[ -n "$line" ]]; do
  line="${line%%#*}"
  line="$(trim_string "$line")"
  [[ -z "$line" || "$line" != *=* ]] && continue
  key="$(trim_string "${line%%=*}")"
  value="$(trim_string "${line#*=}")"
  if [[ "$key" =~ ^SWEEP_([A-Za-z_][A-Za-z0-9_]*)$ ]]; then
    if [[ "${BASH_REMATCH[1]}" == CaseNo ]]; then
      printf 'ERROR: SWEEP_CaseNo is reserved for deterministic case numbering\n' >&2
      exit 2
    fi
    [[ -n "$value" ]] || { printf 'ERROR: %s has no values\n' "$key" >&2; exit 2; }
    SWEEP_VARS+=("${BASH_REMATCH[1]}")
    SWEEP_VALUES+=("$value")
  fi
done < "$SWEEP_FILE"
[[ ${#SWEEP_VARS[@]} -gt 0 ]] || { printf 'ERROR: no SWEEP_* values in %s\n' "$SWEEP_FILE" >&2; exit 2; }

TEMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/hyphal-sweep.XXXXXX")"
trap 'rm -rf "$TEMP_DIR"' EXIT
PARAM_FILES=()
COMBINATION_COUNT=0
CURRENT_CASE=$((10#$CASE_START))

generate_combinations() {
  local depth="$1"
  shift
  local current=("$@")
  local index raw value case_file
  local values=()

  if [[ "$depth" -eq "${#SWEEP_VARS[@]}" ]]; then
    if ((CURRENT_CASE > 9999)); then
      printf 'ERROR: generated CaseNo exceeds 9999\n' >&2
      exit 2
    fi
    case_file="${TEMP_DIR}/case_$(printf '%04d' "$CURRENT_CASE").params"
    cp "$BASE_CONFIG" "$case_file"
    set_param_in_file CaseNo "$(printf '%04d' "$CURRENT_CASE")" "$case_file"
    for index in "${!SWEEP_VARS[@]}"; do
      set_param_in_file "${SWEEP_VARS[$index]}" "${current[$index]}" "$case_file"
    done
    PARAM_FILES+=("$case_file")
    if [[ $DRY_RUN -eq 1 || $VERBOSE -eq 1 ]]; then
      printf 'Case %04d:' "$CURRENT_CASE"
      for index in "${!SWEEP_VARS[@]}"; do
        printf ' %s=%s' "${SWEEP_VARS[$index]}" "${current[$index]}"
      done
      printf '\n'
    fi
    ((COMBINATION_COUNT += 1))
    ((CURRENT_CASE += 1))
    return
  fi

  IFS=',' read -r -a values <<< "${SWEEP_VALUES[$depth]}"
  for raw in "${values[@]}"; do
    value="$(trim_string "$raw")"
    [[ -n "$value" ]] || continue
    if [[ ${#current[@]} -gt 0 ]]; then
      generate_combinations "$((depth + 1))" "${current[@]}" "$value"
    else
      generate_combinations "$((depth + 1))" "$value"
    fi
  done
}

generate_combinations 0
((COMBINATION_COUNT > 0)) || { printf 'ERROR: no combinations generated\n' >&2; exit 2; }
GENERATED_END=$((10#$CASE_START + COMBINATION_COUNT - 1))
if [[ -n "$CASE_END" ]] && ((10#$CASE_END != GENERATED_END)); then
  printf 'ERROR: CASE_START/CASE_END imply %d cases but %d combinations were generated\n' \
    "$((10#$CASE_END - 10#$CASE_START + 1))" "$COMBINATION_COUNT" >&2
  exit 2
fi

printf 'Sweep: %d case(s), %04d..%04d, source=%s\n' \
  "$COMBINATION_COUNT" "$((10#$CASE_START))" "$GENERATED_END" "$EXEC_CODE"
if [[ $DRY_RUN -eq 1 ]]; then
  printf 'Dry run complete; no simulations executed.\n'
  exit 0
fi

SUCCESSFUL=0
FAILED=0
for param_file in "${PARAM_FILES[@]}"; do
  command=(bash "${SCRIPT_DIR}/runSimulation.sh" "$param_file"
           --exec "$EXEC_CODE" --output-root "$OUTPUT_ROOT")
  if [[ $USE_MPI -eq 1 ]]; then
    command+=(--mpi --cpus "$MPI_CPUS")
  fi
  if "${command[@]}"; then
    ((SUCCESSFUL += 1))
  else
    ((FAILED += 1))
  fi
done

printf 'Sweep complete: %d passed, %d failed.\n' "$SUCCESSFUL" "$FAILED"
((FAILED == 0))
