#!/usr/bin/env bash

# Compile and run one canonical hyphal-flow case. Runtime output is isolated
# from source; an existing case is never reused without an explicit,
# input-matched --resume.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=src-local/parse_params.sh
source "${SCRIPT_DIR}/src-local/parse_params.sh"

usage() {
  command cat <<'EOF'
Usage: bash runSimulation.sh [params_file] [OPTIONS]

Options:
  --exec FILE       Simulation source (only hyphal-flow.c is supported)
  --output-root DIR Case-output root (default: simulationCases)
  --compile-only    Compile the case but do not execute it
  --dry-run         Validate and print the intended compile/run without writes
  --resume          Resume only when source and parameter hashes match
  --mpi             Compile and execute with MPI
  --cpus N          MPI process count (alias: --CPUs; default: 4)
  -h, --help        Show this help
EOF
}

EXEC_CODE="hyphal-flow.c"
PARAM_FILE="default.params"
PARAM_FILE_SET=0
OUTPUT_ROOT="${SCRIPT_DIR}/simulationCases"
COMPILE_ONLY=0
DRY_RUN=0
RESUME=0
USE_MPI=0
MPI_CPUS=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --exec) EXEC_CODE="${2:-}"; shift 2 ;;
    --exec=*) EXEC_CODE="${1#*=}"; shift ;;
    --output-root) OUTPUT_ROOT="${2:-}"; shift 2 ;;
    --output-root=*) OUTPUT_ROOT="${1#*=}"; shift ;;
    --compile-only) COMPILE_ONLY=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    --resume) RESUME=1; shift ;;
    --mpi) USE_MPI=1; shift ;;
    --cpus|--CPUs) MPI_CPUS="${2:-}"; shift 2 ;;
    --cpus=*|--CPUs=*) MPI_CPUS="${1#*=}"; shift ;;
    --) shift; break ;;
    -*) printf 'ERROR: unknown option: %s\n' "$1" >&2; usage; exit 2 ;;
    *)
      if [[ $PARAM_FILE_SET -ne 0 ]]; then
        printf 'ERROR: unexpected argument: %s\n' "$1" >&2
        exit 2
      fi
      PARAM_FILE="$1"
      PARAM_FILE_SET=1
      shift
      ;;
  esac
done

if [[ $# -gt 0 ]]; then
  printf 'ERROR: unexpected trailing arguments: %s\n' "$*" >&2
  exit 2
fi
if [[ "$EXEC_CODE" != "hyphal-flow.c" ]]; then
  printf 'ERROR: only simulationCases/hyphal-flow.c accepts the runtime contract; got %s\n' "$EXEC_CODE" >&2
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

[[ "$PARAM_FILE" = /* ]] || PARAM_FILE="${SCRIPT_DIR}/${PARAM_FILE}"
[[ "$OUTPUT_ROOT" = /* ]] || OUTPUT_ROOT="${SCRIPT_DIR}/${OUTPUT_ROOT}"
SOURCE_FILE="${SCRIPT_DIR}/simulationCases/${EXEC_CODE}"

[[ -f "$PARAM_FILE" ]] || { printf 'ERROR: parameter file not found: %s\n' "$PARAM_FILE" >&2; exit 2; }
[[ -f "$SOURCE_FILE" ]] || { printf 'ERROR: source file not found: %s\n' "$SOURCE_FILE" >&2; exit 2; }

CASE_NO="$(get_param_value CaseNo "$PARAM_FILE")"
if [[ ! "$CASE_NO" =~ ^[0-9]{4}$ ]] || ((10#$CASE_NO < 1000 || 10#$CASE_NO > 9999)); then
  printf 'ERROR: CaseNo must be a four-digit value from 1000 to 9999; got %s\n' "${CASE_NO:-<missing>}" >&2
  exit 2
fi

# A local configuration may add qcc to PATH, but a clean clone does not need it.
if [[ -f "${SCRIPT_DIR}/.project_config" ]]; then
  # shellcheck disable=SC1091
  source "${SCRIPT_DIR}/.project_config"
fi
if command -v qcc >/dev/null 2>&1; then
  QCC="$(command -v qcc)"
elif [[ -n "${BASILISK:-}" && -x "${BASILISK}/qcc" ]]; then
  QCC="${BASILISK}/qcc"
else
  printf 'ERROR: qcc is not available in PATH or BASILISK/qcc\n' >&2
  exit 2
fi

QCC_DIR="$(dirname "$QCC")"
if grep -q 'void set_prolongation' "${QCC_DIR}/grid/multigrid-common.h" 2>/dev/null; then
  BASILISK_API=modern
else
  BASILISK_API=legacy
fi

if [[ $USE_MPI -eq 1 ]]; then
  command -v mpicc >/dev/null 2>&1 || { printf 'ERROR: mpicc is required for --mpi\n' >&2; exit 2; }
  command -v mpirun >/dev/null 2>&1 || { printf 'ERROR: mpirun is required for --mpi\n' >&2; exit 2; }
fi
command -v shasum >/dev/null 2>&1 || { printf 'ERROR: shasum is required for input manifests\n' >&2; exit 2; }

CASE_DIR="${OUTPUT_ROOT}/${CASE_NO}"
EXECUTABLE="${CASE_DIR}/hyphal-flow"
SOURCE_HASH="$(shasum -a 256 "$SOURCE_FILE" "${SCRIPT_DIR}"/src-local/*.h | awk '{print $1}' | shasum -a 256 | awk '{print $1}')"
SOURCE_FILE_HASH="$(shasum -a 256 "$SOURCE_FILE" | awk '{print $1}')"
PARAM_HASH="$(shasum -a 256 "$PARAM_FILE" | awk '{print $1}')"
MANIFEST="${CASE_DIR}/run-manifest.txt"
if [[ $USE_MPI -eq 1 ]]; then
  RUN_MODE=mpi
  MANIFEST_CPUS="$MPI_CPUS"
else
  RUN_MODE=serial
  MANIFEST_CPUS=1
fi

printf 'Case %s: %s\n' "$CASE_NO" "$CASE_DIR"
printf 'Source: %s\n' "$SOURCE_FILE"
printf 'Parameters: %s\n' "$PARAM_FILE"
printf 'Mode: %s\n' "$([[ $USE_MPI -eq 1 ]] && printf 'MPI (%s ranks)' "$MPI_CPUS" || printf 'serial')"

if [[ $DRY_RUN -eq 1 ]]; then
  printf 'Dry run: would compile with %s and %s the case.\n' "$QCC" "$([[ $COMPILE_ONLY -eq 1 ]] && printf 'not execute' || printf 'execute')"
  exit 0
fi

if [[ -d "$CASE_DIR" && -n "$(find "$CASE_DIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  if [[ $RESUME -ne 1 ]]; then
    printf 'ERROR: case directory is non-empty; use a new CaseNo or an input-matched --resume: %s\n' "$CASE_DIR" >&2
    exit 2
  fi
  [[ -f "$MANIFEST" && -f "${CASE_DIR}/restart" ]] || {
    printf 'ERROR: --resume requires run-manifest.txt and restart\n' >&2
    exit 2
  }
  grep -Fxq "source_sha256=${SOURCE_HASH}" "$MANIFEST" &&
    grep -Fxq "source_file_sha256=${SOURCE_FILE_HASH}" "$MANIFEST" &&
    grep -Fxq "params_sha256=${PARAM_HASH}" "$MANIFEST" &&
    grep -Fxq "mode=${RUN_MODE}" "$MANIFEST" &&
    grep -Fxq "cpus=${MANIFEST_CPUS}" "$MANIFEST" &&
    grep -Fxq "basilisk_api=${BASILISK_API}" "$MANIFEST" &&
    [[ "$(shasum -a 256 "${CASE_DIR}/${EXEC_CODE}" | awk '{print $1}')" == "$SOURCE_FILE_HASH" ]] &&
    [[ "$(shasum -a 256 "${CASE_DIR}/case.params" | awk '{print $1}')" == "$PARAM_HASH" ]] || {
      printf 'ERROR: --resume inputs differ from the recorded case manifest\n' >&2
      exit 2
    }
else
  mkdir -p "$CASE_DIR"
  cp "$PARAM_FILE" "${CASE_DIR}/case.params"
  cp "$SOURCE_FILE" "${CASE_DIR}/${EXEC_CODE}"
  {
    printf 'source_sha256=%s\n' "$SOURCE_HASH"
    printf 'source_file_sha256=%s\n' "$SOURCE_FILE_HASH"
    printf 'params_sha256=%s\n' "$PARAM_HASH"
    printf 'source=%s\n' "$EXEC_CODE"
    printf 'case=%s\n' "$CASE_NO"
    printf 'mode=%s\n' "$RUN_MODE"
    printf 'cpus=%s\n' "$MANIFEST_CPUS"
    printf 'basilisk_api=%s\n' "$BASILISK_API"
  } > "$MANIFEST"
fi

LOG_LINES_BEFORE=0
LAST_TIME_BEFORE=-1
if [[ $RESUME -eq 1 && -f "${CASE_DIR}/log" ]]; then
  LOG_LINES_BEFORE="$(wc -l < "${CASE_DIR}/log")"
  LAST_TIME_BEFORE="$(awk 'NR > 1 { value=$3 } END { print value + 0 }' "${CASE_DIR}/log")"
fi

cd "$CASE_DIR"
compile=("$QCC" "-I${SCRIPT_DIR}/src-local" -Wall -O2 -disable-dimensions
         "$EXEC_CODE" -o hyphal-flow -lm)
if [[ "$BASILISK_API" == legacy ]]; then
  compile+=( -DHYPHAL_LEGACY_BASILISK=1 )
fi
printf 'Compiling %s\n' "$EXEC_CODE"
if [[ $USE_MPI -eq 1 ]]; then
  CC99='mpicc -std=c99 -D_GNU_SOURCE=1' "${compile[@]:0:3}" -D_MPI=1 "${compile[@]:3}"
else
  "${compile[@]}"
fi

if [[ $COMPILE_ONLY -eq 1 ]]; then
  printf 'Compile smoke passed: %s\n' "$EXECUTABLE"
  exit 0
fi

if [[ $USE_MPI -eq 1 ]]; then
  run_command=(mpirun -np "$MPI_CPUS" ./hyphal-flow case.params)
else
  run_command=(./hyphal-flow case.params)
fi
printf 'Running case %s\n' "$CASE_NO"
"${run_command[@]}"

[[ -s log && -s restart && -s final ]] || {
  printf 'ERROR: run exited without non-empty log, restart and final outputs\n' >&2
  exit 1
}
if [[ $RESUME -eq 1 ]]; then
  LOG_LINES_AFTER="$(wc -l < log)"
  LAST_TIME_AFTER="$(awk 'NR > 1 { value=$3 } END { print value + 0 }' log)"
  if ((LOG_LINES_AFTER <= LOG_LINES_BEFORE)) ||
     ! awk -v before="$LAST_TIME_BEFORE" -v after="$LAST_TIME_AFTER" \
       'BEGIN { exit(after > before ? 0 : 1) }'; then
    printf 'ERROR: resumed simulation produced no forward log progress\n' >&2
    exit 1
  fi
fi
awk 'NR > 1 && $3 > 0 { found=1 } END { exit(found ? 0 : 1) }' log || {
  printf 'ERROR: log contains no positive simulation time\n' >&2
  exit 1
}
if grep -Ein 'nan|inf|ERROR: non-' log >/dev/null 2>&1; then
  printf 'ERROR: numerical failure signature found in log\n' >&2
  exit 1
fi

printf 'Simulation smoke passed: %s\n' "$CASE_DIR"
