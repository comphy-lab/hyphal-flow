#!/bin/bash
# parse_params.sh
#
# Shell helper library for hyphal-flow parameter files and sweeps.
# Source this file in scripts that need parameter parsing utilities.
#
# Usage:
#   source src-local/parse_params.sh
#
# Provides:
# - parse_param_file <file>
# - get_param <key> [default]
# - set_param_in_file <key> <value> <file>
# - generate_sweep_cases <sweep_file>
# - validate_required_params <key> [key...]
# - print_params

trim_string() {
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf '%s' "$s"
}

clear_loaded_params() {
  local var
  for var in $(compgen -v PARAM_); do
    unset "$var"
  done
}

# Parse key=value parameters and export as PARAM_<key> environment variables.
# Usage: parse_param_file <file>
parse_param_file() {
  local param_file="$1"
  local line key value

  if [[ ! -f "$param_file" ]]; then
    echo "ERROR: Parameter file not found: $param_file" >&2
    return 1
  fi

  clear_loaded_params

  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%%#*}"
    line="$(trim_string "$line")"
    [[ -z "$line" ]] && continue
    [[ "$line" != *=* ]] && continue

    key="${line%%=*}"
    value="${line#*=}"
    key="$(trim_string "$key")"
    value="$(trim_string "$value")"

    [[ -z "$key" ]] && continue
    [[ -z "$value" ]] && continue

    export "PARAM_${key}=${value}"
  done < "$param_file"

  return 0
}

# Get loaded parameter value with optional default.
# Usage: get_param <key> [default]
get_param() {
  local key="$1"
  local default="${2:-}"
  local var_name="PARAM_${key}"
  printf '%s\n' "${!var_name:-$default}"
}

# Set/update key=value inside a parameter file.
# Usage: set_param_in_file <key> <value> <file>
set_param_in_file() {
  local key="$1"
  local value="$2"
  local file="$3"

  if [[ ! -f "$file" ]]; then
    echo "ERROR: Parameter file not found: $file" >&2
    return 1
  fi

  if grep -q "^${key}=" "$file"; then
    sed -i'.bak' "s|^${key}=.*|${key}=${value}|" "$file"
  else
    printf '%s=%s\n' "$key" "$value" >> "$file"
  fi
  rm -f "${file}.bak"
  return 0
}

# Generate sweep case files from a sweep config.
# Usage: generate_sweep_cases <sweep_file>
# Output: prints temp directory path containing case_*.params and cases.list
generate_sweep_cases() {
  local sweep_file="$1"
  local config_dir base_config case_start case_end expected_count
  local line var_name var_values
  local -a sweep_vars=()
  local -a sweep_values=()
  local -a generated_files=()
  local temp_dir case_num combination_count

  if [[ ! -f "$sweep_file" ]]; then
    echo "ERROR: Sweep file not found: $sweep_file" >&2
    return 1
  fi

  # shellcheck disable=SC1090
  source "$sweep_file"

  config_dir="$(cd "$(dirname "$sweep_file")" && pwd)"
  base_config="${BASE_CONFIG:-default.params}"
  case_start="${CASE_START:-1000}"
  case_end="${CASE_END:-}"

  if [[ "$base_config" != /* ]]; then
    base_config="${config_dir}/${base_config}"
  fi

  if [[ ! -f "$base_config" ]]; then
    echo "ERROR: Base config not found: $base_config" >&2
    return 1
  fi

  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%%#*}"
    line="$(trim_string "$line")"
    [[ -z "$line" ]] && continue

    if [[ "$line" =~ ^SWEEP_([^=]+)=(.*)$ ]]; then
      var_name="${BASH_REMATCH[1]}"
      var_values="$(trim_string "${BASH_REMATCH[2]}")"
      if [[ -z "$var_values" ]]; then
        echo "ERROR: Empty sweep values for SWEEP_${var_name}" >&2
        return 1
      fi
      sweep_vars+=("$var_name")
      sweep_values+=("$var_values")
    fi
  done < "$sweep_file"

  if [[ ${#sweep_vars[@]} -eq 0 ]]; then
    echo "ERROR: No SWEEP_* variables found in $sweep_file" >&2
    return 1
  fi

  temp_dir="$(mktemp -d "${TMPDIR:-/tmp}/hyphal-sweep.XXXXXX")"
  case_num="$case_start"
  combination_count=0

  _generate_recursive() {
    local depth="$1"
    shift || true
    local -a current_values=("$@")
    local case_file i raw_val trimmed_val
    local -a value_list=()

    if [[ "$depth" -eq "${#sweep_vars[@]}" ]]; then
      case_file="${temp_dir}/case_$(printf '%04d' "$case_num").params"
      cp "$base_config" "$case_file"
      set_param_in_file "CaseNo" "$case_num" "$case_file"

      for i in "${!sweep_vars[@]}"; do
        set_param_in_file "${sweep_vars[$i]}" "${current_values[$i]}" "$case_file"
      done

      generated_files+=("$case_file")
      ((case_num += 1))
      ((combination_count += 1))
      return
    fi

    IFS=',' read -r -a value_list <<< "${sweep_values[$depth]}"
    for raw_val in "${value_list[@]}"; do
      trimmed_val="$(trim_string "$raw_val")"
      [[ -z "$trimmed_val" ]] && continue
      if [[ ${#current_values[@]} -gt 0 ]]; then
        _generate_recursive "$((depth + 1))" "${current_values[@]}" "$trimmed_val"
      else
        _generate_recursive "$((depth + 1))" "$trimmed_val"
      fi
    done
  }

  _generate_recursive 0

  if [[ -n "$case_end" ]]; then
    expected_count=$((case_end - case_start + 1))
    if [[ "$expected_count" -ne "$combination_count" ]]; then
      echo "ERROR: CASE_START/CASE_END imply ${expected_count} cases, generated ${combination_count}" >&2
      rm -rf "$temp_dir"
      return 1
    fi
  fi

  printf '%s\n' "${generated_files[@]}" > "${temp_dir}/cases.list"
  printf '%s\n' "$temp_dir"
  return 0
}

# Validate that required parameters exist in the loaded PARAM_* set.
# Usage: validate_required_params <key> [key...]
validate_required_params() {
  local missing=0
  local key var_name

  for key in "$@"; do
    var_name="PARAM_${key}"
    if [[ -z "${!var_name:-}" ]]; then
      echo "ERROR: Required parameter '${key}' not found" >&2
      missing=1
    fi
  done

  return "$missing"
}

# Print currently loaded PARAM_* variables.
print_params() {
  local var key
  echo "Loaded parameters:"
  while IFS= read -r var; do
    key="${var#PARAM_}"
    printf '  %s = %s\n' "$key" "${!var}"
  done < <(compgen -v PARAM_ | sort)
}
