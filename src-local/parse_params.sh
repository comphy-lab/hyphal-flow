#!/usr/bin/env bash

# Shared shell-side parser for trusted key=value parameter and sweep files.
# Keys are allowlisted before becoming environment-variable names; values are
# treated as data and are never evaluated or sourced as shell code.

trim_string() {
  local value="${1:-}"
  value="${value#"${value%%[![:space:]]*}"}"
  value="${value%"${value##*[![:space:]]}"}"
  printf '%s' "$value"
}

valid_param_key() {
  [[ "${1:-}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]]
}

clear_loaded_params() {
  local variable
  while IFS= read -r variable; do
    [[ "$variable" == PARAM_* ]] && unset "$variable"
  done < <(compgen -v PARAM_ || true)
}

parse_param_file() {
  local file="${1:-}"
  local line key value
  [[ -f "$file" ]] || { printf 'ERROR: parameter file not found: %s\n' "$file" >&2; return 1; }
  clear_loaded_params

  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%%#*}"
    line="$(trim_string "$line")"
    [[ -z "$line" || "$line" != *=* ]] && continue
    key="$(trim_string "${line%%=*}")"
    value="$(trim_string "${line#*=}")"
    valid_param_key "$key" || {
      printf 'ERROR: invalid parameter key: %s\n' "$key" >&2
      return 1
    }
    [[ -n "$value" ]] || continue
    export "PARAM_${key}=${value}"
  done < "$file"
}

get_param() {
  local key="${1:-}"
  local default_value="${2:-}"
  local variable="PARAM_${key}"
  valid_param_key "$key" || return 1
  printf '%s\n' "${!variable:-$default_value}"
}

# Print the last value for a key, matching the C parser's duplicate-key rule.
get_param_value() {
  local key="${1:-}"
  local file="${2:-}"
  local default_value="${3:-}"
  valid_param_key "$key" || { printf 'ERROR: invalid parameter key: %s\n' "$key" >&2; return 1; }
  [[ -f "$file" ]] || { printf 'ERROR: parameter file not found: %s\n' "$file" >&2; return 1; }

  awk -F '=' -v wanted="$key" -v fallback="$default_value" '
    {
      line = $0
      sub(/[[:space:]]*#.*/, "", line)
      separator = index(line, "=")
      if (!separator) next
      key = substr(line, 1, separator - 1)
      value = substr(line, separator + 1)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", key)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
      if (key == wanted && value != "") {
        found = 1
        result = value
      }
    }
    END { print found ? result : fallback }
  ' "$file"
}

set_param_in_file() {
  local key="${1:-}"
  local value="${2:-}"
  local file="${3:-}"
  local temporary
  valid_param_key "$key" || { printf 'ERROR: invalid parameter key: %s\n' "$key" >&2; return 1; }
  [[ -f "$file" ]] || { printf 'ERROR: parameter file not found: %s\n' "$file" >&2; return 1; }
  [[ "$value" != *$'\n'* && "$value" != *$'\r'* ]] || {
    printf 'ERROR: parameter values may not contain newlines\n' >&2
    return 1
  }

  temporary="$(mktemp "${TMPDIR:-/tmp}/hyphal-param.XXXXXX")"
  awk -v wanted="$key" -v replacement="$value" '
    BEGIN { replaced = 0 }
    {
      line = $0
      plain = line
      sub(/[[:space:]]*#.*/, "", plain)
      separator = index(plain, "=")
      candidate = separator ? substr(plain, 1, separator - 1) : ""
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", candidate)
      if (candidate == wanted) {
        if (!replaced) print wanted "=" replacement
        replaced = 1
      }
      else print line
    }
    END { if (!replaced) print wanted "=" replacement }
  ' "$file" > "$temporary"
  mv "$temporary" "$file"
}

print_params() {
  local variable
  while IFS= read -r variable; do
    [[ "$variable" == PARAM_* ]] || continue
    printf '%s=%s\n' "${variable#PARAM_}" "${!variable}"
  done < <(compgen -v PARAM_ | sort)
}
