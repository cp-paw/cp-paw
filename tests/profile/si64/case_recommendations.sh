#!/usr/bin/env bash

cppaw_load_capabilities() {
  local helper

  if [[ "${CPPAW_GPU_CAPABILITIES_LOADED:-no}" == yes ]]; then
    [[ -n "${CPPAW_GPU_CAPABILITIES_OUTPUT:-}" ]]
    return $?
  fi

  CPPAW_GPU_CAPABILITIES_OUTPUT=""
  if [[ -n "${CPPAW_GPU_CAPABILITIES_FILE:-}" ]]; then
    if [[ -f "${CPPAW_GPU_CAPABILITIES_FILE}" ]]; then
      CPPAW_GPU_CAPABILITIES_OUTPUT=$(cat "${CPPAW_GPU_CAPABILITIES_FILE}")
    fi
  else
    helper=${CPPAW_GPU_CAPABILITIES_CMD:-"${ROOT}/src/Tools/Scripts/paw_gpu_capabilities.sh"}
    if [[ -x "${helper}" ]]; then
      CPPAW_GPU_CAPABILITIES_OUTPUT=$("${helper}" 2>/dev/null || true)
    fi
  fi
  CPPAW_GPU_CAPABILITIES_LOADED=yes

  [[ -n "${CPPAW_GPU_CAPABILITIES_OUTPUT}" ]]
}

cppaw_capability_value() {
  local key=$1
  cppaw_load_capabilities || return 1
  printf '%s\n' "${CPPAW_GPU_CAPABILITIES_OUTPUT}" \
    | awk -F= -v key="${key}" '$1 == key { sub(/^[^=]*=/, ""); value=$0 } END { print value }'
}

cppaw_dedup_space_list() {
  local item
  local norm=""
  for item in "$@"; do
    [[ -n "${item}" ]] || continue
    case " ${norm} " in
      *" ${item} "*) ;;
      *) norm="${norm:+${norm} }${item}" ;;
    esac
  done
  echo "${norm}"
}

cppaw_append_case() {
  local list=$1
  local item=$2
  cppaw_dedup_space_list ${list} "${item}"
}

cppaw_normalize_capability_cases() {
  local cases=$1
  if [[ "${cases}" == none ]]; then
    echo ""
  else
    echo "${cases}"
  fi
}

cppaw_resolve_recommended_cases() {
  local selected=$1
  local key=$2
  local fallback=${3:-}
  local value

  case "${selected}" in
    auto|recommended)
      if value=$(cppaw_capability_value "${key}"); then
        if [[ -n "${value}" ]]; then
          cppaw_normalize_capability_cases "${value}"
        else
          echo "${fallback}"
        fi
      else
        echo "${fallback}"
      fi
      ;;
    *)
      echo "${selected}"
      ;;
  esac
}

cppaw_resolve_recommended_gpu_exploration_cases() {
  local selected=$1
  local fallback=${2:-}
  local gpu_cases diagnostic_cases raw_gpu_cases raw_diagnostic_cases

  case "${selected}" in
    auto|recommended)
      if cppaw_load_capabilities; then
        raw_gpu_cases=$(cppaw_capability_value recommended_gpu_cases || true)
        raw_diagnostic_cases=$(cppaw_capability_value recommended_gpu_diagnostic_cases || true)
        if [[ -z "${raw_gpu_cases}" && -z "${raw_diagnostic_cases}" ]]; then
          echo "${fallback}"
          return 0
        fi
        gpu_cases=$(cppaw_normalize_capability_cases "${raw_gpu_cases}")
        diagnostic_cases=$(cppaw_normalize_capability_cases "${raw_diagnostic_cases}")
        cppaw_dedup_space_list ${gpu_cases} ${diagnostic_cases}
      else
        echo "${fallback}"
      fi
      ;;
    *)
      echo "${selected}"
      ;;
  esac
}

cppaw_write_capabilities_file() {
  local path=$1
  cppaw_load_capabilities || return 1
  printf '%s\n' "${CPPAW_GPU_CAPABILITIES_OUTPUT}" > "${path}"
}

cppaw_prepend_path_value() {
  local name=$1
  local value=$2
  local current

  [[ -n "${value}" && "${value}" != none ]] || return 0
  current=${!name:-}
  case ":${current}:" in
    *":${value}:"*) ;;
    *) export "${name}=${value}${current:+:${current}}" ;;
  esac
}

cppaw_apply_capability_env() {
  local path

  path=$(cppaw_capability_value host_fftw_pkg_config_path || true)
  cppaw_prepend_path_value PKG_CONFIG_PATH "${path}"

  path=$(cppaw_capability_value host_fftw_ld_library_path || true)
  cppaw_prepend_path_value LD_LIBRARY_PATH "${path}"
}
