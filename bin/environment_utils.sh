#!/bin/bash

# Environment detection and activation utilities for BenchmarkDA
# This script provides functions to detect and activate the correct environment
# regardless of the user's setup

# Function to detect available package managers
detect_package_manager() {
    # Prefer user-configured executables to avoid IT placeholder scripts
    if [ -n "${MAMBA_EXE}" ] && [ -x "${MAMBA_EXE}" ]; then
        echo "mamba_exe"
    elif [ -n "${CONDA_EXE}" ] && [ -x "${CONDA_EXE}" ]; then
        echo "conda_exe"
    elif command -v micromamba &> /dev/null; then
        echo "micromamba"
    else
        echo ""
    fi
}

# Function to detect benchmarkda environment
detect_benchmarkda_environment() {
    local pkg_manager=$(detect_package_manager)

    if [ -z "$pkg_manager" ]; then
        echo ""
        return 1
    fi

    # Get the proper command based on package manager
    local list_cmd=""
    case "$pkg_manager" in
        "mamba_exe")
            list_cmd="${MAMBA_EXE} env list"
            ;;
        "conda_exe")
            list_cmd="${CONDA_EXE} env list"
            ;;
        "micromamba")
            list_cmd="micromamba env list"
            ;;
        *)
            echo ""
            return 1
            ;;
    esac

    # Check for standard benchmarkda environment
    if eval "$list_cmd" 2>/dev/null | grep -q "benchmarkda"; then
        echo "benchmarkda"
        return 0
    fi

    # Check for common alternative names
    for env_name in "diffabundance" "kompot_v1" "da_benchmark" "benchmark_da"; do
        if eval "$list_cmd" 2>/dev/null | grep -q "$env_name"; then
            echo "$env_name"
            return 0
        fi
    done

    echo ""
    return 1
}

# Function to activate the benchmarkda environment using shell hook
activate_benchmarkda_environment() {
    local pkg_manager=$(detect_package_manager)
    local env_name=$(detect_benchmarkda_environment)

    if [ -z "$pkg_manager" ]; then
        echo "ERROR: No conda-compatible package manager found (MAMBA_EXE/CONDA_EXE or micromamba)" >&2
        echo "Please set MAMBA_EXE or CONDA_EXE environment variable, or install micromamba" >&2
        echo "Then create the benchmarkda environment using: bash setup_environment.sh" >&2
        return 1
    fi

    if [ -z "$env_name" ]; then
        echo "ERROR: BenchmarkDA environment not found" >&2
        echo "Please create it using: bash setup_environment.sh" >&2
        echo "Or manually activate your DA environment before running scripts" >&2
        return 1
    fi

    echo "Detected package manager: $pkg_manager"
    echo "Detected environment: $env_name"

    # First, completely deactivate any existing mamba/conda environments
    # This ensures module paths take precedence when we activate
    echo "Deactivating any existing conda/mamba environments..."
    case "$pkg_manager" in
        "mamba_exe")
            eval "$(${MAMBA_EXE} shell hook --shell bash)" 2>/dev/null
            while [ ! -z "$CONDA_PREFIX" ]; do
                mamba deactivate 2>/dev/null || break
            done
            ;;
        "conda_exe")
            eval "$(${CONDA_EXE} shell hook --shell bash)" 2>/dev/null
            while [ ! -z "$CONDA_PREFIX" ]; do
                conda deactivate 2>/dev/null || break
            done
            ;;
        "micromamba")
            eval "$(micromamba shell hook --shell bash)" 2>/dev/null
            while [ ! -z "$CONDA_PREFIX" ]; do
                micromamba deactivate 2>/dev/null || break
            done
            ;;
    esac

    # Now activate the benchmarkda environment
    # Module paths (if loaded) will be at the front of PATH
    case "$pkg_manager" in
        "mamba_exe")
            mamba activate "$env_name" 2>/dev/null
            ;;
        "conda_exe")
            conda activate "$env_name" 2>/dev/null
            ;;
        "micromamba")
            micromamba activate "$env_name" 2>/dev/null
            ;;
    esac

    # Verify activation succeeded
    if [ -z "$CONDA_PREFIX" ]; then
        echo "ERROR: Failed to activate environment $env_name" >&2
        return 1
    fi

    echo "Successfully activated environment: $env_name"
    echo "CONDA_PREFIX: $CONDA_PREFIX"
    return 0
}

# Function to get the proper run command for executing in environment
get_environment_run_command() {
    local pkg_manager=$(detect_package_manager)
    local env_name=$(detect_benchmarkda_environment)

    if [ -z "$pkg_manager" ] || [ -z "$env_name" ]; then
        echo ""
        return 1
    fi

    case "$pkg_manager" in
        "mamba_exe")
            echo "${MAMBA_EXE} run -n $env_name"
            ;;
        "conda_exe")
            echo "${CONDA_EXE} run -n $env_name"
            ;;
        "micromamba")
            echo "micromamba run -n $env_name"
            ;;
        *)
            echo ""
            return 1
            ;;
    esac
}

# Function to get environment activation command for script generation
get_environment_activation_command() {
    local pkg_manager=$(detect_package_manager)
    local env_name=$(detect_benchmarkda_environment)

    if [ -z "$pkg_manager" ] || [ -z "$env_name" ]; then
        echo "# Environment detection failed, assuming correct environment is active"
        return
    fi

    cat << EOF
# Detect and activate benchmarkda environment
if command -v $pkg_manager &> /dev/null; then
    if $pkg_manager env list | grep -q "$env_name"; then
        $pkg_manager deactivate 2>/dev/null || true
        $pkg_manager activate $env_name 2>/dev/null || echo "Using existing environment"
    else
        echo "WARNING: Environment $env_name not found, using current environment"
    fi
else
    echo "WARNING: $pkg_manager not found, using current environment"
fi
EOF
}

# Note: Functions are defined and can be sourced by other scripts
# Shell compatibility: export -f works in bash but not zsh
# These functions will be available after sourcing this file