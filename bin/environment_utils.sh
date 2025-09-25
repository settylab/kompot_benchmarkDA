#!/bin/bash

# Environment detection and activation utilities for BenchmarkDA
# This script provides functions to detect and activate the correct environment
# regardless of the user's setup

# Function to detect available package managers
detect_package_manager() {
    if command -v mamba &> /dev/null; then
        echo "mamba"
    elif command -v micromamba &> /dev/null; then
        echo "micromamba"
    elif command -v conda &> /dev/null; then
        echo "conda"
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

    # Check for standard benchmarkda environment
    if $pkg_manager env list | grep -q "benchmarkda"; then
        echo "benchmarkda"
        return 0
    fi

    # Check for common alternative names
    for env_name in "diffabundance" "kompot_v1" "da_benchmark" "benchmark_da"; do
        if $pkg_manager env list | grep -q "$env_name"; then
            echo "$env_name"
            return 0
        fi
    done

    echo ""
    return 1
}

# Function to activate the benchmarkda environment
activate_benchmarkda_environment() {
    local pkg_manager=$(detect_package_manager)
    local env_name=$(detect_benchmarkda_environment)

    if [ -z "$pkg_manager" ]; then
        echo "ERROR: No conda-compatible package manager found (conda/mamba/micromamba)" >&2
        echo "Please install one and create the benchmarkda environment using:" >&2
        echo "  bash setup_environment.sh" >&2
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

    # Use MAMBA_EXE if available for more reliable execution
    if [ -n "$MAMBA_EXE" ]; then
        echo "Using MAMBA_EXE for environment execution"
        return 0
    fi

    # Fallback: Set up shell hook and activate environment
    case "$pkg_manager" in
        "micromamba")
            eval "$(micromamba shell hook --shell bash 2>/dev/null)" || true
            micromamba deactivate 2>/dev/null || true
            micromamba activate "$env_name" 2>/dev/null || {
                echo "WARNING: Failed to activate $env_name, using MAMBA_EXE fallback" >&2
                return 0
            }
            ;;
        "mamba")
            eval "$(conda shell.bash hook 2>/dev/null)" || true
            mamba activate "$env_name" 2>/dev/null || {
                echo "WARNING: Failed to activate $env_name, using MAMBA_EXE fallback" >&2
                return 0
            }
            ;;
        "conda")
            eval "$(conda shell.bash hook 2>/dev/null)" || true
            conda activate "$env_name" 2>/dev/null || {
                echo "WARNING: Failed to activate $env_name, using MAMBA_EXE fallback" >&2
                return 0
            }
            ;;
    esac

    echo "Successfully activated environment: $env_name"
    return 0
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