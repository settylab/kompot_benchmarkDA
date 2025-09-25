#!/bin/bash
set -e

# BenchmarkDA: Main execution script
# Orchestrates the complete differential abundance benchmarking workflow
# This script provides a clear, readable overview of what gets executed

# =============================================================================
# CONFIGURATION
# =============================================================================

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd "$SCRIPT_DIR"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m'

print_header() {
    echo -e "${BOLD}${BLUE}========================================${NC}"
    echo -e "${BOLD}${BLUE}$1${NC}"
    echo -e "${BOLD}${BLUE}========================================${NC}"
}

print_step() {
    echo -e "${GREEN}[STEP]${NC} $1"
}

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

# =============================================================================
# MAIN BENCHMARKING WORKFLOW
# =============================================================================

print_header "BenchmarkDA: Complete Benchmarking Pipeline"

echo ""
print_info "This script runs the complete differential abundance benchmarking pipeline:"
print_info "1. Environment setup and validation"
print_info "2. Dataset preprocessing with embeddings (DM + PCA)"
print_info "3. Synthetic condition label generation"
print_info "4. All methods benchmarking (Python + R)"
echo ""

print_step "DATASETS TO PROCESS:"
echo "  • Synthetic: linear, branch, cluster"
echo "  • Real: covid19-pbmc, pancreas, bcr-xl, levine32 (if available)"
echo ""

print_step "METHODS TO RUN:"
echo "  • Python: mellon, meld, cna, kompot (+ variants)"
echo "  • R: milo, daseq, cydar, louvain"
echo ""

print_step "EMBEDDINGS:"
echo "  • Diffusion Maps (DM): Dataset-specific components"
echo "  • Principal Components (PCA): Standard dimensionality"
echo ""

# Prompt user for confirmation
echo -e "${YELLOW}This will run the complete benchmarking pipeline.${NC}"
echo "Continue? [Y/n]"
read -r response
if [[ "$response" =~ ^[Nn]$ ]]; then
    echo "Benchmarking cancelled."
    exit 0
fi

# =============================================================================
# EXECUTE COMPLETE PIPELINE
# =============================================================================

print_header "Executing Complete Benchmarking Pipeline"

print_step "Running all steps: setup → preprocess → labels → benchmark"
print_info "Command: ./cli.sh all"

# Execute the complete pipeline using our modular CLI
if ! ./cli.sh all; then
    print_error "Pipeline execution failed!"
    print_info "For debugging, try individual steps:"
    print_info "  ./cli.sh setup"
    print_info "  ./cli.sh preprocess"
    print_info "  ./cli.sh labels"
    print_info "  ./cli.sh benchmark"
    exit 1
fi

# =============================================================================
# COMPLETION SUMMARY
# =============================================================================

print_header "Benchmarking Pipeline Completed Successfully"

echo ""
print_success "WHAT WAS EXECUTED:"
echo "  ✓ Environment setup and validation"
echo "  ✓ Dataset preprocessing (DM + PCA embeddings)"
echo "  ✓ Synthetic condition label generation"
echo "  ✓ Python method benchmarks (mellon, meld, cna, kompot)"
echo "  ✓ R method benchmarks (milo, daseq, cydar, louvain)"
echo ""

print_success "RESULTS LOCATION:"
echo "  • Complete results: benchmark/"
echo "  • Individual results: results/individual/"
echo "  • Generated scripts: benchmark_scripts/"
echo "  • Execution logs: SlurmLog/"
echo ""

print_step "NEXT STEPS:"
echo "• Analyze results in the benchmark/ directory"
echo "• Check individual method results in results/individual/"
echo "• Review generated scripts in benchmark_scripts/"
if command -v sbatch &> /dev/null; then
    echo "• Monitor any running Slurm jobs: squeue -u \$USER"
fi
echo ""

print_step "INDIVIDUAL EXECUTION EXAMPLES:"
echo "# Run specific components only:"
echo "./cli.sh --datasets linear,branch preprocess"
echo "./cli.sh --methods python benchmark"
echo "./cli.sh --embeddings dm --methods python benchmark"
echo ""
echo "# Run individual methods:"
echo "python run_da_method.py kompot linear --embedding dm"
echo "python run_da_method.py --list  # Show all available methods"
echo ""

print_step "TROUBLESHOOTING:"
echo "# Test your setup:"
echo "python test_setup.py"
echo ""
echo "# Dry run to see what would execute:"
echo "./cli.sh --dry-run all"
echo ""
echo "# Get help:"
echo "./cli.sh --help"

print_header "All Done! Happy Analyzing! 🎉"