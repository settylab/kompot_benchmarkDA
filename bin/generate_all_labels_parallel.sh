#!/bin/bash
# Parallel label generation for all datasets
# Runs all datasets in parallel, with covid19-pbmc split into ~8 parameter combinations

set -e

export PYTHONUNBUFFERED=1 # prevent python output buffering

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

echo "============================================================"
echo "Parallel Label Generation for All Datasets"
echo "============================================================"
echo ""

# Array to store background job PIDs and descriptions
declare -a PIDS
declare -a DESCS
declare -a LOG_FILES

# Function to run label generation in background using CLI
run_labels() {
    local dataset="$1"
    local filters="$2"
    local desc="$3"
    local delay="$4"  # Delay in seconds before starting

    local log_file="logs/labels_${dataset}_${desc//[^a-zA-Z0-9]/_}.log"

    echo "Starting: $desc (delay: ${delay}s)"

    # Build filter arguments for CLI
    local filter_args=""
    if [ -n "$filters" ]; then
        filter_args="$filters"
    fi

    # Run via CLI with delay
    (sleep $delay && ./cli.sh labels --datasets "$dataset" $filter_args > "$log_file" 2>&1) &

    PIDS+=($!)
    DESCS+=("$desc")
    LOG_FILES+=("$log_file")
}

# Create logs directory
mkdir -p logs

# Start label generation for each dataset in parallel with staggered delays
echo "Launching label generation jobs with 60s delays..."
echo ""

# Synthetic datasets (small, run as single jobs)
run_labels "linear" "" "linear-all" 0
run_labels "branch" "" "branch-all" 0
run_labels "cluster" "" "cluster-all" 0

# Real datasets
run_labels "bcr-xl" "" "bcr-xl-all" 0
run_labels "levine32" "" "levine32-all" 0
run_labels "pancreas" "" "pancreas-all" 0

# COVID19-PBMC: Split into parallel jobs by population
# Actual populations: PB, CD14_Monocyte, CD8_T, CD4_T, Platelet, NK, Granulocyte, CD16_Monocyte, gd_T, pDC, DC
# Split into groups for parallel execution (11 populations -> split into manageable groups)
run_labels "covid19-pbmc" "--populations PB" "covid19-pbmc-PB" 0
run_labels "covid19-pbmc" "--populations CD14_Monocyte,CD16_Monocyte" "covid19-pbmc-Monocytes" 0
run_labels "covid19-pbmc" "--populations CD4_T,CD8_T" "covid19-pbmc-T-cells" 0
run_labels "covid19-pbmc" "--populations NK" "covid19-pbmc-NK" 0
run_labels "covid19-pbmc" "--populations DC,pDC" "covid19-pbmc-DC" 0
run_labels "covid19-pbmc" "--populations Granulocyte" "covid19-pbmc-Granulocyte" 0
run_labels "covid19-pbmc" "--populations gd_T" "covid19-pbmc-gd-T" 0
run_labels "covid19-pbmc" "--populations Platelet" "covid19-pbmc-Platelet" 0

echo ""
echo "============================================================"
echo "Started ${#PIDS[@]} parallel jobs"
echo "============================================================"
echo ""

# Display job information
for i in "${!PIDS[@]}"; do
    echo "[$((i+1))] ${DESCS[$i]} (PID: ${PIDS[$i]})"
    echo "    Log: ${LOG_FILES[$i]}"
done

echo ""
echo "============================================================"
echo "Monitoring all logs (Ctrl+C to stop monitoring)..."
echo "============================================================"
echo ""

# Single cleanup function that handles everything
cleanup() {
    local exit_code=$?
    echo ""
    echo "Cleanup: Stopping log monitoring..."

    # Kill tail process if running
    if [ -n "$TAIL_PID" ] && kill -0 $TAIL_PID 2>/dev/null; then
        kill $TAIL_PID 2>/dev/null
    fi

    # On interrupt, give summary of what was running
    if [ $exit_code -eq 130 ]; then
        echo ""
        echo "Script interrupted. Background jobs are still running."
        echo "Active PIDs: ${PIDS[*]}"
        echo "To kill all jobs: pkill -f 'cli.sh.*labels'"
    fi

    exit $exit_code
}

# Set up single trap for all exit scenarios
trap cleanup EXIT INT TERM HUP

# Tail all logs in parallel
tail -f "${LOG_FILES[@]}" 2>/dev/null &
TAIL_PID=$!

# Wait for all background jobs and track failures
echo "Waiting for all jobs to complete..."
echo "(Press Ctrl+C to stop monitoring logs - jobs will continue in background)"
echo ""

FAILED=0
SUCCESS=0

for i in "${!PIDS[@]}"; do
    PID=${PIDS[$i]}
    if wait $PID 2>/dev/null; then
        SUCCESS=$((SUCCESS + 1))
        echo "[OK] ${DESCS[$i]} completed (PID: $PID)"
    else
        EXIT_CODE=$?
        if [ $EXIT_CODE -ne 0 ] && [ $EXIT_CODE -ne 127 ]; then
            FAILED=$((FAILED + 1))
            echo "[FAIL] ${DESCS[$i]} failed with exit code $EXIT_CODE (PID: $PID)"
        fi
    fi
done

# Clean exit - disable trap and stop tail
trap - EXIT INT TERM HUP
kill $TAIL_PID 2>/dev/null || true

echo ""
echo "============================================================"
echo "All jobs completed!"
echo "  Success: $SUCCESS"
echo "  Failed: $FAILED"
echo "============================================================"
echo ""
echo "Review individual logs:"
for i in "${!LOG_FILES[@]}"; do
    echo "  ${DESCS[$i]}: ${LOG_FILES[$i]}"
done

exit $FAILED
