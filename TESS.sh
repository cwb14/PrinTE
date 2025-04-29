#!/bin/bash

# Script base directory — assumes TESS.sh is in the root of your project
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Tool Paths
RANDSEQGEN_DIR="$SCRIPT_DIR/RandSeqGen"
LTR_SIMULATOR_DIR="$SCRIPT_DIR/LTR_simulator"
TEGENOMESIMULATOR_DIR="$SCRIPT_DIR/TEgenomeSimulator"

# Global Flags
VERBOSE=0
FORCE_UPDATE=0

# Helper Functions
log() {
    if [ "$VERBOSE" -eq 1 ]; then
        echo "INFO: $1" >&2
    fi
}

error() {
    echo "ERROR: $1" >&2
    exit 1
}

usage() {
    cat << EOF
TEsim (Tools for simulation of Transposable Element evolution)
Version: 1.00

Usage:   $0 <command> [options]

Commands:

     RandSeqGen             Simulate TIR insertion
     LTR_simulator          Simulate LTR mutation, insertion, and evolution
     TEgenomeSimulator      Simulate TE mutation, insertion, and tracking

Use "$0 <command> -h" for information about individual tools.
EOF
    exit 1
}

clone_TEgenomeSimulator() {
    local target_dir="$1"

    log "Cloning TEgenomeSimulator into $target_dir"

    if [ -d "$target_dir" ]; then
        if [ "$FORCE_UPDATE" -eq 1 ]; then
            log "Removing old TEgenomeSimulator directory due to --force-update"
            rm -rf "$target_dir"
        else
            log "Using cached TEgenomeSimulator at $target_dir"
            return 0
        fi
    fi

    git clone --depth=1 https://github.com/Plant-Food-Research-Open/TEgenomeSimulator "$target_dir"
    if [ $? -ne 0 ]; then
        error "Failed to clone TEgenomeSimulator repository."
    else
        log "Successfully cloned TEgenomeSimulator into $target_dir"
    fi
}

run_TEgenomeSimulator() {
    local args=("$@")
    local PYTHON_EXEC=$(command -v python3 || command -v python || error "Python is required but not installed.")

    # Parse flags
    FORCE_UPDATE=0
    VERBOSE=0
    while [[ "$#" -gt 0 ]]; do
        case "$1" in
            --force-update)
                FORCE_UPDATE=1
                ;;
            --verbose)
                VERBOSE=1
                ;;
            *)
                break
                ;;
        esac
        shift
    done

    local TE_SIMULATOR_DIR="${SCRIPT_DIR}/TEgenomeSimulator"
    clone_TEgenomeSimulator "$TE_SIMULATOR_DIR"

    # Decompress input FASTA if needed
    local INPUT_GZ="$TE_SIMULATOR_DIR/test/input/combined_curated_TE_lib_ATOSZM_selected.fasta.gz"
    local INPUT_FASTA="$TE_SIMULATOR_DIR/test/input/combined_curated_TE_lib_ATOSZM_selected.fasta"
    if [ -f "$INPUT_GZ" ] && [ ! -f "$INPUT_FASTA" ]; then
        log "Decompressing $INPUT_GZ"
        gunzip "$INPUT_GZ" || error "Failed to gunzip $INPUT_GZ"
    elif [ ! -f "$INPUT_FASTA" ]; then
        error "Input FASTA missing: $INPUT_FASTA"
    fi

    # Navigate to root of cloned repo
    cd "$TE_SIMULATOR_DIR" || error "Failed to enter $TE_SIMULATOR_DIR"

    # Convert all file paths to absolute paths
    local cleaned_args=()
    for arg in "$@"; do
        if [[ "$arg" == *TEgenomeSimulator/* ]]; then
            abs_path="${arg#*TEgenomeSimulator/}"
            abs_path="$TE_SIMULATOR_DIR/$abs_path"
            echo "DEBUG: Resolved path: $abs_path" >&2
            cleaned_args+=("$abs_path")
        else
            cleaned_args+=("$arg")
        fi
    done

    # Reconstruct the full command for clarity
    local CMD_ARGS=("${cleaned_args[@]}")
    echo "DEBUG: Running in $(pwd)" >&2
    echo "DEBUG: Full command: $PYTHON_EXEC TEgenomeSimulator/TEgenomeSimulator.py ${CMD_ARGS[*]}" >&2

    # Execute the main simulator
    $PYTHON_EXEC TEgenomeSimulator/TEgenomeSimulator.py "${CMD_ARGS[@]}"
}

# Main Execution
if [ $# -lt 1 ]; then
    usage
fi

COMMAND="$1"
shift

case "$COMMAND" in
    RandSeqGen)
        if [ ! -d "$RANDSEQGEN_DIR" ]; then
            error "RandSeqGen directory not found: $RANDSEQGEN_DIR"
        fi
        python "$RANDSEQGEN_DIR/RandSeqGen.py" "$@"
        ;;
    LTR_simulator)
        if [ ! -d "$LTR_SIMULATOR_DIR" ]; then
            error "LTR_simulator directory not found: $LTR_SIMULATOR_DIR"
        fi
        bash "$LTR_SIMULATOR_DIR/LTR_simulator.sh" "$@"
        ;;
    TEgenomeSimulator)
        run_TEgenomeSimulator "$@"
        ;;
    -h|--help|help)
        usage
        ;;
    *)
        echo "Unknown command: $COMMAND"
        usage
        ;;
esac
