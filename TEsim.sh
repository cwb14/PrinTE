#!/bin/bash

# TEsim.sh: Wrapper for TE simulation tools

TEsim_version="1.00"

# Dynamically determine the directory of the current script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

# Paths to the tools relative to the repository's structure
RANDSEQGEN_DIR="$SCRIPT_DIR/RandSeqGen"
LTR_SIMULATOR_DIR="$SCRIPT_DIR/LTR_simulator"
TEGENOMESIMULATOR_DIR="$SCRIPT_DIR/TEgenomeSimulator"

print_main_help() {
    cat << EOF
TEsim (Tools for simulation Transposable element evolution)
Version: $TEsim_version

Usage:   TEsim <command> [options]

Commands:

     RandSeqGen             TIR insertion
     LTR_simulator          LTR mutation, insertion, and evolution
     TEgenomeSimulator      TE mutation, insertion, and tracking

Use "TEsim <command> -h" for more information about a command.
EOF
}

if [ $# -lt 1 ]; then
    print_main_help
    exit 1
fi

COMMAND="$1"
shift

case "$COMMAND" in
    RandSeqGen)
        python "$RANDSEQGEN_DIR/RandSeqGen.py" "$@"
        ;;
    LTR_simulator)
        bash "$LTR_SIMULATOR_DIR/LTR_simulator.sh" "$@"
        ;;
    TEgenomeSimulator)
        python "$TEGENOMESIMULATOR_DIR/TEgenomeSimulator.py" "$@"
        ;;
    -h|--help)
        print_main_help
        ;;
    *)
        echo "Unknown command: $COMMAND" >&2
        print_main_help
        exit 1
        ;;
esac
