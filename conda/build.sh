#!/bin/bash
set -euo pipefail

# The mutator is compiled here rather than by PrinTE.sh at first run, so an installed
# package needs no compiler and no writable install prefix.
mkdir -p "${PREFIX}/bin"
${CXX} -std=c++17 -O3 -fopenmp src/printe/cpp/ltr_mutator.cpp -o "${PREFIX}/bin/ltr_mutator"

${PYTHON} -m pip install . -vv --no-deps --no-build-isolation
