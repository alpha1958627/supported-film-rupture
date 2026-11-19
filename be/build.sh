#!/bin/bash
set -euo pipefail
: "${BASILISK:?Set BASILISK to the Basilisk src directory}"
export PATH="$BASILISK:$PATH"
mkdir -p build_tmp
export TMPDIR="$PWD/build_tmp"
qcc -disable-dimensions -D_GNU_SOURCE -Wall -O2 -g -fopenmp -fPIC -shared wall_native.c -o wall_native.so -lm
