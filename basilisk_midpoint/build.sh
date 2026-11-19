#!/bin/bash
set -euo pipefail
: "${BASILISK:?Set BASILISK to Basilisk src}"
: "${CC99:?Set CC99 to the MPI C compiler}"
export PATH="$BASILISK:$PATH"
mkdir -p compiled build_tmp
export TMPDIR="$PWD/build_tmp"
qcc -I"$BASILISK" -D_GNU_SOURCE -D_MPI=1 -O3 -fno-fast-math -ffp-contract=off -shared -fPIC basilisk_midpoint.c -o compiled/basilisk_midpoint.so -lm
sha256sum basilisk_midpoint.c compiled/basilisk_midpoint.so
