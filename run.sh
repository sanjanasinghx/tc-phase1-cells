#!/usr/bin/env bash
# compile the channels (first run only), then run the phase 1 gate.
# exits non-zero if any check fails.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
if [ ! -d arm64 ] && [ ! -d x86_64 ]; then nrnivmodl mod > nrnivmodl.log 2>&1; fi
mkdir -p results
$PY phase1_gate.py | tee results/phase1_gate.log
