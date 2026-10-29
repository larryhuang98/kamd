#!/bin/bash
# KAMD smoke test for pmemd.cuda. Usage: AMBERHOME=/path/to/amber26_kamd ./run_test.sh [SPFP|DPFP]
set -euo pipefail
cd "$(dirname "$0")"
PREC="${1:-SPFP}"
PMEMD="${AMBERHOME:?set AMBERHOME to the KAMD-patched Amber 26 install}/bin/pmemd.cuda_${PREC}"
"$PMEMD" -O -i md_test.in -p prmtop -c restrt -o md_test.out -r md_test.rst
n=$(grep -c "KAMD: velocity reset, group 1," md_test.out || true)
first=$(grep -m1 "KAMD: velocity reset, group 1," md_test.out | awk '{print $NF}')
echo "KAMD groups parsed:   $(grep -c 'KAMD refer_mask' md_test.out)  (expected 4)"
echo "KAMD resets in 3000 steps: ${n}  (first at step ${first:-none}; typically 1-2 resets)"
[ "$n" -ge 1 ] && echo "PASS" || { echo "FAIL: no KAMD reset in 3000 steps"; exit 1; }
