#!/bin/bash
#SBATCH -J kamd_omm_val
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 03:00:00
#SBATCH -o slurm_%j.log
# Site-specific settings come from environment variables:
#   PMEMD        KAMD-patched pmemd.cuda_SPFP          AMBER_TEST  amber/test of this repository
#   KAMD_OPENMM  openmm/ of this repository            PY          python with openmm, parmed, scipy
#   KAMD_GMX_DIR directory holding install/ (KAMD-patched GROMACS) and validation/
# Load your compiler/CUDA/Amber runtime libraries (LD_LIBRARY_PATH) before submitting.
# Validation of kamd_openmm against pmemd.cuda KAMD on the 29,467-atom test system.
set -u
for i in $(seq 30); do [ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ] && break; sleep 10; done
echo "host $(hostname)  $(nvidia-smi --query-gpu=name --format=csv,noheader)  $(date)"

ROOT=${KAMD_OPENMM:?set KAMD_OPENMM to the openmm/ folder of this repository}
SRC=${AMBER_TEST:?set AMBER_TEST to amber/test of this repository}
PMEMD=${PMEMD:?set PMEMD to the KAMD-patched pmemd.cuda_SPFP}
PY=${PY:-python3}
export PYTHONPATH=$ROOT
W=$ROOT/validation/amber29k; mkdir -p $W; cd $W
ln -sf $SRC/prmtop prmtop; ln -sf $SRC/restrt restrt

mdin() {  # $1 nstlim  $2 ntwr  $3 ntxo
cat <<EOF
KAMD validation (same settings as test/kamd/Run_test.sh)
&cntrl
   imin=0, irest=1, ntx=5,
   ntpr=500, ntwx=0, ntwr=$2, ntxo=$3,
   nstlim=$1, iwrap=1,
   dt=0.002, ntt=3, gamma_ln=5,
   temp0=355, ig=1060,
   ntp=1, ntc=2, ntf=2,
   ntb=2, taup=2.0, cut=8.0,
/
&kamd
   pseudobd_cycles = 250,
   shell_cut       = 5.0,
   surface_cut     = 3.0,
   nzero_crossings = 2,
   refer_mask(1) = "@20-96",
   refer_mask(2) = "@285-339",
   refer_mask(3) = "@366-482",
   refer_mask(4) = "@568-611",
   solvent_mask  = "@623-29467",
/
EOF
}

echo "=== 1. pmemd.cuda KAMD_DIAG run (stops at the first reset)"
mdin 20000 20000 1 > diagA.in
KAMD_DIAG=1 $PMEMD -O -i diagA.in -p prmtop -c restrt -o diagA.out -r diagA.rst
S=$(grep "first-reset diagnostic at step" diagA.out | awk '{print $NF}')
echo "first reset at step $S"; grep "KAMD-DIAG" diagA.out

echo "=== 2. same run stopped at completed step S-1, ASCII restart"
mdin $((S-1)) $((S-1)) 1 > stateB.in
$PMEMD -O -i stateB.in -p prmtop -c restrt -o stateB.out -r stateB.rst7

echo "=== 3. first-shell parity"
$PY $ROOT/validation/parity_first_shell.py --prmtop prmtop --mdin stateB.in --restart stateB.rst7 \
    --gpu-diag kamd_diag_gpu.txt --json parity_first_shell.json

NST=${NST:-250000}
echo "=== 4. pmemd.cuda KAMD, $NST steps"
mdin $NST $NST 1 > amber_long.in
/usr/bin/time -f "amber_wall %e" $PMEMD -O -i amber_long.in -p prmtop -c restrt -o amber_long.out -r amber_long.rst7
grep "velocity reset, group 1," amber_long.out | awk '{print $NF}' > amber_long.resets
echo "amber resets: $(wc -l < amber_long.resets)"; grep "ns/day" amber_long.out | tail -1

echo "=== 5. OpenMM + KAMD, $NST steps"
$PY $ROOT/validation/run_openmm_kamd.py --prmtop prmtop --restart restrt --mdin amber_long.in \
    --nsteps $NST --out omm_kamd
echo "=== 6. OpenMM without KAMD, $NST steps (overhead baseline)"
$PY $ROOT/validation/run_openmm_kamd.py --prmtop prmtop --restart restrt --mdin amber_long.in \
    --nsteps $NST --no-kamd --out omm_plain
echo "done $(date)"
