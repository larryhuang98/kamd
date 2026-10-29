#!/bin/bash
#SBATCH -J kamd_bench
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 00:40:00
#SBATCH -o slurm_%j.log
# Site-specific settings come from environment variables:
#   PMEMD        KAMD-patched pmemd.cuda_SPFP          AMBER_TEST  amber/test of this repository
#   KAMD_OPENMM  openmm/ of this repository            PY          python with openmm, parmed, scipy
#   KAMD_GMX_DIR directory holding install/ (KAMD-patched GROMACS) and validation/
# Load your compiler/CUDA/Amber runtime libraries (LD_LIBRARY_PATH) before submitting.
for i in $(seq 30); do [ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ] && break; sleep 10; done
echo "host $(hostname) $(nvidia-smi --query-gpu=name --format=csv,noheader)"
ROOT=${KAMD_OPENMM:?set KAMD_OPENMM to the openmm/ folder of this repository}
PMEMD=${PMEMD:?set PMEMD to the KAMD-patched pmemd.cuda_SPFP}
PY=${PY:-python3}
export PYTHONPATH=$ROOT
cd $ROOT/validation/bench
ln -sf ${AMBER_TEST:?set AMBER_TEST}/prmtop prmtop; ln -sf ${AMBER_TEST}/restrt restrt
N=50000
mk() { # $1 ntb  $2 ntp  $3 barostat
cat <<EOF
bench
&cntrl
   imin=0, irest=1, ntx=5, ntpr=5000, ntwx=0, ntwr=$N, nstlim=$N,
   dt=0.002, ntt=3, gamma_ln=5, temp0=355, ig=1060,
   ntc=2, ntf=2, cut=8.0, ntb=$1, ntp=$2, barostat=$3, taup=2.0, mcbarint=25,
/
EOF
}
mk 2 1 1 > amber_npt_berendsen.in
mk 2 1 2 > amber_npt_mc.in
mk 1 0 1 > amber_nvt.in
for c in amber_npt_berendsen amber_npt_mc amber_nvt; do
  $PMEMD -O -i $c.in -p prmtop -c restrt -o $c.out -r $c.rst7
  echo "$c $(grep 'ns/day' $c.out | tail -1)"
done
# OpenMM: NPT (MC barostat every 25 steps) and NVT, no KAMD
$PY $ROOT/validation/run_openmm_kamd.py --prmtop prmtop --restart restrt --mdin amber_npt_mc.in --nsteps $N --no-kamd --out omm_npt_mc | grep ns_per_day | sed "s/^/omm_npt_mc /"
$PY $ROOT/validation/run_openmm_kamd.py --prmtop prmtop --restart restrt --mdin amber_nvt.in --nsteps $N --no-kamd --out omm_nvt | grep ns_per_day | sed "s/^/omm_nvt /"
