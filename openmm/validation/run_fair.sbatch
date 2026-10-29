#!/bin/bash
#SBATCH -J kamd_fair
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 01:00:00
#SBATCH -o slurm_%j.log
# Site-specific settings come from environment variables:
#   PMEMD        KAMD-patched pmemd.cuda_SPFP          AMBER_TEST  amber/test of this repository
#   KAMD_OPENMM  openmm/ of this repository            PY          python with openmm, parmed, scipy
#   KAMD_GMX_DIR directory holding install/ (KAMD-patched GROMACS) and validation/
# Load your compiler/CUDA/Amber runtime libraries (LD_LIBRARY_PATH) before submitting.
# Amber vs OpenMM, KAMD on/off, matched conditions:
#   MC barostat every 25 steps, Amber's Ewald coefficient + 64^3 grid in OpenMM, COM removal every 1000 steps.
for i in $(seq 30); do [ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ] && break; sleep 10; done
echo "host $(hostname) $(nvidia-smi --query-gpu=name --format=csv,noheader)"
ROOT=${KAMD_OPENMM:?set KAMD_OPENMM to the openmm/ folder of this repository}
PMEMD=${PMEMD:?set PMEMD to the KAMD-patched pmemd.cuda_SPFP}
PY=${PY:-python3}
export PYTHONPATH=$ROOT
mkdir -p $ROOT/validation/fair && cd $ROOT/validation/fair
ln -sf ${AMBER_TEST:?set AMBER_TEST}/prmtop prmtop; ln -sf ${AMBER_TEST}/restrt restrt
N=${N:-100000}
mk() { # $1 = 1 for KAMD
cat <<EOF
matched-condition benchmark
&cntrl
   imin=0, irest=1, ntx=5, ntpr=5000, ntwx=0, ntwr=$N, ntxo=1, nstlim=$N,
   dt=0.002, ntt=3, gamma_ln=5, temp0=355, ig=1060,
   ntc=2, ntf=2, cut=8.0, ntb=2, ntp=1, barostat=2, mcbarint=25, pres0=1.0, nscm=1000,
/
EOF
[ "$1" = 1 ] && cat <<EOF
&kamd
   pseudobd_cycles = 250, shell_cut = 5.0, surface_cut = 3.0, nzero_crossings = 2,
   refer_mask(1) = "@20-96", refer_mask(2) = "@285-339",
   refer_mask(3) = "@366-482", refer_mask(4) = "@568-611",
   solvent_mask  = "@623-29467",
/
EOF
}
mk 0 > plain.in; mk 1 > kamd.in
for c in plain kamd; do
  $PMEMD -O -i $c.in -p prmtop -c restrt -o amber_$c.out -r amber_$c.rst7
  echo "amber_$c $(grep 'ns/day' amber_$c.out | tail -1)"
done
grep "velocity reset, group 1," amber_kamd.out | awk '{print $NF}' > amber_kamd.resets
echo "amber_kamd resets: $(wc -l < amber_kamd.resets)"
A=$(grep "Ewald Coefficient" amber_plain.out | head -1 | awk '{print $NF}')
G=$(grep "NFFT1" amber_plain.out | head -1 | awk '{print $3, $6, $9}')
echo "matching OpenMM PME to Amber: alpha=$A /A grid=$G"
for c in plain kamd; do
  extra=""; [ $c = plain ] && extra="--no-kamd"
  $PY $ROOT/validation/run_openmm_kamd.py --prmtop prmtop --restart restrt --mdin $c.in --nsteps $N $extra \
      --pme-alpha $A --pme-grid $G --cmm-freq 1000 --out omm_$c | grep -E "ns_per_day|n_resets|T_mean|density" | tr -d '\n' | sed "s/^/omm_$c /"; echo
done
