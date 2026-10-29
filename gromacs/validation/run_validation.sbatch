#!/bin/bash
#SBATCH -J kamd_gmx_val
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
# GROMACS KAMD validation on the 29,467-atom Amber KAMD test system (converted with ParmEd).
for i in $(seq 30); do [ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ] && break; sleep 10; done
echo "host $(hostname) $(nvidia-smi --query-gpu=name --format=csv,noheader) cores $(nproc) $(date)"
G=${KAMD_GMX_DIR:?set KAMD_GMX_DIR (see header)}
V=$G/validation
source $G/install/bin/GMXRC
PY=${PY:-python3}
NT=${NT:-8}
cd $V

echo "=== 1. Langevin (sd) run with KAMD, 250000 steps, first-reset diagnostic on"
mkdir -p sd_kamd && cd sd_kamd
gmx grompp -f ../sd.mdp -c ../system.gro -p ../system.top -o sd.tpr -maxwarn 2 > grompp.log 2>&1 || { tail -30 grompp.log; exit 1; }
sed -e "s/^log .*/log = kamd.log/" ../kamd.dat > kamd.dat; echo "diag = yes" >> kamd.dat; ln -sf ../index.ndx .
gmx mdrun -s sd.tpr -kamd kamd.dat -ntmpi 1 -ntomp $NT -nb gpu -pme gpu -update cpu -deffnm sd > mdrun.out 2>&1 || { tail -40 mdrun.out; exit 1; }
grep -c reset kamd.log; grep Performance sd.log
echo "=== 2. first-shell parity vs kamd_openmm reference"
PYTHONPATH=${KAMD_OPENMM:?} $PY ../parity_vs_reference.py kamd_diag.txt index.ndx
cd $V

echo "=== 3. speed: md + v-rescale, 100000 steps"
mkdir -p speed && cd speed
gmx grompp -f ../md.mdp -c ../system.gro -p ../system.top -o md.tpr -maxwarn 2 > grompp.log 2>&1 || { tail -30 grompp.log; exit 1; }
ln -sf ../index.ndx .; cp ../kamd.dat kamd.dat
for upd in gpu cpu; do
  for k in plain kamd; do
    extra=""; [ $k = kamd ] && extra="-kamd kamd.dat"
    sed -i "s/^log .*/log = kamd_$upd.log/" kamd.dat
    gmx mdrun -s md.tpr $extra -ntmpi 1 -ntomp $NT -nb gpu -pme gpu -bonded gpu -update $upd -resethway -deffnm ${upd}_$k > ${upd}_$k.out 2>&1 || { tail -30 ${upd}_$k.out; }
    echo "update=$upd $k $(grep Performance ${upd}_$k.log) resets=$( [ $k = kamd ] && grep -c reset kamd_$upd.log )"
  done
done
echo "done $(date)"
