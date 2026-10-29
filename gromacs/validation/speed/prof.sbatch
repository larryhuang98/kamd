#!/bin/bash
#SBATCH -J kamd_gmx_prof
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 00:40:00
#SBATCH -o prof_%j.log
# Site-specific settings come from environment variables:
#   PMEMD        KAMD-patched pmemd.cuda_SPFP          AMBER_TEST  amber/test of this repository
#   KAMD_OPENMM  openmm/ of this repository            PY          python with openmm, parmed, scipy
#   KAMD_GMX_DIR directory holding install/ (KAMD-patched GROMACS) and validation/
# Load your compiler/CUDA/Amber runtime libraries (LD_LIBRARY_PATH) before submitting.
for i in $(seq 30); do [ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ] && break; sleep 10; done
G=${KAMD_GMX_DIR:?set KAMD_GMX_DIR (see header)}
source $G/install/bin/GMXRC
cd $G/validation/speed
echo "host $(hostname)"
sed -e "s/^nzero-crossings.*/nzero-crossings = 100000000/" -e "s/^log .*/log = kamd_mononly.log/" kamd.dat > kamd_mononly.dat
sed -e "s/^log .*/log = kamd_full.log/" kamd.dat > kamd_full.dat
run() { gmx mdrun -s md.tpr $2 -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update gpu -resethway -deffnm p_$1 > p_$1.out 2>&1; echo "$1 $(grep Performance p_$1.log | awk '{print $2}') ns/day"; }
run plainA ""
run mononly "-kamd kamd_mononly.dat"
run full "-kamd kamd_full.dat"
run plainB ""
run full2 "-kamd kamd_full.dat"
