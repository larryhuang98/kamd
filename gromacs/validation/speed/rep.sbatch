#!/bin/bash
#SBATCH -J kamd_gmx_rep
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 00:50:00
#SBATCH -o rep_%j.log
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
sed -e "s/^log .*/log = kamd_rep.log/" kamd.dat > kamd_rep.dat
for upd in gpu cpu; do for r in 1 2 3; do for k in plain kamd; do
  extra=""; [ $k = kamd ] && extra="-kamd kamd_rep.dat"
  gmx mdrun -s md.tpr $extra -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update $upd -resethway -deffnm r_${upd}_$k > /dev/null 2>&1
  echo "$upd $k $r $(grep Performance r_${upd}_$k.log | awk '{print $2}')"
done; done; done
