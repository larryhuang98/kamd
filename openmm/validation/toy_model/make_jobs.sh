#!/bin/bash
# Writes jobs.txt: one run_toy.py command per line (run 4 at a time with xargs -P 4).
cd "$(dirname "$0")"
mkdir -p runs
J=jobs.txt; : > $J
# main matrix: 8 seeds x 2.5 ns, gamma = 1/ps (interleaved by seed)
for s in 1 2 3 4 5 6 7 8; do
  echo "python3 run_toy.py --mode plain  --seed $s --ns 2.5 --out runs/plain_s$s.npz" >> $J
  echo "python3 run_toy.py --mode kamd   --seed $s --ns 2.5 --out runs/kamd_def_s$s.npz" >> $J
  echo "python3 run_toy.py --mode kamd   --pbd 20 --nz 1 --seed $s --ns 2.5 --out runs/kamd_freq_s$s.npz" >> $J
  echo "python3 run_toy.py --mode mbonly --pbd 20 --nz 1 --seed $s --ns 2.5 --out runs/mbonly_s$s.npz" >> $J
done
# friction scan (conditions 1 and 3) and separation-vector monitor: 4 seeds x 1.2 ns
for s in 1 2 3 4; do
  for g in 0.1 10; do
    echo "python3 run_toy.py --mode plain --friction $g --seed $s --ns 1.2 --out runs/plain_g${g}_s$s.npz" >> $J
    echo "python3 run_toy.py --mode kamd --pbd 20 --nz 1 --friction $g --seed $s --ns 1.2 --out runs/kamd_freq_g${g}_s$s.npz" >> $J
  done
  echo "python3 run_toy.py --mode sepmon --pbd 20 --nz 1 --seed $s --ns 1.2 --out runs/sepmon_freq_s$s.npz" >> $J
done
