#!/usr/bin/env bash
DirName=$1
# echo "=== Removing old figure files ==="
# rm ../GeneratedOutput/animation/${DirName}/*
# rm ../GeneratedOutput/figures/${DirName}/*

python3 ../Analysis/visualise_master.py \
          -S -1 \
          --sus-vis-radius-factor 0.7 \
          --data-dir "SimOutput/${DirName}" \
          --fig-format pdf \
          --anim-dir "figures/${DirName}"

# Uncomment to make a movie
python3 ../Analysis/visualise_master.py \
          -S 0  \
          -E -1 \
          --sus-vis-radius-factor 0.7 \
          --data-dir "SimOutput/${DirName}" \
          --anim-dir "animation/${DirName}" \
         --max-cores 6
