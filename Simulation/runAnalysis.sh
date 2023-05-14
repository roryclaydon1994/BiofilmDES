#!/usr/bin/env bash
DirName=$1
# echo "=== Removing old figure files ==="
# rm ../GeneratedOutput/animation/${DirName}/*
# rm ../GeneratedOutput/figures/${DirName}/*
#

# Compare for a single sweep through the colony
# python3 ../Analysis/analysis_master.py \
#         -S -1 \
#         --max-cores 12 \
#         --compare-only \
#         --data-dir "SimOutput/${DirName}" \
#         --analysis-dir "AnalysisResultsTarly/${DirName}"
# python3 ../Analysis/analysis_master.py \
#           -S -1 \
#           --data-dir "SimOutput/${DirName}" \
#           --analysis-dir "AnalysisResultsTarly/${DirName}"

# find the individual descriptors and then compare all on a phase diagram
python3 ../Analysis/analysis_master.py \
          -E -1 \
          --phase\
          --max-cores 12 \
          --data-dir "SimOutput/${DirName}" \
          --analysis-dir "AnalysisResultsTarly/${DirName}"

# This is only for the energy plots
# python3 ../Analysis/analysis_master.py \
#           -S 0 \
#           -E -1 \
#           --data-dir "SimOutput/${DirName}" \
#           --analysis-dir "AnalysisResults/${DirName}"
