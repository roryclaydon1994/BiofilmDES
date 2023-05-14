#!/usr/bin/env bash

repo_dir="/storage/datastore-personal/s1836309/BiofilmDES/SimOutput/SCOWF_Perturbation_K0"
local_dir="/media/rory/Elements/ScalingOrientation/"
# rsync -abzuvi --exclude="*/radius_*_width_*_linking_prob_*_kappa_1_B_1/biofilm_*.dat" \
#                $TARLY:$dir ../GeneratedOutput/SimOutput/
rsync -abzuvim --exclude="*/*width_-*" $TARLY:$repo_dir $local_dir

# To copy the whole of the AnalysisResults in datastore, leave run_name blank
# run_name=ChainingBucklingTransition
# dir="/storage/datastore-personal/s1836309/BiofilmDES/AnalysisResults/${run_name}"
# rsync -abzuvi $TARLY:$dir ../GeneratedOutput/AnalysisResultsTarly/

# To copy the whole of the SimOutput in datastore
# dir="/storage/datastore-personal/s1836309/BiofilmDES/SimOutput/"
# rsync -abzuvi $TARLY:$dir ../GeneratedOutput/SimOutput/
