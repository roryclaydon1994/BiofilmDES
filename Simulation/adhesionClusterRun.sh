#!/bin/bash
#
# Run adhesion parameter sweeps Test 1
#
#SBATCH --job-name=MyJobName
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --mem=500M
#SBATCH --array=[1-5]%2
II=$SLURM_ARRAY_TASK_ID
echo "Array job: $II"

kappa_depletion=1e-2
# kappa_depletion=$(bc -l <<< "(0.001*4^($II-1))")

Rd=1

Ri=2
# Ri=$(bc -l <<< "2.5 + 0.5*$II")

rep_strength=1
# rep_strength=$(bc -l <<< "(2^$II)")

# Set OMP_NUM_THREADS to match --cpus-per-tasks
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

build_name=${SLURM_JOB_NAME}_build_$II
out_dir_name=${SLURM_JOB_NAME}/run$II

echo "build_name=$build_name"
echo "out_dir_name=$out_dir_name"

echo "== Clearing old build files =="
rm -rf ./build # clear old build directory (if necessary)
rm -rf ./$build_name
echo "== Creating new build directory =="
mkdir $build_name
cd ./$build_name
echo "== Preparing MakeFile with CMake =="
cmake ..
echo "== Making =="
make -j $OMP_NUM_THREADS
echo "== Running sim =="
./Main/main.out $out_dir_name $kappa_depletion $Rd $Ri $rep_strength
echo "== Clearing build file =="
cd ..
rm -rf ./$build_name
