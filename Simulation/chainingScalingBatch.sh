#!/bin/bash
#
# Chaining transition sweep
#
#SBATCH --job-name=ScalingExpGrwthRate
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END
#SBATCH --partition=long
#SBATCH --time=24:00:00
#SBATCH --mem=500M
#SBATCH --array=[1-26]
#SBATCH --output=ScalingExpGrwthRate_%j.out
JJ=$SLURM_ARRAY_TASK_ID

# Set up simulation parameters
# Kappa=$( bc -l <<< "( ($JJ-1) / 10 )" )
# BendRig=$( bc -l <<< "( ($JJ-1) / 10 )" )

Kappa=1
BendRig=1
GrwthRate=$( bc -l <<< "( 10 * $JJ / 26 )" )

# Set OMP_NUM_THREADS to match --cpus-per-tasks
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

build_name=${SLURM_JOB_NAME}_build_${JJ}
out_dir_name=${SLURM_JOB_NAME}/

echo "build_name=${build_name}"
echo "out_dir_name=$out_dir_name"

echo "== Clearing old build files =="
rm -rf ./build # clear old build directory (if necessary)
rm -rf ./${build_name}
echo "== Creating new build directory =="
mkdir ${build_name}
cd ./${build_name}
echo "== Preparing MakeFile with CMake =="
# cmake ..
~/cmake-3.23.1/bin/cmake .. -DON_CLUSTER=True
echo "== Making =="
make -j $OMP_NUM_THREADS
echo "== Running sim =="
echo "./Main/findScaling.out ${out_dir_name} $Kappa $BendRig $GrwthRate"
./Main/findScaling.out ${out_dir_name} $Kappa $BendRig $GrwthRate
echo "== Clearing build file =="
cd ..
rm -rf ./${build_name}
