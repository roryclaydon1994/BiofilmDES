#!/bin/bash
#
# Chaining transition sweep
#
#SBATCH --job-name=ChainingBucklingTransition
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END
#SBATCH --partition=long
#SBATCH --time=24:00:00
#SBATCH --mem=500M
#SBATCH --array=[1-420]%20
JJ=$SLURM_ARRAY_TASK_ID
reps=20

# Set the maximum number of runs below
MaxRuns=420
divs=$(( ${MaxRuns}/${reps} ))

# First run num to use
start_from=1

# Get the 0 indexed row and column
col=$(( (${JJ}-1)/${divs} ))
row=$(( (${JJ}-1)%${divs} ))

# Create a run number and repeat
II=$(( $start_from + $row ))
REPEAT=$(( 1 + $col ))
echo "Array job: $JJ Run: $II Repeat: $REPEAT"

# Set up simulation parameters
Kappa=1
BendRig=1
# LinkingProb=1e-2
LinkingProb=$( bc -l <<< "( 0+0.05*( $II-1 ) )" )
ForceThresh=100

# Set OMP_NUM_THREADS to match --cpus-per-tasks
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

build_name=${SLURM_JOB_NAME}_build_${II}_repeat_${REPEAT}
out_dir_name=${SLURM_JOB_NAME}/run${II}/repeat${REPEAT}

echo "build_name=${build_name}"
echo "out_dir_name=$out_dir_name"

echo "== Clearing old build files =="
rm -rf ./build # clear old build directory (if necessary)
rm -rf ./${build_name}
echo "== Creating new build directory =="
mkdir ${build_name}
cd ./${build_name} || exit
echo "== Preparing MakeFile with CMake =="
# cmake ..
~/cmake-3.23.1/bin/cmake .. -DON_CLUSTER=True
echo "== Making =="
make -j $OMP_NUM_THREADS
echo "== Running sim =="
echo "./Main/main.out ${out_dir_name} $Kappa $BendRig $LinkingProb $ForceThresh"
./Main/main.out ${out_dir_name} $Kappa $BendRig $LinkingProb $ForceThresh
echo "== Clearing build file =="
cd ..
rm -rf ./${build_name}
