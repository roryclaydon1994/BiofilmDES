#!/bin/bash
#
# Run adhesion parameter sweeps Test 1
#
#SBATCH --job-name=AG43SweepsLowForceThreshold
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END
#SBATCH --partition=long
#SBATCH --time=24:00:00
#SBATCH --mem=500M
#SBATCH --array=[1-30]%16
JJ=$SLURM_ARRAY_TASK_ID
reps=3

# First run num to use
start_from=1

# Get the 0 indexed row and column
row=$(( (${JJ}-1)/${reps} ))
col=$(( (${JJ}-1)%${reps} ))

# Create a run number and repeat
II=$(( $start_from + $row ))
REPEAT=$(( 1 + $col ))
echo "Array job: $JJ Run: $II Repeat: $REPEAT"

# Set up simulation parameters
# Kappa=$(bc -l <<< "(0.6*$II)")
# Kappa=$(bc -l <<< "(0.5+0.1*$row)")
DivLen=5
incr=$(bc -l <<< "( 1/3 )")
GrwthRate=2e-4
Kappa=$(bc -l <<< "( $row*0.05 )")
ForceThresh=$(bc -l <<< "($Kappa * 10^(-3))")

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
cd ./${build_name}
echo "== Preparing MakeFile with CMake =="
cmake ..
echo "== Making =="
make -j $OMP_NUM_THREADS
echo "== Running sim =="
echo "./Main/main.out ${out_dir_name} $DivLen $GrwthRate $Kappa $ForceThresh"
./Main/main.out ${out_dir_name} $DivLen $GrwthRate $Kappa $ForceThresh
echo "== Clearing build file =="
cd ..
rm -rf ./${build_name}
