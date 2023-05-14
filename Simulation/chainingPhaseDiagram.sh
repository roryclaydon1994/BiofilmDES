#!/bin/bash
#
# Chaining transition sweep
#
#SBATCH --job-name=AspectGrowthSweep
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END
#SBATCH --partition=long
#SBATCH --time=00:10:00
#SBATCH --mem=500M
#SBATCH --output=R1_AspectGrowthSweep_repeat0_%j.out
#SBATCH --array=[1-100]%20
JJ=$SLURM_ARRAY_TASK_ID

# Set the maximum number of runs below
Ny=10 # number of bendings
Nx=10 # number of pls
Nz=1

idx=$(( ${JJ}-1 )) # index 0 indexed
z=$(( ${idx}/(${Ny}*${Nx})  ))  # repeat number
tmp=$(( ${idx}-${z}*(${Ny}*${Nx}) ))
y=$(( ${tmp}/${Nx} )) # bending rigidity index
x=$(( ${tmp}%${Nx} )) # pl index

# First run num to use
start_from=1

# Create a run number and repeat
II=$(( $x + ${Nx}*${y} ))
# REPEAT=$(( $z + 2 ))
REPEAT=$(( 0 + $z ))
echo "Array job: $JJ Run: $II Repeat: $REPEAT"

### Set up simulation parameters in simulation units
# chaining
Kappa=1
BendRig=1
LinkingProb=0
ForceThresh=100

AspectRatio=$( bc -l <<< "( 1.5+0.5*$y )" )
GrowthRate=$( bc -l <<< "( 0.0002*( $x + 1 ) )" )

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
# cmake ..
~/cmake-3.23.1/bin/cmake .. -DON_CLUSTER=True
echo "== Making =="
make -j $OMP_NUM_THREADS
echo "== Running sim =="
echo "./Main/main.out ${out_dir_name} $Kappa $BendRig $LinkingProb $ForceThresh $AspectRatio $GrowthRate"
./Main/main.out ${out_dir_name} $Kappa $BendRig $LinkingProb $ForceThresh $AspectRatio $GrowthRate
echo "== Clearing build file =="
cd ..
rm -rf ./${build_name}
