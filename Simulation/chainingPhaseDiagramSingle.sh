#!/bin/bash
#
# Chaining transition sweep
#
#SBATCH --job-name=ChainingPhaseDiagramFinal
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END
#SBATCH --partition=long
#SBATCH --time=36:00:00
#SBATCH --mem=500M
#SBATCH --output=repeat0_%j.out
JJ=379

# Set the maximum number of runs below
Ny=19 # number of bendings
Nx=20 # number of pls
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

# Set up simulation parameters
Kappa=1
BendRig=$( bc -l <<< "( 0.02*( 1.5^$y ) )" )
LinkingProb=$( bc -l <<< "( 0.05*( $x + 1 ) )" )
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
cd ./${build_name}
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
