#!/bin/bash
#
# Chaining transition sweep
#
#SBATCH --job-name=SCOWF_Perturbation
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END
#SBATCH --partition=long
#SBATCH --time=168:00:00
#SBATCH --mem=500M
#SBATCH --array=[1-1000]
#SBATCH --output=SCOWF1_%j.out
JJ=$SLURM_ARRAY_TASK_ID

NWdths=10
NRadii=100
R_indx=$(( (${JJ}-1)/${NWdths} ))
W_indx=$(( (${JJ}-1)%${NWdths} ))

P_indx=$1

# Set up simulation parameters
radius=10
bend_rig=1
kappa=1

tmp_radius=$( bc -l <<< "( ${radius}*1.03^(${R_indx})  )" )
tmp_lnk_prb=$( bc -l <<< "( ${P_indx}*0.05 )" )
tmp_width=$( bc -l <<< "( 5+( 4 / ( ${NWdths}-1 ) ) * ( ${W_indx} ) )" )

# Set OMP_NUM_THREADS to match --cpus-per-tasks
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

build_name=${SLURM_JOB_NAME}_build_${JJ}_Pindx_${P_indx}_Windx_${W_indx}_Rindx_${R_indx}_R_${tmp_radius}
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
echo "./Main/findOrientationalStiffness.out ${out_dir_name} ${tmp_radius} ${tmp_width} ${tmp_lnk_prb} ${bend_rig} ${kappa}"
./Main/findOrientationalStiffness.out ${out_dir_name} ${tmp_radius} ${tmp_width} ${tmp_lnk_prb} ${bend_rig} ${kappa}
echo "== Clearing build file =="
cd ..
rm -rf ./${build_name}
