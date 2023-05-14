#!/usr/bin/env bash

NUM_THREADS=$(nproc --all)
echo "== Clearing old build files =="
rm -rf ./build # clear old build directory (if necessary)
echo "== Creating new build directory =="
mkdir build
cd ./build
echo "== Preparing MakeFile with CMake =="
cmake ..
echo "== Making =="
make -j $NUM_THREADS
echo "== Running test =="
echo "Change back to actually run tests"

NUM_THREADS=4
dir=ScalingOrientationWidthRelaxPackingFraction
radius=20
linking_prob=0.0
bend_rig=1
kappa=1

for jj in {1..1}
do
  for ii in {1..55..5}
  do
    tmp_radius=$( bc -l <<< "( 20*1.1^( ${ii}-1 ) )" )
    tmp_lnk_prb=$( bc -l <<< "( 0 )" )
    tmp_width=$( bc -l <<< "( 7+2*( ${jj}-1 )/2 )" )
    echo "./Main/findOrientationalStiffness.out ${dir} ${tmp_radius} ${tmp_width} ${tmp_lnk_prb} ${bend_rig} ${kappa}"
    OMP_NUM_THREADS=$NUM_THREADS ./Main/findOrientationalStiffness.out \
                                  ${dir} ${tmp_radius} ${tmp_width} ${tmp_lnk_prb} \
                                  ${bend_rig} ${kappa}
  done
done
