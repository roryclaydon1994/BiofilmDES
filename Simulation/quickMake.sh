#!/bin/bash
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
./Main/findOrientationalStiffness.out testStiffness 30 7 0.0 0 0
# OMP_NUM_THREADS=$NUM_THREADS ./Main/main.out VerletDefault/run0 1e-2 1 2 1
# OMP_NUM_THREADS=$NUM_THREADS ./Main/main.out ChainingDefault/run0 1e-2 1e-2 0.7 100
# OMP_NUM_THREADS=1 ./Main/main.out AdhesionHashTest/run3 5 5e-4 1e-2 2e-4
# OMP_NUM_THREADS=1 ./Main/main.out ChainingPhaseDiagramFinal/run367/repeat4 1e-2 1e-2 0.7 100
