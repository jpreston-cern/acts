#!/bin/bash
set -e

echo "Sourcing..."
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a_ATLAS_3/x86_64-el9-gcc13-opt/setup.sh

echo "Cmake building..."
cmake -B build -S . -DACTS_BUILD_FATRAS=on -DACTS_BUILD_EXAMPLES_PYTHIA8=on  -DACTS_BUILD_EXAMPLES_PYTHON_BINDINGS=on

echo "building files..."
cmake --build build -j8 &> compile.log
code compile.log

echo "source python binds"
source build/python/setup.sh

echo "run gbtsv2 seeding"
Examples/Scripts/Python/full_chain_itk_Gbts.py &> runtime.log
code runtime.log