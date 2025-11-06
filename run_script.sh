
echo "sourcing"
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a_ATLAS_3/x86_64-el9-gcc14-opt/setup.sh

echo "compiling ACTS"
cmake -B build -S . -DACTS_BUILD_FATRAS=on -DACTS_BUILD_EXAMPLES_PYTHIA8=on  -DACTS_BUILD_EXAMPLES_PYTHON_BINDINGS=on

echo "building acts"
code compile.log
cmake --build build -j8 &> compile.log

echo "sourcing"
source build/python/setup.sh

echo "running GBTS seeding"
code runtime.log
python3 Examples/Scripts/Python/full_chain_itk_Gbts.py &> runtime.log