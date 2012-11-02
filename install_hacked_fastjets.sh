#!/bin/bash

# copy local fastjets to build directory
echo "**********************************************"
echo "This is a dumb script, read it before running!"
echo -e "**********************************************\n"

cp ./include/Rivet/Projections/FastJets.hh ../../build/rivet/include/Rivet/Projections/FastJets.hh
cp ./src/Projections/FastJets.cc ../../build/rivet/src/Projections/FastJets.cc
cd ../../build/rivet
echo "I'm rebuilding FastJets for you"
make -j 2 && make install 
echo "Done, returning you to your old pwd"
cd -

