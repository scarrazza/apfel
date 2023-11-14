#!/bin/sh -l
set -x
brew install wget coreutils gcc swig
brew tap davidchall/hep
brew install lhapdf
cmake -S . -B BUILD -DCMAKE_INSTALL_PREFIX=$(pwd)/../INSTALL -DCMAKE_Fortran_COMPILER=gfortran-13  -DCMAKE_CXX_COMPILER=g++-13 -DCMAKE_C_COMPILER=gcc-13 || exit 1
cmake --build BUILD -j 2 || exit 1
cmake --install BUILD || exit 1
ctest --test-dir BUILD -j 2 || exit 1
