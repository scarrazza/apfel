#!/bin/sh -l
set -x
brew install wget coreutils gcc swig
brew tap davidchall/hep
brew install lhapdf
cmake -S . -B BUILD -DCMAKE_INSTALL_PREFIX=$(pwd)/../INSTALL
cmake --build BUILD -j 2
cmake --install BUILD
ctest --test-dir BUILD -j 2
