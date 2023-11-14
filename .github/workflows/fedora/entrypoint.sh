#!/bin/sh -l
set -x
uname -a 
cat /etc/issue
yum -y install  gcc gcc-c++ gcc-gfortran make wget which cmake cmake-data cmake-filesystem lhapdf lhapdf-devel
out=0
cmake -S . -B BUILD -DCMAKE_INSTALL_PREFIX=$(pwd)/../INSTALL || out=1
cmake --build BUILD -j 2 || out=1
cmake --install BUILD || out=1
ctest --test-dir BUILD -j 2 || out=1
echo ::set-output fineapfel=out::$out
