#!/bin/bash

PRFX=${HOME/home/work}
GMT_LABEL=gmt
GMT_VERSION=6.3.0


module -q restore
module -q load PrgEnv-gnu

module -q load cmake

module -q load cray-hdf5
module -q load cray-netcdf
module -q load cray-fftw


cd ${PRFX}/eCSE04-8/utils/gmt/src

mkdir build
cd build

export CC=cc
export CXX=CC
export FC=ftn

cmake ../ -DCMAKE_INSTALL_PREFIX=${PRFX}/utils/${GMT_LABEL}/${GMT_VERSION}
make -j 8 install

cd ..
rm -rf build
