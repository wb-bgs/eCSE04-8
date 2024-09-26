#!/bin/bash


function set_compile_options {
  MAKEFILE=$1
  if [[ "${BUILD}" == "release" ]]; then
    sed -i "s:FFLAGS =:FFLAGS = -O3 -cpp -fopenmp -foffload=nvptx-none -Wno-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8:g" ${MAKEFILE}
  elif [[ "${BUILD}" == "debug" ]]; then
    sed -i "s:FFLAGS =:FFLAGS = -g -O0 -cpp -fopenmp -foffload=nvptx-none -Wno-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8 -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace:g" ${MAKEFILE}
  fi
}


BUILD=$1
VERSION=5.0
GCC_VERSION=12.3.0
NVHPC_VERSION=24.5
ERRMSG="Invalid syntax: build-cirrus-gnu-offload.sh release|debug"


if [[ "${BUILD}" != "release" && "${BUILD}" != "debug" ]]; then
  echo ${ERRMSG}
  exit
fi


module -s load nvidia/nvhpc-nompi/${NVHPC_VERSION}
module -s load openmpi/4.1.6-cuda-12.4 
module -s swap -f gcc gcc/${GCC_VERSION}-offload


PRFX=${HOME/home/work}
GLOBLIB_LABEL=globlibi
GLOBLIB_VERSION=${VERSION}
GLOBLIB_NAME=${GLOBLIB_LABEL}-${GLOBLIB_VERSION}
GLOBLIB_BUILD_ROOT=${PRFX}/projects/eCSE04-8/libs/${GLOBLIB_LABEL}/src
GLOBLIB_INSTALL_ROOT=${PRFX}/libs/${GLOBLIB_LABEL}/${GLOBLIB_VERSION}/GNU/${GCC_VERSION}-offload
GLOBLIB_INSTALL_PATH=${GLOBLIB_INSTALL_ROOT}/${BUILD}/lib

echo -e "\n\nBuilding ${GLOBLIB_LABEL} ${GLOBLIB_VERSION} (${BUILD}) using GNU ${GCC_VERSION}-offload...\n\n"


cd ${GLOBLIB_BUILD_ROOT}

cp makefile.cirrus.gnu.offload makefile
sed -i "s:libdir =:libdir = ${GLOBLIB_INSTALL_PATH}:g" ./makefile

set_compile_options ./makefile

make
make install
make clean
