#!/bin/bash


function set_compile_options {
  MAKEFILE=$1
  if [[ "${BUILD}" == "release" ]]; then
    sed -i "s:FFLAGS =:FFLAGS = -I${MPI_HOME}/include -O3 -cpp -fopenmp -DOMP_OFFLOAD -DOMP_OFFLOAD_CPTP -DOMP_OFFLOAD_SSQGH -gpu=cuda12.4 -Minfo=mp -mp=gpu -target=gpu -tp=cascadelake -r8:g" ${MAKEFILE}
  elif [[ "${BUILD}" == "debug" ]]; then
    sed -i "s:FFLAGS =:FFLAGS = -I${MPI_HOME}/include -g -O0 -cpp -fopenmp -DOMP_OFFLOAD -DOMP_OFFLOAD_CPTP -DOMP_OFFLOAD_SSQGH -gpu=cuda12.4 -Minfo=mp -mp=gpu -target=gpu -tp=cascadelake -r8 -C -Mnobounds -ffpe-trap=invalid,zero,overflow,underflow,inexact -Ktrap=divz,denorm,inexact,inv,ovf,unf:g" ${MAKEFILE}
  fi
}


BUILD=$1
VERSION=5.0
GCC_VERSION=10.2.0
NVHPC_VERSION=24.5
ERRMSG="Invalid syntax: build-cirrus-nvfortran.sh release|debug"


if [[ "${BUILD}" != "release" && "${BUILD}" != "debug" ]]; then
  echo ${ERRMSG}
  exit
fi


module -s load nvidia/nvhpc-nompi/${NVHPC_VERSION}
module -s load openmpi/4.1.6-cuda-12.4-nvfortran


PRFX=${HOME/home/work}
GLOBLIB_LABEL=globlibi
GLOBLIB_VERSION=${VERSION}
GLOBLIB_NAME=${GLOBLIB_LABEL}-${GLOBLIB_VERSION}
GLOBLIB_BUILD_ROOT=${PRFX}/projects/eCSE04-8/libs/${GLOBLIB_LABEL}/src
GLOBLIB_INSTALL_ROOT=${PRFX}/libs/${GLOBLIB_LABEL}/${GLOBLIB_VERSION}/GNU/${GCC_VERSION}-nvfortran
GLOBLIB_INSTALL_PATH=${GLOBLIB_INSTALL_ROOT}/${BUILD}/lib

echo -e "\n\nBuilding ${GLOBLIB_LABEL} ${GLOBLIB_VERSION} (${BUILD}) using nvfortran (NVHPC ${NVHPC_VERSION}) with GNU ${GCC_VERSION}...\n\n"


cd ${GLOBLIB_BUILD_ROOT}

cp makefile.cirrus.nvfortran makefile
sed -i "s:libdir =:libdir = ${GLOBLIB_INSTALL_PATH}:g" ./makefile

set_compile_options ./makefile

make
make install
make clean
