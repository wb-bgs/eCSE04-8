#!/bin/bash


function set_compile_options {
  MAKEFILE=$1
  if [[ "${BUILD}" == "release" ]]; then
    sed -i "s:FFLAGS =:FFLAGS = -O3 -fopenmp:g" ${MAKEFILE}
  elif [[ "${BUILD}" == "debug" ]]; then
    sed -i "s:FFLAGS =:FFLAGS = -g -O0 -fopenmp:g" ${MAKEFILE}
  fi
}


BUILD=$1
VERSION=5.0
ERRMSG="Invalid syntax: build-intel.sh release|debug"


if [[ "${BUILD}" != "release" && "${BUILD}" != "debug" ]]; then
  echo ${ERRMSG}
  exit
fi


module -s load intel-20.4/compilers
module -s load intel-20.4/mpi


PRFX=${HOME/home/work}
GLOBLIB_LABEL=globlibi
GLOBLIB_VERSION=${VERSION}
GLOBLIB_NAME=${GLOBLIB_LABEL}-${GLOBLIB_VERSION}
GLOBLIB_BUILD_ROOT=${PRFX}/projects/eCSE04-8/libs/${GLOBLIB_LABEL}/src
GLOBLIB_INSTALL_ROOT=${PRFX}/libs/${GLOBLIB_LABEL}/${GLOBLIB_VERSION}/INTEL/20.4
GLOBLIB_INSTALL_PATH=${GLOBLIB_INSTALL_ROOT}/${BUILD}/lib


echo -e "\n\nBuilding ${GLOBLIB_LABEL} ${GLOBLIB_VERSION} (${BUILD})...\n\n"
 

cd ${GLOBLIB_BUILD_ROOT}

cp makefile.cirrus.intel makefile
sed -i "s:libdir =:libdir = ${GLOBLIB_INSTALL_PATH}:g" ./makefile

set_compile_options ./makefile

make
make install
make clean
