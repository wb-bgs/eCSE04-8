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
GLOBLIBI_VERSION=5.0
SLATEC_VERSION=4.1
ERRMSG="Invalid syntax: build-cirrus-intel.sh release|debug"


if [[ "${BUILD}" != "release" && "${BUILD}" != "debug" ]]; then
  echo ${ERRMSG}
  exit
fi

if [[ "${VERSION}" == "" ]]; then
  echo ${ERRMSG}
  exit
fi


module -s load intel-20.4/compilers
module -s load intel-20.4/mpi


PRFX=${HOME/home/work}
WMAM_LABEL=WMAM
WMAM_VERSION=${VERSION}
WMAM_NAME=${WMAM_LABEL}-${WMAM_VERSION}
WMAM_BUILD_ROOT=${PRFX}/projects/eCSE04-8/apps/${WMAM_LABEL}/src
WMAM_INSTALL_ROOT=${PRFX}/apps/${WMAM_LABEL}/${WMAM_VERSION}/INTEL/20.4
WMAM_INSTALL_PATH=${WMAM_INSTALL_ROOT}/${BUILD}


echo -e "\n\nBuilding ${WMAM_LABEL} ${WMAM_VERSION} with globlibi ${GLOBLIBI_VERSION} (${BUILD})...\n\n"
 
SLATEC_ROOT=${PRFX}/libs/slatec/${SLATEC_VERSION}/INTEL/20.4/${BUILD}
GLOBLIBI_ROOT=${PRFX}/libs/globlibi/${GLOBLIBI_VERSION}/INTEL/20.4/${BUILD}

cd ${WMAM_BUILD_ROOT}

cp makefile.cirrus.intel makefile
set_compile_options ./makefile

LIBS_MAKEFILE_LINE="${GLOBLIBI_ROOT}/lib/libgloblibi.a ${SLATEC_ROOT}/lib/libslatec.a"
sed -i "s:LIBS =:LIBS = ${LIBS_MAKEFILE_LINE}:g" ./makefile

rm -f ${WMAM_BUILD_ROOT}/*.o
rm -f ${WMAM_BUILD_ROOT}/*.mod
make

mkdir -p ${WMAM_INSTALL_PATH}/bin

mv ${WMAM_BUILD_ROOT}/mod_wmam_020 ${WMAM_INSTALL_PATH}/bin/
