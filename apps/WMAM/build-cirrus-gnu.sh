#!/bin/bash


function set_compile_options {
  MAKEFILE=$1
  if [[ "${BUILD}" == "release" ]]; then
    sed -i "s:FFLAGS =:FFLAGS = -O3 -cpp -fopenmp -Wno-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8:g" ${MAKEFILE}
  elif [[ "${BUILD}" == "debug" ]]; then
    sed -i "s:FFLAGS =:FFLAGS = -g -O0 -cpp -fopenmp -Wno-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8 -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace:g" ${MAKEFILE}
  fi
}


BUILD=$1
VERSION=5.0
SLATEC_VERSION=4.1
GCC_VERSION=10.2.0
ERRMSG="Invalid syntax: build-cirrus-gnu.sh release|debug"


if [[ "${BUILD}" != "release" && "${BUILD}" != "debug" ]]; then
  echo ${ERRMSG}
  exit
fi


module -s load openmpi/4.1.6


PRFX=${HOME/home/work}
WMAM_LABEL=WMAM
WMAM_VERSION=${VERSION}
WMAM_NAME=${WMAM_LABEL}-${WMAM_VERSION}
WMAM_BUILD_ROOT=${PRFX}/projects/eCSE04-8/apps/${WMAM_LABEL}/src
WMAM_INSTALL_ROOT=${PRFX}/apps/${WMAM_LABEL}/${WMAM_VERSION}/GNU/${GCC_VERSION}
WMAM_INSTALL_PATH=${WMAM_INSTALL_ROOT}/${BUILD}

echo -e "\n\nBuilding ${WMAM_LABEL} ${WMAM_VERSION} (${BUILD}) using GNU ${GCC_VERSION}...\n\n"
  

SLATEC_ROOT=${PRFX}/libs/slatec/${SLATEC_VERSION}/GNU/${GCC_VERSION}/${BUILD}

cd ${WMAM_BUILD_ROOT}

cp makefile.cirrus.gnu makefile
set_compile_options ./makefile

LIBS_MAKEFILE_LINE="${SLATEC_ROOT}/lib/libslatec.a"
sed -i "s:LIBS =:LIBS = ${LIBS_MAKEFILE_LINE}:g" ./makefile


find ${WMAM_BUILD_ROOT} -name '*.o' -delete
find ${WMAM_BUILD_ROOT} -name '*.mod' -delete
make
find ${WMAM_BUILD_ROOT} -name '*.o' -delete
find ${WMAM_BUILD_ROOT} -name '*.mod' -delete

mkdir -p ${WMAM_INSTALL_PATH}/bin

mv ${WMAM_BUILD_ROOT}/mod_wmam_020 ${WMAM_INSTALL_PATH}/bin/
