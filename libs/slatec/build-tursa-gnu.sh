#!/bin/bash


function set_compile_options {
  MAKEFILE=$1
  if [[ "${BUILD}" == "release" ]]; then
    sed -i "s:FFLAGS= -O2:FFLAGS= -O3 -std=legacy -fdefault-real-8 -fdefault-double-8:g" ${MAKEFILE}
  elif [[ "${BUILD}" == "debug" ]]; then
    sed -i "s:FFLAGS= -O2:FFLAGS= -g -O0 -std=legacy -fdefault-real-8 -fdefault-double-8 -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace:g" ${MAKEFILE}
  fi
}


BUILD=$1
GCC_VERSION=12.2.0
ERRMSG="Invalid syntax: build-tursa-gnu.sh release|debug"


if [[ "${BUILD}" != "release" && "${BUILD}" != "debug" ]]; then
  echo ${ERRMSG}
  exit
fi


module load gcc/${GCC_VERSION}

PRFX=${HOME}
SLATEC_LABEL=slatec
SLATEC_VERSION=4.1
SLATEC_BUILD_ROOT=${PRFX}/projects/eCSE04-8/libs/${SLATEC_LABEL}/src
SLATEC_INSTALL_ROOT=${PRFX}/libs/${SLATEC_LABEL}/${SLATEC_VERSION}/GNU/${GCC_VERSION}
SLATEC_INSTALL_PATH=${SLATEC_INSTALL_ROOT}/${BUILD}

echo -e "\n\nBuilding ${SLATEC_LABEL} ${SLATEC_VERSION} (${BUILD}) using GNU ${GCC_VERSION}...\n\n"
  

cd ${SLATEC_BUILD_ROOT}

cp makefile.tursa.gnu makefile
sed -i "s:/usr/local:${SLATEC_INSTALL_PATH}:g" ${SLATEC_BUILD_ROOT}/makefile
sed -i "s:ldconfig:/sbin/ldconfig -C ${SLATEC_INSTALL_PATH}/etc/ld.so.cache:g" ${SLATEC_BUILD_ROOT}/makefile

cp ./static/makefile.sav ./static/makefile
cp ./dynamic/makefile.sav ./dynamic/makefile

set_compile_options ${SLATEC_BUILD_ROOT}/static/makefile
set_compile_options ${SLATEC_BUILD_ROOT}/dynamic/makefile


rm -rf ${SLATEC_INSTALL_PATH}
mkdir -p ${SLATEC_INSTALL_PATH}/lib
mkdir -p ${SLATEC_INSTALL_PATH}/etc
mkdir -p ${SLATEC_INSTALL_PATH}/man/man1

FC=gfortran make
make install
make clean
