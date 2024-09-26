#!/bin/bash


function set_compile_options {
  MAKEFILE=$1
  if [[ "${BUILD}" == "release" ]]; then
    sed -i "s:FFLAGS= -O2:FFLAGS= -O3:g" ${MAKEFILE}
  elif [[ "${BUILD}" == "debug" ]]; then
    sed -i "s:FFLAGS= -O2:FFLAGS= -g -O0:g" ${MAKEFILE}
  fi
}


BUILD=$1
ERRMSG="Invalid syntax: build-cirrus-intel.sh release|debug"


if [[ "${BUILD}" != "release" && "${BUILD}" != "debug" ]]; then
  echo ${ERRMSG}
  exit
fi


module -s load intel-20.4/compilers


PRFX=${HOME/home/work}
SLATEC_LABEL=slatec
SLATEC_VERSION=4.1
SLATEC_BUILD_ROOT=${PRFX}/projects/eCSE04-8/libs/${SLATEC_LABEL}/src
SLATEC_INSTALL_ROOT=${PRFX}/libs/${SLATEC_LABEL}/${SLATEC_VERSION}/INTEL/20.4
SLATEC_INSTALL_PATH=${SLATEC_INSTALL_ROOT}/${BUILD}

echo -e "\n\nBuilding ${SLATEC_LABEL} ${SLATEC_VERSION} (${BUILD}) using Intel 20.4...\n\n"


cd ${SLATEC_BUILD_ROOT}

cp makefile.cirrus.intel makefile
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

FC=ifort make
make install
make clean
