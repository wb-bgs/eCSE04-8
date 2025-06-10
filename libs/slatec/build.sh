#!/bin/bash


function set_compile_options {
  MAKEFILE=$1
  if [[ "${BUILD}" == "release" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -O3:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -O3 -fallow-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -O1:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "debug" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g -O0:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g -O0 -fallow-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8 -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g -O0:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "craypat" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g -O3 -h profile_generate:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g -O3 -fallow-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8 -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g -O1:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "armmap" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -G2 -O3 -h ipa0:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g1 -O3 -fallow-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8 -fno-inline -fno-optimize-sibling-calls:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g1 -O1:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "scorep" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -G2 -O3 -h ipa0:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g1 -O3 -fallow-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8 -fno-inline -fno-optimize-sibling-calls:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g1 -O1:g" ${MAKEFILE}
    fi
  fi
}


PE_RELEASE=$1
PRGENV=$2
BUILD=$3
ERRMSG="Invalid syntax: build.sh <CPE release> cray|gnu|aocc release|debug|craypat|armmap|scorep"


if [[ "${PRGENV}" != "cray" && "${PRGENV}" != "gnu" && "${PRGENV}" != "aocc" ]]; then
  echo ${ERRMSG}
  exit
fi

if [[ "${BUILD}" != "release" && "${BUILD}" != "debug" && "${BUILD}" != "craypat" && "${BUILD}" != "armmap" && "${BUILD}" != "scorep" ]]; then
  echo ${ERRMSG}
  exit
fi


PRFX=${HOME/home/work}
SLATEC_LABEL=slatec
SLATEC_VERSION=4.1
SLATEC_BUILD_ROOT=${PRFX}/projects/eCSE04-8/libs/${SLATEC_LABEL}/src
SLATEC_INSTALL_ROOT=${PRFX}/libs/${SLATEC_LABEL}/${SLATEC_VERSION}


echo -e "\n\nBuilding ${SLATEC_LABEL} ${SLATEC_VERSION} (${BUILD}) using ${PRGENV} (CPE ${PE_RELEASE}) programming environment...\n\n"
  
module -q restore
module -q load cpe/${PE_RELEASE}
module -q load PrgEnv-${PRGENV}

if [[ "${BUILD}" == "craypat" ]]; then
  module -q load perftools-base
  module -q load perftools
elif [[ "${BUILD}" == "scorep" ]]; then
  module -q load other-software
  module -q unload perftools-base
  if [[ "${PRGENV}" == "cray" ]]; then
    module -q load scalasca/2.6.1-cray
  elif [[ "${PRGENV}" == "gnu" ]]; then
    module -q load scalasca/2.6.1-gcc11
  elif [[ "${PRGENV}" == "aocc" ]]; then
    module -q load scalasca/2.6.1-aocc
  fi
fi

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}


PE_NAME=${PE_MPICH_FIXED_PRGENV}
SLATEC_INSTALL_PATH=${SLATEC_INSTALL_ROOT}/${PE_NAME}/${PE_RELEASE}/${BUILD}

cd ${SLATEC_BUILD_ROOT}

cp makefile.ARCHER2 makefile
sed -i "s:/usr/local:${SLATEC_INSTALL_PATH}:g" ${SLATEC_BUILD_ROOT}/makefile
sed -i "s:ldconfig:/sbin/ldconfig -C ${SLATEC_INSTALL_PATH}/etc/ld.so.cache:g" ${SLATEC_BUILD_ROOT}/makefile

cp ./static/makefile.sav ./static/makefile
cp ./dynamic/makefile.sav ./dynamic/makefile

set_compile_options ${SLATEC_BUILD_ROOT}/static/makefile
set_compile_options ${SLATEC_BUILD_ROOT}/dynamic/makefile

if [[ "${BUILD}" == "scorep" ]]; then
  sed -i "s:\${FC}:scorep --user \${FC}:g" ${SLATEC_BUILD_ROOT}/static/makefile
  sed -i "s:\${FC}:scorep --user \${FC}:g" ${SLATEC_BUILD_ROOT}/dynamic/makefile
fi

rm -rf ${SLATEC_INSTALL_PATH}
mkdir -p ${SLATEC_INSTALL_PATH}/lib
mkdir -p ${SLATEC_INSTALL_PATH}/etc
mkdir -p ${SLATEC_INSTALL_PATH}/man/man1

FC=ftn make
make install
make clean
