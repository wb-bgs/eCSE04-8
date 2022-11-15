#!/bin/bash


function set_compile_options {
  MAKEFILE=$1
  if [[ "${BUILD}" == "release" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -O3:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -O3 -fallow-argument-mismatch:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -O3:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "debug" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O0:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O0 -fallow-argument-mismatch -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O0:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "craypat" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g -O3 -h profile_generate:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g -O3 -fallow-argument-mismatch -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g -O3:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "armmap" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -G2 -O3 -h ipa0:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g1 -O3 -fallow-argument-mismatch -fno-inline -fno-optimize-sibling-calls:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g1 -O3:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "scorep" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -G2 -O3 -h ipa0:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g1 -O3 -fallow-argument-mismatch -fno-inline -fno-optimize-sibling-calls:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g1 -O3:g" ${MAKEFILE}
    fi
  fi

  if [[ "${PRGENV}" == "gnu" ]]; then
    if [[ "${PE_RELEASE}" == "21.04" ]]; then
      sed -i "s:-fallow-argument-mismatch::g" ${MAKEFILE}
    fi
  fi
}


PE_RELEASE=21.09
PRGENV=$1
BUILD=$2
VERSION=3.5
ERRMSG="Invalid syntax: build.sh cray|gnu|aocc release|debug|craypat|armmap|scorep"

if [[ "${PRGENV}" != "cray" && "${PRGENV}" != "gnu" && "${PRGENV}" != "aocc" ]]; then
  echo ${ERRMSG}
  exit
fi

if [[ "${BUILD}" != "release" && "${BUILD}" != "debug" && "${BUILD}" != "craypat" && "${BUILD}" != "armmap" && "${BUILD}" != "scorep" ]]; then
  echo ${ERRMSG}
  exit
fi

if [[ "${VERSION}" == "" ]]; then
  echo ${ERRMSG}
  exit
fi

PRFX=${HOME/home/work}
GLOBLIB_LABEL=globlibi
GLOBLIB_VERSION=${VERSION}
GLOBLIB_NAME=${GLOBLIB_LABEL}-${GLOBLIB_VERSION}
GLOBLIB_BUILD_ROOT=${PRFX}/eCSE04-8/libs/${GLOBLIB_LABEL}/src
GLOBLIB_INSTALL_ROOT=${PRFX}/libs/${GLOBLIB_LABEL}/${GLOBLIB_VERSION}


echo -e "\n\nBuilding ${GLOBLIB_LABEL} ${GLOBLIB_VERSION} (${BUILD}) using ${PRGENV} (CPE ${PE_RELEASE}) programming environment...\n\n"
  
module -q restore
module -q load cpe/${PE_RELEASE}
module -q load PrgEnv-${PRGENV}

if [[ "${BUILD}" == "craypat" ]]; then
  module -q load perftools-base
  module -q load perftools
elif [[ "${BUILD}" == "scorep" ]]; then
  module -q use /work/y23/shared/scalasca/modulefiles
  if [[ "${PRGENV}" == "cray" ]]; then
    module -q load scalasca/2.6-cce
  elif [[ "${PRGENV}" == "gnu" ]]; then
    module -q load scalasca/2.6-gcc10
  else
    echo "Error, ${PRGENV} not supported by scorep, please try either cray or gnu."
    exit
  fi
fi

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}


PE_NAME=${PE_MPICH_FIXED_PRGENV}
GLOBLIB_INSTALL_PATH=${GLOBLIB_INSTALL_ROOT}/${PE_NAME}/${PE_RELEASE}/${BUILD}/lib


cd ${GLOBLIB_BUILD_ROOT}

cp makefile.ARCHER2 makefile
sed -i "s:libdir =:libdir = ${GLOBLIB_INSTALL_PATH}:g" ./makefile

set_compile_options ./makefile

if [[ "${BUILD}" == "scorep" ]]; then
  sed -i "s:\$(FC):scorep --user \$(FC):g" ./makefile
fi

make
make install
make clean
