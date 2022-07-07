#!/bin/bash


function set_compile_options {
  MAKEFILE=$1
  if [[ "${BUILD}" == "release" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -O3:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -O3 -fallow-argument-mismatch:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -O1:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "debug" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g -O0:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g -O0 -fallow-argument-mismatch -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g -O0:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "craypat" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g -O3 -h profile_generate:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g -O3 -fallow-argument-mismatch -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g -O1:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "armmap" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -G2 -O3 -h ipa0:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g1 -O3 -fallow-argument-mismatch -fno-inline -fno-optimize-sibling-calls:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g1 -O1:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "scalasca" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -G2 -O3 -h ipa0:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g1 -O3 -fallow-argument-mismatch -fno-inline -fno-optimize-sibling-calls:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS= -O2:FFLAGS= -g1 -O1:g" ${MAKEFILE}
    fi
  fi
}


PE_RELEASE=21.09
PRGENV=$1
BUILD=$2
ERRMSG="Invalid syntax: build.sh cray|gnu|aocc release|debug|craypat|armmap|scalasca"

if [[ "${PRGENV}" != "cray" && "${PRGENV}" != "gnu" && "${PRGENV}" != "aocc" ]]; then
  echo ${ERRMSG}
  exit
fi

if [[ "${BUILD}" != "release" && "${BUILD}" != "debug" && "${BUILD}" != "craypat" && "${BUILD}" != "armmap" && "${BUILD}" != "scalasca" ]]; then
  echo ${ERRMSG}
  exit
fi


PRFX=${HOME/home/work}
SLATEC_LABEL=slatec
SLATEC_VERSION=4.1
SLATEC_BUILD_ROOT=${PRFX}/eCSE04-8/libs/${SLATEC_LABEL}/src
SLATEC_INSTALL_ROOT=${PRFX}/libs/${SLATEC_LABEL}/${SLATEC_VERSION}


echo -e "\n\nBuilding ${SLATEC_LABEL} ${SLATEC_VERSION} (${BUILD}) using ${PRGENV} programming environment...\n\n"
  
module -q restore
module -q load cpe/${PE_RELEASE}
module -q load PrgEnv-${PRGENV}

if [[ "${BUILD}" == "craypat" ]]; then
  module -q load perftools-base
  module -q load perftools
elif [[ "${BUILD}" == "scalasca" ]]; then
  module -q use /work/y23/shared/scalasca/modulefiles
  if [[ "${PRGENV}" == "cray" ]]; then
    module -q load scalasca/2.6-cce
  elif [[ "${PRGENV}" == "gnu" ]]; then
    module -q load scalasca/2.6-gcc10
  else
    echo "Error, ${PRGENV} not supported by scalasca, please try either cray or gnu."
    exit
  fi
fi

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}


PE_NAME=${PE_MPICH_FIXED_PRGENV}
PE_VERSION=$(eval echo "\${PE_MPICH_GENCOMPILERS_${PE_NAME}}")
SLATEC_INSTALL_PATH=${SLATEC_INSTALL_ROOT}/${PE_NAME}/${PE_VERSION}/${BUILD}

cd ${SLATEC_BUILD_ROOT}

cp makefile.ARCHER2 makefile
sed -i "s:/usr/local:${SLATEC_INSTALL_PATH}:g" ${SLATEC_BUILD_ROOT}/makefile
sed -i "s:ldconfig:/sbin/ldconfig -C ${SLATEC_INSTALL_PATH}/etc/ld.so.cache:g" ${SLATEC_BUILD_ROOT}/makefile

cp ./static/makefile.sav ./static/makefile
cp ./dynamic/makefile.sav ./dynamic/makefile

set_compile_options ${SLATEC_BUILD_ROOT}/static/makefile
set_compile_options ${SLATEC_BUILD_ROOT}/dynamic/makefile

if [[ "${BUILD}" == "scalasca" ]]; then
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
