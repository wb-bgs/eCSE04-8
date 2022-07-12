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
      sed -i "s:FFLAGS =:FFLAGS = -g -O0 :g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O0 -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O0:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "craypat" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O3 -h profile_generate:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O3 -fallow-argument-mismatch:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O3:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "armmap" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -G2 -O3 -h ipa0:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g1 -O3 -fno-inline -fno-optimize-sibling-calls:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g1 -O3:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "scalasca" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -G2 -O3 -h ipa0:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g1 -O3 -fno-inline -fno-optimize-sibling-calls:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g1 -O3:g" ${MAKEFILE}
    fi
  fi
}


PE_RELEASE=21.09
PRGENV=$1
BUILD=$2
VERSION=2.0
GLOBLIBI_VERSION=2.0
SLATEC_VERSION=4.1
ERRMSG="Invalid syntax: build.sh cray|gnu|aocc release|debug|craypat|armmap|scalasca <version> <globlibi version>"

if [[ "${PRGENV}" != "cray" && "${PRGENV}" != "gnu" && "${PRGENV}" != "aocc" ]]; then
  echo ${ERRMSG}
  exit
fi

if [[ "${BUILD}" != "release" && "${BUILD}" != "debug" && "${BUILD}" != "craypat" && "${BUILD}" != "armmap" && "${BUILD}" != "scalasca" ]]; then
  echo ${ERRMSG}
  exit
fi

if [[ "${VERSION}" == "" ]]; then
  echo ${ERRMSG}
  exit
fi

PRFX=${HOME/home/work}
WMAM_LABEL=WMAM
WMAM_VERSION=${VERSION}
WMAM_NAME=${WMAM_LABEL}-${WMAM_VERSION}
WMAM_BUILD_ROOT=${PRFX}/eCSE04-8/apps/${WMAM_LABEL}/src
WMAM_INSTALL_ROOT=${PRFX}/apps/${WMAM_LABEL}/${WMAM_VERSION}


echo -e "\n\nBuilding ${WMAM_LABEL} ${WMAM_VERSION} with globlibi ${GLOBLIBI_VERSION} (${BUILD}) using ${PRGENV} programming environment...\n\n"
  
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
elif [[ "${BUILD}" == "armmap" ]]; then
  module -q load arm/forge/22.0.2  
fi

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

PE_NAME=${PE_MPICH_FIXED_PRGENV}
PE_VERSION=$(eval echo "\${PE_MPICH_GENCOMPILERS_${PE_NAME}}")
WMAM_INSTALL_PATH=${WMAM_INSTALL_ROOT}/${PE_NAME}/${PE_VERSION}/${BUILD}

SLATEC_ROOT=${PRFX}/libs/slatec/${SLATEC_VERSION}/${PE_NAME}/${PE_VERSION}/${BUILD}
GLOBLIBI_ROOT=${PRFX}/libs/globlibi/${GLOBLIBI_VERSION}/${PE_NAME}/${PE_VERSION}/${BUILD}

cd ${WMAM_BUILD_ROOT}

cp makefile.ARCHER2 makefile
set_compile_options ./makefile

LIBS_MAKEFILE_LINE="${GLOBLIBI_ROOT}/lib/libgloblibi.a ${SLATEC_ROOT}/lib/libslatec.a"
if [[ "${BUILD}" == "armmap" ]]; then
  ARM_MAPLIB_PATH=${FORGE_ROOT}/map/lib/${PE_NAME,,}/${PE_VERSION}
  LIBS_MAKEFILE_LINE="${LIBS_MAKEFILE_LINE} -L${ARM_MAPLIB_PATH} -lmap-sampler-pmpi -lmap-sampler -Wl,--eh-frame-hdr -Wl,-rpath=${ARM_MAPLIB_PATH}"
fi
sed -i "s:LIBS =:LIBS = ${LIBS_MAKEFILE_LINE}:g" ./makefile

if [[ "${BUILD}" == "scalasca" ]]; then
  sed -i "s:FC = ftn:FC = scorep --user ftn:g" ./makefile
fi

rm -f ${WMAM_BUILD_ROOT}/*.o
rm -f ${WMAM_BUILD_ROOT}/*.mod
make

mkdir -p ${WMAM_INSTALL_PATH}/bin

if [[ "${BUILD}" == "craypat" ]]; then
  rm -f ${WMAM_INSTALL_PATH}/bin/mod_wmam_020+pat
  pat_build -o ${WMAM_INSTALL_PATH}/bin/mod_wmam_020+pat ${WMAM_BUILD_ROOT}/mod_wmam_020
else
  mv ${WMAM_BUILD_ROOT}/mod_wmam_020 ${WMAM_INSTALL_PATH}/bin/
fi

if [[ "${BUILD}" != "craypat" ]]; then
  rm -f ${WMAM_BUILD_ROOT}/*.o
  rm -f ${WMAM_BUILD_ROOT}/*.mod
fi
