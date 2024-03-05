#!/bin/bash


function set_compile_options {
  MAKEFILE=$1
  if [[ "${BUILD}" == "release" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -O3 -h omp:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -O3 -fopenmp -fallow-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -O3 -fopenmp:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "debug" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O0 -h omp:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O0 -fopenmp -fallow-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8 -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O0 -fopenmp:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "craypat" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O3 -DCRAYPAT -h omp -h profile_generate:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O3 -fopenmp -DCRAYPAT -fallow-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS = -g -O3 -fopenmp -DCRAYPAT:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "armmap" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -G2 -O3 -h omp -h ipa0:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g1 -O3 -fopenmp -fallow-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8 -fno-inline -fno-optimize-sibling-calls:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g1 -O3 -fopenmp:g" ${MAKEFILE}
    fi
  elif [[ "${BUILD}" == "scorep" ]]; then
    if [[ "${PRGENV}" == "cray" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -G2 -O3 -h omp -h ipa0:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "gnu" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g1 -O3 -fopenmp -fallow-argument-mismatch -std=legacy -fdefault-real-8 -fdefault-double-8 -fno-inline -fno-optimize-sibling-calls:g" ${MAKEFILE}
    elif [[ "${PRGENV}" == "aocc" ]]; then
      sed -i "s:FFLAGS =:FFLAGS= -g1 -O3 -fopenmp:g" ${MAKEFILE}
    fi
  fi
}


PE_RELEASE=$1
PRGENV=$2
BUILD=$3
VERSION=5.0
GLOBLIBI_VERSION=5.0
SLATEC_VERSION=4.1
ERRMSG="Invalid syntax: build.sh <CPE release> cray|gnu|aocc release|debug|craypat|armmap|scorep"


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
WMAM_LABEL=WMAM
WMAM_VERSION=${VERSION}
WMAM_NAME=${WMAM_LABEL}-${WMAM_VERSION}
WMAM_BUILD_ROOT=${PRFX}/projects/eCSE04-8/apps/${WMAM_LABEL}/src
WMAM_INSTALL_ROOT=${PRFX}/apps/${WMAM_LABEL}/${WMAM_VERSION}


echo -e "\n\nBuilding ${WMAM_LABEL} ${WMAM_VERSION} with globlibi ${GLOBLIBI_VERSION} (${BUILD}) using ${PRGENV} (CPE ${PE_RELEASE}) programming environment...\n\n"
  
module -q restore
module -q load cpe/${PE_RELEASE}
module -q load PrgEnv-${PRGENV}

if [[ "${BUILD}" == "craypat" ]]; then
  module -q load perftools-base
  module -q load perftools
elif [[ "${BUILD}" == "scorep" ]]; then
  module -q use /work/y23/shared/scalasca/modulefiles
  if [[ "${PRGENV}" == "cray" ]]; then
    module -q load scalasca/2.6.1-cray
  elif [[ "${PRGENV}" == "gnu" ]]; then
    module -q load scalasca/2.6.1-gcc11
  elif [[ "${PRGENV}" == "aocc" ]]; then
    module -q load scalasca/2.6.1-aocc
  fi
elif [[ "${BUILD}" == "armmap" ]]; then
  module -q load arm/forge
fi

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

PE_NAME=${PE_MPICH_FIXED_PRGENV}
WMAM_INSTALL_PATH=${WMAM_INSTALL_ROOT}/${PE_NAME}/${PE_RELEASE}/${BUILD}

SLATEC_ROOT=${PRFX}/libs/slatec/${SLATEC_VERSION}/${PE_NAME}/${PE_RELEASE}/${BUILD}
GLOBLIBI_ROOT=${PRFX}/libs/globlibi/${GLOBLIBI_VERSION}/${PE_NAME}/${PE_RELEASE}/${BUILD}

cd ${WMAM_BUILD_ROOT}

cp makefile.ARCHER2 makefile
set_compile_options ./makefile

LIBS_MAKEFILE_LINE="${GLOBLIBI_ROOT}/lib/libgloblibi.a ${SLATEC_ROOT}/lib/libslatec.a"
if [[ "${BUILD}" == "armmap" ]]; then
  ARM_MAPLIB_PATH=${FORGE_DIR}/map/libs/cpe-${PE_RELEASE}/${PRGENV}/ofi
  LIBS_MAKEFILE_LINE="${LIBS_MAKEFILE_LINE} -L${ARM_MAPLIB_PATH} -lmap-sampler-pmpi -lmap-sampler -Wl,--eh-frame-hdr -Wl,-rpath=${ARM_MAPLIB_PATH}"
fi
sed -i "s:LIBS =:LIBS = ${LIBS_MAKEFILE_LINE}:g" ./makefile

if [[ "${BUILD}" == "scorep" ]]; then
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
