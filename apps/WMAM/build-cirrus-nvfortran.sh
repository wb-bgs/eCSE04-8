#!/bin/bash


function set_compile_options {
  MAKEFILE=$1
  if [[ "${BUILD}" == "release" ]]; then
    sed -i "s:FFLAGS =:FFLAGS = -I${MPI_HOME}/include -O3 -cpp -r8 -tp=cascadelake -target=gpu -mp=gpu -Minfo=mp -gpu=cuda12.4,cc70,mem\:unified -cuda:g" ${MAKEFILE}
  elif [[ "${BUILD}" == "debug" ]]; then
    sed -i "s:FFLAGS =:FFLAGS = -I${MPI_HOME}/include -g -O0 -cpp -r8 -tp=cascadelake -target=gpu -mp=gpu -Minfo=mp -gpu=cuda12.4,cc70 -cuda -C -Mnobounds -ffpe-trap=invalid,zero,overflow,underflow,inexact -Ktrap=divz,denorm,inexact,inv,ovf,unf:g" ${MAKEFILE}
  fi
  sed -i "s:LIBS =:LIBS = -L${MPI_HOME}/lib -lmpi_mpifh:g" ${MAKEFILE}
}


BUILD=$1
VERSION=5.0.cuda
SLATEC_VERSION=4.1
GCC_VERSION=10.2.0
NVHPC_VERSION=24.5
ERRMSG="Invalid syntax: build-cirrus-nvfortran.sh release|debug"


if [[ "${BUILD}" != "release" && "${BUILD}" != "debug" ]]; then
  echo ${ERRMSG}
  exit
fi


module -s load nvidia/nvhpc-nompi/${NVHPC_VERSION}
module -s load openmpi/4.1.6-cuda-12.4-nvfortran


PRFX=${HOME/home/work}
WMAM_LABEL=WMAM
WMAM_VERSION=${VERSION}
WMAM_NAME=${WMAM_LABEL}-${WMAM_VERSION}
WMAM_BUILD_ROOT=${PRFX}/projects/eCSE04-8/apps/${WMAM_LABEL}/src
WMAM_INSTALL_ROOT=${PRFX}/apps/${WMAM_LABEL}/${WMAM_VERSION}/GNU/${GCC_VERSION}-nvfortran
WMAM_INSTALL_PATH=${WMAM_INSTALL_ROOT}/${BUILD}

echo -e "\n\nBuilding ${WMAM_LABEL} ${WMAM_VERSION} (${BUILD}) using GNU ${GCC_VERSION}-nvfortran...\n\n"
  
SLATEC_BUILD=${BUILD}
SLATEC_ROOT=${PRFX}/libs/slatec/${SLATEC_VERSION}/GNU/${GCC_VERSION}-nvfortran/${SLATEC_BUILD}

cd ${WMAM_BUILD_ROOT}

cp makefile.cirrus.nvfortran makefile
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
