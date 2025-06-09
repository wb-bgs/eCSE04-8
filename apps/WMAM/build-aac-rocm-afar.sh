#!/bin/bash


function set_compile_options {
  MAKEFILE=$1
  if [[ "${BUILD}" == "release" ]]; then
    #sed -i "s:FFLAGS =:FFLAGS = -I${MPI_PATH}/include -O0 -cpp -DOMP_OFFLOAD -DOMP_OFFLOAD_CPTP -DOMP_OFFLOAD_SSQGH -fopenmp --offload-arch=gfx942:g" ${MAKEFILE}
    sed -i "s:FFLAGS =:FFLAGS = -I${MPI_PATH}/include -O3 -cpp -DOMP_OFFLOAD -DOMP_OFFLOAD_CPTP -DOMP_OFFLOAD_SSQGH -fopenmp --offload-arch=gfx942:g" ${MAKEFILE}
  elif [[ "${BUILD}" == "debug" ]]; then
    sed -i "s:FFLAGS =:FFLAGS = -I${MPI_PATH}/include -g -O0 -cpp -fdefault-real-8 -DOMP_OFFLOAD -DOMP_OFFLOAD_CPTP -DOMP_OFFLOAD_SSQGH -fopenmp --offload-arch=gfx942 -C -Mnobounds -ffpe-trap=invalid,zero,overflow,underflow,inexact -Ktrap=divz,denorm,inexact,inv,ovf,unf:g" ${MAKEFILE}
  fi
  #sed -i "s:LIBS =:LIBS = -L${MPI_PATH}/lib -lmpi_mpifh:g" ${MAKEFILE}
}


BUILD=$1
VERSION=5.0
SLATEC_VERSION=4.1
ROCM_VERSION=6.4.0
ROCM_AFAR_VERSION=6.1.0
OMPI_VERSION=5.0.7
ERRMSG="Invalid syntax: build-acc-rocm-afar.sh release|debug"


if [[ "${BUILD}" != "release" && "${BUILD}" != "debug" ]]; then
  echo ${ERRMSG}
  exit
fi


module load rocm/${ROCM_VERSION}
module load amdflang-new/rocm-afar-${ROCM_AFAR_VERSION}
module load openmpi-amdflang-${ROCM_AFAR_VERSION}/${OMPI_VERSION}-ucc1.3.0-ucx1.18.0


PRFX=${HOME}
WMAM_LABEL=WMAM
WMAM_VERSION=${VERSION}
WMAM_NAME=${WMAM_LABEL}-${WMAM_VERSION}
WMAM_BUILD_ROOT=${PRFX}/projects/eCSE04-8/apps/${WMAM_LABEL}/src
WMAM_INSTALL_ROOT=${PRFX}/apps/${WMAM_LABEL}/${WMAM_VERSION}/ROCM-AFAR/${ROCM_AFAR_VERSION}-amdflang
WMAM_INSTALL_PATH=${WMAM_INSTALL_ROOT}/${BUILD}

echo -e "\n\nBuilding ${WMAM_LABEL} ${WMAM_VERSION} (${BUILD}) using ROCM AFAR ${ROCM_AFAR_VERSION}-amdflang...\n\n"
  
SLATEC_BUILD=${BUILD}
SLATEC_ROOT=${PRFX}/libs/slatec/${SLATEC_VERSION}/ROCM-AFAR/${ROCM_AFAR_VERSION}-amdflang/${SLATEC_BUILD}

cd ${WMAM_BUILD_ROOT}

cp makefile.aac.rocm-afar makefile
set_compile_options ./makefile

#LIBS_MAKEFILE_LINE="${SLATEC_ROOT}/lib/libslatec.a"
LIBS_MAKEFILE_LINE="-lflang_rt.hostdevice"
sed -i "s:LIBS =:LIBS = ${LIBS_MAKEFILE_LINE}:g" ./makefile


find ${WMAM_BUILD_ROOT} -name '*.o' -delete
find ${WMAM_BUILD_ROOT} -name '*.mod' -delete
make
find ${WMAM_BUILD_ROOT} -name '*.o' -delete
find ${WMAM_BUILD_ROOT} -name '*.mod' -delete

mkdir -p ${WMAM_INSTALL_PATH}/bin

mv ${WMAM_BUILD_ROOT}/mod_wmam_020 ${WMAM_INSTALL_PATH}/bin/
