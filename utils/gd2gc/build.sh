#!/bin/bash

PRFX=${HOME/home/work}

GD2GC_LABEL=gd2gc
GD2GC_VERSION=1.0
GD2GC_NAME=${GD2GC_LABEL}-${GD2GC_VERSION}
GD2GC_BUILD_ROOT=${PRFX}/eCSE04-8/utils/${GD2GC_LABEL}/src
GD2GC_INSTALL_ROOT=${PRFX}/utils/${GD2GC_LABEL}/${GD2GC_VERSION}


echo -e "\n\nBuilding ${GD2GC_LABEL} ${GD2GC_VERSION} using the gnu programming environment...\n\n"
  
module -q restore
module -q load PrgEnv-gnu


cd ${GD2GC_BUILD_ROOT}

make

mkdir -p ${GD2GC_INSTALL_ROOT}/bin
mv gd2gc_cmdln ${GD2GC_INSTALL_ROOT}/bin/
rm *.o
