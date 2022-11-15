#!/bin/bash

declare -a shdegrees=("200" "300" "720" "1440" "2000" "tests")

for shd in "${shdegrees[@]}"; do

#  find ./${shd}/*/*/*.ll -type f -exec sed -i 's:APP_VERSION=3.0:APP_VERSION=3.1:g' {} \;
#  find ./${shd}/*/*/*.ll -type f -exec sed -i 's:SCHEME=1:SCHEME=1\nDAMPFAC=5.0:g' {} \;
#  find ./${shd}/*/*/*.ll -type f -exec sed -i 's:\${SCHEME}:\${SCHEME} \${DAMPFAC}:g' {} \;
#  find ./${shd}/*.ll -type f -exec sed -i 's:PE_VERSION=.*:PE_RELEASE=21.09:g' {} \;
#  find ./${shd}/*.ll -type f -exec sed -i 's:PE_VERSION:PE_RELEASE:' {} \;
#  find ./${shd}/*.ll -type f -exec sed -i 's:results/\${DEGREE}/\${APP_COMPILER_LABEL}:results/\${PE_RELEASE}/\${APP_COMPILER_LABEL}:g' {} \;
#  find ./${shd}/*.ll -type f -exec sed -i 's:(\${APP_MPI_LABEL}-\${APP_COMPILER_LABEL}):(\${PE_RELEASE},\${APP_COMPILER_LABEL},\${APP_MPI_LABEL}-\${APP_COMMS_LABEL}):g' {} \;

  find ./${shd}/*.ll -type f -exec sed -i 's:results/\${PE_RELEASE}/\${APP_COMPILER_LABEL}:results/\${DEGREE}/\${PE_RELEASE}/\${APP_COMPILER_LABEL}:g' {} \;

done
