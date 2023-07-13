SLATEC
======

The SLATEC Common Mathematical Library, Version 4.1, July 1993
is a comprehensive software library containing over 1400 general
purpose mathematical and statistical routines written in Fortran 77.

SLATEC is required by the World Magnetic Anomaly Model (WMAM) code
and by the iterative version of the global model fitting library (GlobLibI).

The SLATEC source is stored in this GitLab repo (alongside WMAM and
GlobLibI) for convenience. The [build.sh](build.sh) script builds
the SLATEC library for [ARCHER2](https://www.archer2.ac.uk/).

Simply run `./build.sh <CPE release> <compiler environment> <build type>` to build the SLATEC library.

As of July 2023, the Cray Programming Environment (CPE) release installed on ARCHER2 is `22.12`.
The available CPE version may change in the future - this can be checked by running `module avail cpe`.
If the chosen CPE release is not installed on ARCHER2, the default CPE release will be used instead.

The compiler environment can be `cray`, `gnu` or `aocc`; these are the three compiler environments
provided by the CPE release.

The build type can be `release`, `debug`, `craypat`, `armmap` or `scalasca`; the last three
build types refer to profiling tools, [Cray PAT](https://docs.archer2.ac.uk/user-guide/profile/#craypat), \
[Arm MAP](https://docs.archer2.ac.uk/data-tools/arm-forge/) and [Scalasca](https://www.archer2.ac.uk/training/courses/230822-scalasca/).

Once the build has completed the library file(s) will be written to
`${HOME/home/work}/libs/slatec/`.


The instructions below show how to obtain the SLATEC source
code from [netlib.org](http://www.netlib.org/slatec/).

```bash
PRFX=${HOME/home/work}/eCSE04-8/libs
SLATEC_LABEL=slatec
SLATEC_VERSION=4.1
SLATEC_VERSION_MAJOR=`echo ${SLATEC_VERSION} | cut -d'.' -f1`
SLATEC_NAME=${SLATEC_LABEL}-${SLATEC_VERSION}
SLATEC_ARC_MAKE=${SLATEC_LABEL}${SLATEC_VERSION_MAJOR}linux.tgz
SLATEC_ARC_SRC=${SLATEC_LABEL}_src.tgz

mkdir -p ${PRFX}/${SLATEC_LABEL}/src
cd ${PRFX}/${SLATEC_LABEL}/src

wget http://www.netlib.org/${SLATEC_LABEL}/${SLATEC_ARC_MAKE}
tar -xvzf ${SLATEC_ARC_MAKE}
rm ${SLATEC_ARC_MAKE}

cp makefile makefile.ARCHER2
cp ./static/makefile ./static/makefile.sav
cp ./dynamic/makefile ./dynamic/makefile.sav

wget http://www.netlib.org/${SLATEC_LABEL}/${SLATEC_ARC_SRC}
tar -xvzf ${SLATEC_ARC_SRC}
rm ${SLATEC_ARC_SRC}
mv ./src/* ./
rmdir src
```
