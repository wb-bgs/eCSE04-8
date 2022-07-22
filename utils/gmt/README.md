Generic Mapping Tools (GMT)
===========================

GMT is an open source collection of about 90 command-line tools
for manipulating geographic and Cartesian data sets.

GMT is used to help generate the input data for the various WMAM
test cases.

The GMT source is stored in this GitLab repo for convenience.
The [build.sh](build.sh) script builds the GMT suite for [ARCHER2](https://www.archer2.ac.uk/).


The instructions below show how to obtain the GMT source
code from [github.com](https://github.com/GenericMappingTools/gmt).

```bash
PRFX=${HOME/home/work}/eCSE04-8/utils
GMT_LABEL=gmt
GMT_VERSION=6.3.0
GMT_NAME=${GMT_LABEL}-${GMT_VERSION}

mkdir -p ${PRFX}/${GMT_LABEL}
cd ${PRFX}/${GMT_LABEL}

wget https://github.com/GenericMappingTools/${GMT_LABEL}/releases/download/${GMT_VERSION}/${GMT_NAME}-src.tar.gz
tar -xzf ${GMT_NAME}-src.tar.gz
rm ${GMT_NAME}-src.tar.gz
mv ${GMT_NAME} src
```