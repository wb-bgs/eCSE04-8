#!/bin/bash

PYPP_ROOT=${HOME/home/work}/utils/pypp

module -q load cray-python

python -m venv --system-site-packages ${PYPP_ROOT}

source ${PYPP_ROOT}/bin/activate

python -m pip install --upgrade pip
python -m pip install numpy
python -m pip install scipy
python -m pip install matplotlib
python -m pip install pyqt5

deactivate
