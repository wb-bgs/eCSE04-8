#!/bin/bash
  
module -q load cray-python


PYTHON_ROOT=${HOME/home/work}/utils/pypp
PYTHON_VER=`echo ${CRAY_PYTHON_LEVEL} | cut -d'.' -f1-2`
PYTHON_BIN=${PYTHON_ROOT}/${CRAY_PYTHON_LEVEL}/bin
PYTHON_ACTIVATE=${PYTHON_BIN}/activate
PYTHON_DEACTIVATE=${PYTHON_BIN}/deactivate


mkdir -p ${PYTHON_BIN}

echo -e '#!/bin/bash\n' > ${PYTHON_ACTIVATE}
echo -e "module load cray-python/${CRAY_PYTHON_LEVEL}\n" >> ${PYTHON_ACTIVATE}
echo -e "PYTHON_ROOT=${HOME/home/work}/utils/pypp\n" >> ${PYTHON_ACTIVATE}
echo -e "export PIP_CACHE_DIR=\${PYTHON_ROOT}/.cache/pip" >> ${PYTHON_ACTIVATE}
echo -e "export MPLCONFIGDIR=\${PYTHON_ROOT}/.config/matplotlib\n" >> ${PYTHON_ACTIVATE}
echo -e "export PYTHONUSERBASE=\${PYTHON_ROOT}/${CRAY_PYTHON_LEVEL}" >> ${PYTHON_ACTIVATE}
echo -e "export PATH=\${PYTHONUSERBASE}/bin:${PATH}" >> ${PYTHON_ACTIVATE}
echo -e "export PYTHONPATH=\${PYTHONUSERBASE}/lib/python${PYTHONVER}/site-packages:${PYTHONPATH}" >> ${PYTHON_ACTIVATE}

echo -e '#!/bin/bash\n' > ${PYTHON_DEACTIVATE}
echo -e "export PIP_CACHE_DIR=${PIP_CACHE_DIR}" >> ${PYTHON_DEACTIVATE}
echo -e "export MPLCONFIGDIR=${MPLCONFIGDIR}\n" >> ${PYTHON_DEACTIVATE}
echo -e "export PATH=${PATH}" >> ${PYTHON_DEACTIVATE}
echo -e "export PYTHONPATH=${PYTHONPATH}\n" >> ${PYTHON_DEACTIVATE}
echo -e "module unload cray-python" >> ${PYTHON_DEACTIVATE}

chmod 700 ${PYTHON_ACTIVATE}
chmod 700 ${PYTHON_DEACTIVATE}


. ${PYTHON_ACTIVATE}

pip install --user --upgrade pip
pip install --user numpy
pip install --user scipy
pip install --user matplotlib
pip install --user pyqt5

. ${PYTHON_DEACTIVATE}
