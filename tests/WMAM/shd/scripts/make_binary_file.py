#!/usr/bin/env python 

# python make_binary_file.py /work/ecsead08/ecsead08/mrbad08/tests/WMAM/shd/300 coef_1990_15.dat
# python make_binary_file.py /work/ecsead08/ecsead08/mrbad08/tests/WMAM/shd/300 model.in
# python make_binary_file.py /work/ecsead08/ecsead08/mrbad08/tests/WMAM/shd/300 wdmam_geocentric.dat

import sys
import os
import struct
from decimal import Decimal


def convertRefModelFile(fpath):
  with open(fpath, 'r') as fa, open(fpath+".bin", 'wb') as fb:
    while True:
      ln = fa.readline()
      if not ln or 0 == len(ln):
        break    

      if ' ' == ln[0]:
        # convert reference year
        fb.write(struct.pack('d',float(ln)))
      elif 'x' == ln[0]:
        # convert coefficient
        lnp = ln.split()
        fb.write(struct.pack('d',float(lnp[3])))
      else:
        continue

  return


def convertModelFile(fpath):
  with open(fpath, 'r') as fa, open(fpath+".bin", 'wb') as fb:
    while True:
      ln = fa.readline()
      if not ln or 0 == len(ln):
        break

      if ' ' == ln[0]:
        lnp = ln.split()
        if 3 != len(lnp):
          continue
        # convert coefficent 
        fb.write(struct.pack('d',float(lnp[1])))
      else:
        continue

  return


def convertDataFile(fpath):
  with open(fpath, 'r') as fa, open(fpath+".bin", 'wb') as fb:
    while True:
      ln = fa.readline()
      if not ln or 0 == len(ln):
        break

      lnp = ln.split()
      
      # convert longitude
      fb.write(struct.pack('d',float(lnp[0])))
      # convert colatitude
      fb.write(struct.pack('d',float(lnp[1])))
      # convert radius
      fb.write(struct.pack('d',float(lnp[2])))
      # convert magscalar
      fb.write(struct.pack('d',float(lnp[3])))

  return


workdir = sys.argv[1]
ascii_file = sys.argv[2]
fpath = workdir+'/'+ascii_file
if ascii_file == "coef_1990_15.dat":
  convertRefModelFile(fpath)
elif ascii_file == "model.in":
  convertModelFile(fpath)
elif ascii_file == "wdmam_geocentric.dat":
  convertDataFile(fpath)
else:
  print("Unrecognized file: ", ascii_file)  
