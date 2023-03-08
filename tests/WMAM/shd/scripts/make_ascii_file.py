#!/usr/bin/env python 

# python make_binary_file.py /path/to/results fit_No_P.out.bin
# python make_binary_file.py /path/to/results fit_damp.out.bin

import sys
import os
import struct


def fortranFormat(d):
  sn = '{:1.8E}'.format(d)
  expi = sn.find('E')
  
  msign = sn[0]
  if msign == '-':
    man = sn[1]+sn[3:expi]
  else:
    msign = '+'
    man = sn[0]+sn[2:expi]
  man = '0.' + man

  esign = sn[expi+1]
  exp = (int(sn[expi+1:]) if esign == '-' else int(sn[expi+2:]))
  exp += 1
  if exp < 0:
    exp *= -1
  exp = '{:02d}'.format(exp)
  exp = 'E' + esign + exp
  
  spacing = (' ' if msign == '-' else '  ')
  msign = ('-' if msign == '-' else '')

  return spacing + msign + man + exp


def convertFitDataFile(fpath):
  with open(fpath, 'rb') as fb, open(fpath[:-4], 'w') as fa:
    fa.write('##########\n')
    while True:
      # assume ten 8-byte reals per line
      bb = fb.read(10*8)
      if not bb:
        break

      vals = struct.unpack('dddddddddd', bb)
      if not vals or 10 != len(vals):
        break

      ln = ''
      for v in vals:
        ln += fortranFormat(v)

      fa.write(ln+'\n')
  return


workdir = sys.argv[1]
bin_file = sys.argv[2]
fpath = workdir+'/'+bin_file
if bin_file == "fit_No_P.out.bin" or bin_file == "fit_damp.out.bin" :
  print('Converting ' + bin_file + ' to ASCII format...')
  convertFitDataFile(fpath)
  print('Done')
else:
  print("Unrecognized file: ", bin_file)