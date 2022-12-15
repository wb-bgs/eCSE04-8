#!/usr/bin/env python 

import sys
import os
import math
import numpy
from datetime import datetime
from decimal import Decimal



debug = 1


def output_lines(refLn, outLn):
  print("    Ref: "+refLn)
  print("    Out: "+outLn)
  print("")


def verify_output(refFilename,outFilenames):
  
  if 0 < debug: print('Reading '+refFilename+'...')
  refData = os.popen('cat ' + refFilename).read()
  refLines = refData.split('\n')

  lmax_ref = 0
  ln = refLines[1]
  if '#lmax=' in ln:
    lmax_ref = int(ln.split()[1])
  if 0 == lmax_ref:
    print("Error, illegal spherical harmonic degree value, "+str(lmax_ref)+", in reference file.")

  days_ref = float(refLines[3])

  for outfn in outFilenames:
    if 0 < debug: print('Verifying '+outfn+'...')

    max_pval = 0.0
    max_pval_lnum = -1
    max_perr = 0.0
    max_perr_lnum = -1

    outData = os.popen('cat ' + outfn).read()
    outLines = outData.split('\n')

    lmax = 0
    ln = outLines[1]
    if '#lmax=' in ln:
      lmax = int(ln.split()[1])
    if lmax != lmax_ref:
      print("Error, lmax="+str(lmax)+" is not equal to expected spherical harmonic degree, "+str(lmax_ref)+".")

    days = float(outLines[3])
    if days != days_ref:
      print("Error, days="+str(days)+" is not equal to expected days value, "+str(days_ref)+".")

    outLnCnt = len(outLines)
    refLnCnt = len(refLines)
    if outLnCnt != refLnCnt:
      if outLnCnt-5 == 2*(refLnCnt-5):
        print("Each coefficient result is split over two lines - this can happen if running an aocc-built code.")
        outLines2 = []
        for ln in outLines[4::2]:
          outLines2.append(ln)
        for i,ln in enumerate(outLines[5::2]):
          outLines2[i] = outLines2[i] + " " + ln
        outLines = outLines[:4] + outLines2
      else:
        print("Error, length of output file does not match reference file ("+str(len(outLines))+"!="+str(len(refLines))+").")


    for i,ln in enumerate(outLines[4:]):
      lnParts = ln.split()
      if 0 == len(lnParts):
        continue

      # strip commas
      for j, lp in enumerate(lnParts):
        if lp[-1] == ',':
          lnParts[j] = lp[:-1]

      model_type = lnParts[0]

      if len(lnParts) == 4:
        # order is one field
        orderParts = lnParts[1].split('*')
        sh_degree_order = orderParts[-1] + " " + orderParts[-1]
      else:
        sh_degree_order = lnParts[1] + " " + lnParts[2]

      sh_coeff_value = float(lnParts[-2])
      sh_coeff_error = float(lnParts[-1])


      lnRef = refLines[i+4]
      lnRefParts = lnRef.split()

      # strip commas
      for j, lp in enumerate(lnRefParts):
        if lp[-1] == ',':
          lnRefParts[j] = lp[:-1]

      model_type_ref = lnRefParts[0]

      if len(lnRefParts) == 4:
        # order is one field
        orderRefParts = lnRefParts[1].split('*')
        sh_degree_order_ref = orderRefParts[-1] + " " + orderRefParts[-1]
      else:
        sh_degree_order_ref = lnRefParts[1] + " " + lnRefParts[2]

      sh_coeff_value_ref = float(lnRefParts[-2])
      sh_coeff_error_ref = float(lnRefParts[-1])


      if model_type != model_type_ref:
        print("Error, model type mismatch at line "+str(i+5)+".")
        output_lines(refLines[i+4],outLines[i+4])
        continue

      if sh_degree_order != sh_degree_order_ref:
        print("Error, spherical harmonic degree order mismatch at line "+str(i+5)+".")
        output_lines(refLines[i+4],outLines[i+4])
        continue

 
      if sh_coeff_value_ref == 0.0:
        if sh_coeff_value != 0.0:
          print("Error, spherical harmonic degree value mismatched sign at line "+str(i+5)+".")
          output_lines(refLines[i+4],outLines[i+4])
          continue
        pval = 0.0
      else:
        pval = (sh_coeff_value / sh_coeff_value_ref)
      if pval < 0.0:
        print("Error, spherical harmonic degree value mismatched sign at line "+str(i+5)+".")
        output_lines(refLines[i+4],outLines[i+4])
        continue
      elif pval < 1.0:
        pval = 1.0 - pval
      else:
        pval = pval - 1.0
      if pval > max_pval:
        max_pval = pval
        max_pval_lnum = i+5

      perr = (sh_coeff_error / sh_coeff_error_ref)
      if perr < 0.0:
        print("Error, spherical harmonic degree error mismatched sign at line "+str(i+5)+".")
        output_lines(refLines[i+4],outLines[i+4])
        continue
      elif perr < 1.0:
        perr = 1.0 - perr
      else:
        perr = perr - 1.0
      if perr > max_perr:
        max_perr = perr
        max_perr_lnum = i+5


    if -1 < max_pval_lnum:
      print("Maximum spherical harmonic degree value difference of "+"{:.3e}".format(max_pval)+"% at line "+str(max_pval_lnum)+".")
      #output_lines(refLines[max_pval_lnum-1],outLines[max_pval_lnum-1])

    if -1 < max_perr_lnum:
      print("Maximum spherical harmonic degree error percentage difference of "+"{:.3e}".format(max_perr*100.0)+"% at line "+str(max_perr_lnum)+".")
      #output_lines(refLines[max_perr_lnum-1],outLines[max_perr_lnum-1])


outFilenameMask = sys.argv[1]
refFilename = sys.argv[2]

refData = os.popen('cat ' + refFilename).read()
refLines = refData.split('\n')

outFilenameListing = os.popen('ls ' + outFilenameMask).read()
outFilenames = outFilenameListing.split()

verify_output(refFilename,outFilenames)
