#!/usr/bin/env python 

import sys
import os
import math
import numpy
from datetime import datetime
from decimal import Decimal
from scipy.integrate import quad
from pylab import *


WMAM_NAME = 'WMAM'
ARCHER2_CORES_PER_NODE = 128

# the standard deviation limit is expressed
# as a percentage of the mean
STD_DEV_LIMIT = 10.0

debug = 0


def update_running_stat_dict(rsDict, x):
  n = rsDict['n']
  min = rsDict['min']
  max = rsDict['max']
  m = rsDict['m']

  n1 = n
  t = [0.0, 0.0, 0.0, 0.0]

  n += 1.0;

  if 1.0 < n:
    if x < min: min = x
    if x > max: max = x
  else:
    min = max = x

  t[0] = x - m[0]
  t[1] = t[0] / n
  t[2] = t[1] * t[1]
  t[3] = t[0] * t[1] * n1

  m[0] += t[1]
  m[3] += t[3]*t[2]*(n*n - 3.0*n + 3.0) + 6.0*t[2]*m[1] - 4.0*t[1]*m[2]
  m[2] += t[3]*t[1]*(n-2.0) - 3.0*t[1]*m[1]
  m[1] += t[3]

  rsDict['n'] = n
  rsDict['min'] = min
  rsDict['max'] = max
  rsDict['m'] = m


def get_count(rsDict):
  return int(rsDict['n'])

def get_mean(rsDict):
  return rsDict['m'][0]

def get_stddev(rsDict):
  n = rsDict['n']
  m = rsDict['m']
  return math.sqrt(m[1]/(n-1.0) if n > 1.0 else 0.0)

def get_width(rsDict):
  mean = abs(get_mean(rsDict))
  stdev = get_stddev(rsDict)
  return stdev/mean*100.0

def get_skewness(rsDict):
  n = rsDict['n']
  m = rsDict['m']
  return math.sqrt(n)*m[2]/math.pow(m[1],1.5) if 0.0 != m[1] else 0.0

def get_kurtosis(rsDict):
  n = rsDict['n']
  m = rsDict['m']
  return (n*m[3])/(m[1]*m[1]) - 3.0 if 0.0 != m[1] else 0.0

def is_stddev_within_limit(rsDict):
  mean = abs(get_mean(rsDict))
  stddev = get_stddev(rsDict)
  width = stddev/mean*100.0
  return width <= STD_DEV_LIMIT

def serialise_running_stat_dict(rsDict):
  cnt = get_count(rsDict)
  mean = get_mean(rsDict)
  stddev = get_stddev(rsDict)
  min = rsDict['min']
  max = rsDict['max']
  width = get_width(rsDict)
  skew = get_skewness(rsDict)
  kurt = get_kurtosis(rsDict)

  rsStr = str(cnt) + ' '
  rsStr += '%.3E' % Decimal(mean) + ' '
  rsStr += '%.3E' % Decimal(stddev) + ' '
  rsStr += str(round(width,2)) + '% '
  rsStr += '[' + '%.3E' % Decimal(min) + ','
  rsStr += '%.3E' % Decimal(max) + '] '
  rsStr += '%.3E' % Decimal(skew) + ' '
  rsStr += '%.3E' % Decimal(kurt)

  return rsStr



def get_wmam_metric(outFilenames, coresPerNode):
  if 0 < debug: print('Extracting WMAM performance metric(s)...')
  results = {}
  for outfn in outFilenames:
    if 0 < debug: print('Reading '+outfn+'...')
    lfs_i = outfn.rindex('/')
    pfs_i = outfn[:lfs_i].rindex('/')
    ppfs_i = outfn[:pfs_i].rindex('/')

    outLabel = outfn[ppfs_i+1:pfs_i]
    nc = -1
    if 'n' in outLabel:
      nc = int(outLabel[1:])*coresPerNode
    elif 'c' in outLabel:
      nc = int(outLabel[1:])
    outKey = nc

    runTime = -1.0

    outData = os.popen('cat ' + outfn).read()
    outLines = outData.split('\n')
    for ln in outLines:
      lnParts = ln.split()
      if 'srun time: ' in ln:
        runTime = float(lnParts[2])

    if -1.0 == runTime: continue

    if outKey not in results: 
        results[outKey] = []

    results[outKey].append({ 'runtime': runTime })

  sortedResultKeys = sorted(results.keys())
  baseCoreCount = sortedResultKeys[0]
  baseRunTime = results[baseCoreCount][0]['runtime']
  for k in sortedResultKeys:
    rt = results[k][0]['runtime']
    speedup = baseRunTime / rt
    eff = speedup / (k/baseCoreCount)
    print(str(k)+' '+str(round(rt,2))+' '+str(round(speedup,2))+' '+str(round(eff,2)))



codeName = sys.argv[1]
outFilenameMask = sys.argv[2]

outFilenameListing = os.popen('ls ' + outFilenameMask).read()
outFilenames = outFilenameListing.split()

if WMAM_NAME in codeName:
  get_wmam_metric(outFilenames, ARCHER2_CORES_PER_NODE)
else:
  print("Error, "+sys.argv[1]+" does not match recognised code name.")
