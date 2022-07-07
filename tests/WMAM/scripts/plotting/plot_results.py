#!/usr/bin/env python
  
import sys
import os
import math
from datetime import datetime
from decimal import Decimal
import matplotlib.pyplot as plt

markerLabels = [ 'o', '*', 'v', 's', 'd', '^', 'p' ]
markerSizes = [ 7, 10, 10, 10, 10, 10, 10 ]
colorCodes = ['#1b9e77', '#d95f02', '#7570b3']

clArgCnt = len(sys.argv)

resfn = sys.argv[1]
xi = int(sys.argv[2])
xl = sys.argv[3]
yi = int(sys.argv[4])
yl = sys.argv[5]
yf = float(sys.argv[6])
if 1.0 != yf:
  yf = yf/(3600*24*365)
tt = sys.argv[7]
pfn = sys.argv[8]
lw = float(sys.argv[9])
legpos = sys.argv[10]
inv = True if ('1' == sys.argv[11]) else False
cmp = True if ('1' == sys.argv[12]) else False

labelList = sys.argv[13].split('|')
labelOffsets = sys.argv[14].split(',')
labelSpacing = float(labelOffsets[-1])
labelOrigin = [float(labelOffsets[0]), float(labelOffsets[1])]

resData = os.popen('cat ' + resfn).read()
resLines = resData.split('\n')

results = []

for ln in resLines:
  if 0 == len(ln):
    # empty line
    continue
  elif '#' == ln[0]:
    # comment line
    continue
  else:
    # results for machine macName
    lnParts = ln.split()
    resList = []
    for lnpt in lnParts:
      resList.append(float(lnpt))
    results.append(resList)


fig = plt.gcf()

mi = 0
xpts = []
ypts = []
yptsw = []
yptsb = []
for result in results:
  xp = result[xi]
  yp = result[yi]*yf
  #ypw = result[yi+1]*yf
  #ypb = result[yi+2]*yf
  if inv:
    yp = 1.0/yp
    #ypw = 1.0/ypw
    #ypb = 1.0/ypb
  xpts.append(xp)
  ypts.append(yp)
  #yptsw.append(ypw)
  #yptsb.append(ypb)

plt.plot(xpts, ypts, label='', color=colorCodes[mi], marker=markerLabels[mi], markersize=markerSizes[mi], linestyle='-', linewidth=lw)
#plt.plot(xpts, yptsw, label='', color=colorCodes[mi], marker='_', markersize=10, linestyle='')
#plt.plot(xpts, yptsb, label='', color=colorCodes[mi], marker='_', markersize=10, linestyle='')
mi += 1
if mi >= len(markerLabels):
  mi = 0

#plt.legend(loc=legpos)
plt.xlabel(xl)
plt.ylabel(yl)
plt.title(tt)

if len(labelList) > 0:
  for label in labelList:
    labelPos = tuple(map(float, labelOrigin))
    plt.annotate(label, xy=labelPos, xytext=labelPos, xycoords='axes fraction')
    labelOrigin[-1] += labelSpacing

plt.show()

fig.savefig(pfn+'.eps', format='eps', dpi=1000)
