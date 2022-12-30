#!/usr/bin/env python
  
import numpy as np
import matplotlib.pyplot as plt


fig = plt.gcf()

shdeg  = [ 200,  300,  720, 1440, 2000]
resdeg = [1.00, 0.50, 0.25, 0.10, 0.05]

ny = [int(1.0/rd)*360 for rd in resdeg]
nx = [int(1.0/rd)*180-1 for rd in resdeg]

ndatpts = np.multiply(nx, ny)
nsampts = np.array([(shd+1)*(2*shd+1) for shd in shdeg])
ntotpts = np.add(ndatpts, nsampts)
lognpts = np.log10(ntotpts)
samptfrac = 100.0*(nsampts/ntotpts)
ncoeffs = np.array([shd*(shd+2) for shd in shdeg])
logncoeffs = np.log10(ncoeffs)

shdeg_100 = [int(shd/100) if shd/100 % 1 == 0 else shd/100 for shd in shdeg]


plt.plot(shdeg_100,lognpts,linestyle="",marker="o",color="#fc8d62")
plt.xlabel('Spherical Harmonic Degree / 100')
plt.xticks(shdeg_100, shdeg_100)
plt.ylabel('Log10')
plt.title('WMAM Total Point Count')
plt.show()
fig.savefig('wmam-pt-cnt.eps', format='eps', dpi=1000)

plt.clf()

plt.plot(shdeg_100,samptfrac,linestyle="",marker="o",color="#fc8d62")
plt.xlabel('Spherical Harmonic Degree')
plt.xticks(shdeg_100, shdeg_100)
plt.ylabel('Percentage')
plt.title('Percentage of WMAM Points Sampled')
plt.show()
fig.savefig('wmam-sampt-frac.eps', format='eps', dpi=1000)

plt.clf()

plt.plot(shdeg_100,logncoeffs,linestyle="",marker="o",color="#fc8d62")
plt.xlabel('Spherical Harmonic Degree')
plt.xticks(shdeg_100, shdeg_100)
plt.ylabel('Log10')
plt.title('WMAM Number of Coefficients')
plt.show()
fig.savefig('wmam-coeff-cnt.eps', format='eps', dpi=1000)

plt.clf()

memreq1 = [0.9, 2.6, 12.7, 62.3, 187.7]
memreq2 = [0.2, 0.5, 3.1, 12.4, 24.3]
plt.plot(shdeg_100,memreq1,linestyle="",marker="o",label="<= WMAM v2.1",color="#fc8d62")
plt.plot(shdeg_100,memreq2,linestyle="",marker="P",label=">= WMAM v2.2",color="#8da0cb")
plt.xlabel('Spherical Harmonic Degree')
plt.xticks(shdeg_100, shdeg_100)
plt.ylabel('Percentage')
plt.title('WMAM ARCHER2 Node Memory Requirement')
plt.legend(loc="upper left")
plt.show()
fig.savefig('wmam-a2mem-req.eps', format='eps', dpi=1000)
