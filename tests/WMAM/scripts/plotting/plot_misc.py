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
npts = np.add(ndatpts, nsampts)
npts_growth = npts / npts[0]
lognpts = np.log10(npts)

samptfrac = 100.0*(nsampts/npts)

ncoeffs = np.array([shd*(shd+2) for shd in shdeg])
ncoeffs_growth = ncoeffs / ncoeffs[0]
logncoeffs = np.log10(ncoeffs)


shdeg_100 = [int(shd/100) if shd/100 % 1 == 0 else shd/100 for shd in shdeg]


plt.plot(shdeg_100, lognpts, linestyle='-', marker='o', color='#66c2a5', label='$N_p$')
plt.plot(shdeg_100, logncoeffs, linestyle='-', marker='*', markersize=9, color='#fc8d62', label='$N_c$')

plt.xlabel('Degree / 100')
plt.xticks(shdeg_100, shdeg_100)
plt.ylabel('Log10')
plt.title(r'WMAM Point ($N_p$) and Coefficient ($N_c$) Counts')
plt.legend(loc='upper left')
plt.show()
fig.savefig('wmam-npts-ncoeffs.eps', format='eps', dpi=1000)
plt.show()
plt.clf()
fig = plt.gcf()


plt.plot(shdeg_100, npts_growth, linestyle='-', marker='o', color='#66c2a5', label='$N_p$ / 40400')
plt.plot(shdeg_100, ncoeffs_growth, linestyle='-', marker='*', markersize=9, color='#fc8d62', label='$N_c$ / 145041')

plt.xlabel('Degree / 100')
plt.xticks(shdeg_100, shdeg_100)
plt.ylabel('')
plt.title(r'WMAM Point and Coefficient Growth from Degree 200')
plt.legend(loc='upper left')
plt.show()

fig.savefig('wmam-npts-ncoeffs-growth.eps', format='eps', dpi=1000)
plt.show()
plt.clf()
fig = plt.gcf()


plt.plot(shdeg_100,samptfrac,linestyle='-',marker='o',color='#fc8d62')
plt.xlabel('Degree / 100 ')
plt.xticks(shdeg_100, shdeg_100)
plt.ylabel('Percentage')
plt.title('WMAM Points Sampled Percentage')
plt.show()

fig.savefig('wmam-sampt-frac.eps', format='eps', dpi=1000)
plt.show()
plt.clf()
fig = plt.gcf()

memreq1 = [39.123, 39.524, 42.956, 55.407, 252.025]
memreq2 = [0.248, 0.56, 3.197, 12.858, 25.297]

plt.plot(shdeg_100,memreq1,linestyle='-',marker='o',label='WMAM v1.2',color='#fc8d62')
plt.plot(shdeg_100,memreq2,linestyle='-',marker='P',label='WMAM v3.7',color='#8da0cb')
plt.xlabel('Degree')
plt.xticks(shdeg_100, shdeg_100)
plt.ylabel('Percentage')
plt.title('WMAM Memory Requirement as percentage of 256 GB')
plt.legend(loc='upper left')
plt.show()

fig.savefig('wmam-mem-req.eps', format='eps', dpi=1000)
plt.show()
plt.clf()
fig = plt.gcf()
