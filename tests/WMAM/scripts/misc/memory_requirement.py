#!/usr/bin/env python 

# python memory_requirement.py 1 256 128 200  1.00
# 0.231%
# python memory_requirement.py 1 256 128 300  0.50
# 0.522%
# python memory_requirement.py 1 256 128 720  0.25
# 2.978%
# python memory_requirement.py 1 256 128 1440 0.10
# 11.975%
# python memory_requirement.py 1 256 128 2000 0.05
# 23.559%


import sys
import os
import math
import numpy
from datetime import datetime
from decimal import Decimal


KiB = 1024
MiB = 1024*KiB
GiB = 1024*MiB

SIZEOF_DOUBLE=8
SIZEOF_INT=4

nnodes = int(sys.argv[1])

mempn_gib = int(sys.argv[2])
mempn = mempn_gib*GiB

nrankspn = int(sys.argv[3])
nranks = nnodes*nrankspn

shdeg  = int(sys.argv[4])
resdeg = float(sys.argv[5])

nparams=float(shdeg*(shdeg+2))
ny=int(1.0/resdeg)*360
nx=int(1.0/resdeg)*180-1
ndatpts=float(nx*ny)
nsampts=float((shdeg+1)*(2*shdeg+1))
npts=ndatpts+nsampts
ncoeffs=255

print("spherical harmonic degree: "+str(shdeg))
print("resolution degree: "+str(resdeg))
print("coefficients: "+str(ncoeffs))
print("parameters: "+str(int(nparams)))
print("nx points: "+str(int(nx)))
print("ny points: "+str(int(ny)))
print("data points: "+str(int(ndatpts)))
print("sampling points: "+str(int(nsampts)))
print("data+sampling points: "+str(int(npts)))

nlocdatpts=int(ndatpts/nranks)+1
nlocsampts=int(nsampts/nranks)+1
nlocpts=int(npts/nranks)+1

print("local data points: "+str(int(nlocdatpts)))
print("local sampling points: "+str(int(nlocsampts)))
print("local data+sampling points: "+str(int(nlocpts)))

lab = []
mem = []


precalc = nparams*3 + (shdeg+1)

lab.append("mod_wmam_020") # mem[0]
mem.append(SIZEOF_DOUBLE*(1*nparams + 10*nlocpts + precalc) + SIZEOF_INT*(2*(nlocpts+2) + 2*nranks))


lab.append("cpt_dat_vals_p")
mem.append(mem[0] + SIZEOF_DOUBLE*(1*ncoeffs + 1*(shdeg+1)))

lab.append("build_damp_space")
mem.append(mem[0] + SIZEOF_DOUBLE*(10*nlocsampts) + SIZEOF_DOUBLE*(1*ncoeffs + 1*(shdeg+1)))


lab.append("opt_pr_p3") # mem[3]
mem.append(mem[0] + SIZEOF_DOUBLE*(3 + 2*nparams + 3*nparams + 1*nlocpts))


lab.append("cpt_dat_vals_p2")
mem.append(mem[3] + SIZEOF_DOUBLE*(1*nparams + 2*(shdeg+1) + 1*(shdeg)))

lab.append("cptstd_dp")
mem.append(mem[3] + SIZEOF_DOUBLE*(2*nranks + 6))


lab.append("ssqgh_dp")
mem.append(mem[3] + SIZEOF_DOUBLE*(1*nparams) + SIZEOF_DOUBLE*(2*nparams) + SIZEOF_DOUBLE*(9 + 3*nparams))


lab.append("gc_step_p") # mem[7]
mem.append(mem[3] + SIZEOF_DOUBLE*(1*nlocpts))
lab.append("cpt_dat_vals_p2")
mem.append(mem[7] + SIZEOF_DOUBLE*(1*nparams + 2*(shdeg+1) + 1*(shdeg)))

lab.append("gc_step_p") # mem[9]
mem.append(mem[3] + SIZEOF_DOUBLE*(1*nparams))
lab.append("cpt_dat_vals_p2")
mem.append(mem[9] + SIZEOF_DOUBLE*(1*nparams + 2*(shdeg+1) + 1*(shdeg)))
lab.append("cptstd_dp")
mem.append(mem[9] + SIZEOF_DOUBLE*(2*nranks + 6))


lab.append("lsearch_p") # mem[12]
mem.append(mem[3] + SIZEOF_DOUBLE*(1*nparams))
lab.append("cpt_dat_vals_p2")
mem.append(mem[12] + SIZEOF_DOUBLE*(1*nparams + 2*(shdeg+1) + 1*(shdeg)))
lab.append("cptstd_dp")
mem.append(mem[12] + SIZEOF_DOUBLE*(2*nranks + 6))


print("\nEstimated percentage of node memory ("+str(mempn_gib)+" GiB) use for "+str(nrankspn)+" ranks per node.")
for i, l in enumerate(lab):
  mem[i] = mem[i]*nrankspn
  print(l+": "+str(round(100.0*float(mem[i])/float(mempn),3))+"%")
