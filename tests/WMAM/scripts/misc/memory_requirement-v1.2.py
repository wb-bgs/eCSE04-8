#!/usr/bin/env python 

# python memory_requirement-v1.2.py 1 256 128 200  1.00
# 39.123%
# python memory_requirement-v1.2.py 1 256 128 300  0.50
# 39.524%
# python memory_requirement-v1.2.py 1 256 128 720  0.25
# 42.956%
# python memory_requirement-v1.2.py 1 256 128 1440 0.10
# 55.407%
# python memory_requirement-v1.2.py 1 256 128 2000 0.05
# 252.025%


import sys
import os
import math
import numpy
from datetime import datetime
from decimal import Decimal

KB = 1000
MB = KB*KB
GB = KB*MB

KiB = 1024
MiB = KiB*KiB
GiB = KiB*MiB

SIZEOF_DOUBLE=8
SIZEOF_INT=4

nnodes = int(sys.argv[1])

mempn_gb = int(sys.argv[2])
mempn = mempn_gb*GB

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


lab = []
mem = []


lab.append("mod_wmam_020") # mem[0]
if shdeg == 2000:
  mem.append(SIZEOF_DOUBLE*(390*MB + 3) + SIZEOF_INT*(320*MB + 5))
else:
  mem.append(SIZEOF_DOUBLE*(76*MB + 3) + SIZEOF_INT*(42*MB + 5))


lab.append("cpt_dat_vals_p")
mem.append(mem[0] + SIZEOF_DOUBLE*(1*ncoeffs + 1*(shdeg+1)) + SIZEOF_INT*(2*nranks))

lab.append("build_damp_space")
mem.append(mem[0] + SIZEOF_DOUBLE*(10*nsampts) + SIZEOF_DOUBLE*(1*ncoeffs + 1*(shdeg+1)) + SIZEOF_INT*(2*nranks))


lab.append("opt_pr_p3") # mem[3]
mem.append(mem[0] + SIZEOF_DOUBLE*(3 + 5*nparams + 1*npts))


lab.append("cpt_dat_vals_p")
mem.append(mem[3] + SIZEOF_DOUBLE*(4*nparams + 2*(shdeg+1)) + SIZEOF_INT*(2*nranks))

lab.append("cptstd_dp")
mem.append(mem[3] + SIZEOF_DOUBLE*(2*nranks + 6) + SIZEOF_INT*(2*nranks))


lab.append("ssqgh_dp")
mem.append(mem[3] + SIZEOF_DOUBLE*(3*nparams) + SIZEOF_INT*(2*nranks) + SIZEOF_DOUBLE*(9 + 3*nparams))


lab.append("gc_step_p") # mem[7]
mem.append(mem[3] + SIZEOF_DOUBLE*(1*npts))
lab.append("cpt_dat_vals_p")
mem.append(mem[7] + SIZEOF_DOUBLE*(4*nparams + 2*(shdeg+1)) + SIZEOF_INT*(2*nranks))

lab.append("gc_step_p") # mem[9]
mem.append(mem[3] + SIZEOF_DOUBLE*(1*nparams))
lab.append("cpt_dat_vals_p")
mem.append(mem[9] + SIZEOF_DOUBLE*(4*nparams + 2*(shdeg+1)) + SIZEOF_INT*(2*nranks))
lab.append("cptstd_dp")
mem.append(mem[9] + SIZEOF_DOUBLE*(2*nranks + 6) + SIZEOF_INT*(2*nranks))


lab.append("lsearch_p") # mem[12]
mem.append(mem[3] + SIZEOF_DOUBLE*(1*nparams))
lab.append("cpt_dat_vals_p")
mem.append(mem[12] + SIZEOF_DOUBLE*(4*nparams + 2*(shdeg+1)) + SIZEOF_INT*(2*nranks))
lab.append("cptstd_dp")
mem.append(mem[12] + SIZEOF_DOUBLE*(2*nranks + 6) + SIZEOF_INT*(2*nranks))


print("\nEstimated percentage of node memory ("+str(mempn_gb)+" GB) use for "+str(nrankspn)+" ranks per node.")
for i, l in enumerate(lab):
  mem[i] = mem[i]*nrankspn
  print(l+": "+str(round(100.0*float(mem[i])/float(mempn),3))+"%")
