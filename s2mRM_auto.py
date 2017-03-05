#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF')
import dadi
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from dadi import Misc,Spectrum,Numerics,PhiManip,Integration,Demographics1D,Demographics2D
import sys
infile=sys.argv[1]
pop_ids=[sys.argv[2],sys.argv[3]]
projections=[int(sys.argv[4]),int(sys.argv[5])]

import os

# replace this with your appropriate dir name
# os.chdir("/Users/c-monstr/Documents/allRAD_august2015/digitifera/dadi")

dd = Misc.make_data_dict(infile)
#dd = Misc.make_data_dict("5kA_dadi.data")
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=True)
ns=data.sample_sizes
pts=[65,80,95]
np.set_printoptions(precision=3)     

#-------------------
# split with migration, with recent change in migration

def s2mRM(params, ns, pts):
    nu1,nu2,T,M12,M21,M12r,M21r= params
    Tr=0.01
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx,phi)
    phi = Integration.two_pops(phi,xx, T, nu1, nu2, m12=M12, m21=M21)
    phi = Integration.two_pops(phi,xx, Tr, nu1, nu2, m12=M12r, m21=M21r)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
 
S2MRM=Numerics.make_extrap_log_func(s2mRM)

params=array([ 1,1,1,30,30,30,30])
upper_bound = [10, 10, 10, 300,300,300,300]
lower_bound = [0.01, 0.01, 0.01,0.1,0.1,0.1,0.1]
poptg = dadi.Inference.optimize_log(params, data, S2MRM, pts,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(params),
                                   maxiter=15)

model=S2MRM(poptg,ns,pts)
print 's2mRM ',sys.argv[1],sys.argv[2],sys.argv[3],' params : ', poptg
ll_model = dadi.Inference.ll_multinom(model, data)
print 's2mRM ',sys.argv[1],sys.argv[2],sys.argv[3],' logLike :', ll_model
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print 's2mRM ',sys.argv[1],sys.argv[2],sys.argv[3],' theta:', theta

#dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
#                                    pop_ids = pop_ids)
#plt.savefig('s2mRM_'+sys.argv[1]+sys.argv[2]+sys.argv[3]+'.pdf')

