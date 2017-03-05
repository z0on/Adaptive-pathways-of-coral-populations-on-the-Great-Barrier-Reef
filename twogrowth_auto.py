#!/usr/bin/env python

import matplotlib
#matplotlib.use('PDF')
import dadi
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from dadi import Misc,Spectrum,Numerics,PhiManip,Integration,Demographics1D,Demographics2D
import sys
infile=sys.argv[1]
pop_ids=sys.argv[2]
projections=[int(sys.argv[3])]

import os

# replace this with your appropriate dir name
# os.chdir("/Users/c-monstr/Documents/allRAD_august2015/digitifera/dadi")

dd = Misc.make_data_dict(infile)
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=True)
ns=data.sample_sizes
pts=[65,80,95]
np.set_printoptions(precision=3)     

#-------------------
def two_growths((nu1,nu2, T1,T2), (n1,), pts):
      xx = Numerics.default_grid(pts) 
      phi = PhiManip.phi_1D(xx)        
      nu_func1 = lambda t: np.exp(np.log(nu1) * t/T1) 
      phi = Integration.one_pop(phi, xx, T1, nu_func1) 
      nu_func2 = lambda t: nu1*np.exp(np.log(nu2/nu1) * t/T2) 
      phi = Integration.one_pop(phi, xx, T2, nu_func2) 
      sfs = Spectrum.from_phi(phi, (n1,), (xx,)) 
      return sfs 

twog=Numerics.make_extrap_log_func(two_growths)

params=array([ 1,1,1,0.1])
upper_bound = [100,100,100,10]
lower_bound = [0.01,0.01,0.01,0.001]
poptg = dadi.Inference.optimize_log(params, data, twog, pts,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(params),
                                   maxiter=20)

model=twog(poptg,ns,pts)

print 'twog ',sys.argv[1],sys.argv[2],' params : ', poptg
ll_model = dadi.Inference.ll_multinom(model, data)
print 'twog ',sys.argv[1],sys.argv[2],' logLike :', ll_model
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print 'twog ',sys.argv[1],sys.argv[2],' theta:', theta

#dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
#                                    pop_ids = pop_ids)
#plt.savefig('twog_'+sys.argv[1]+sys.argv[2]+sys.argv[3]+'.pdf')

