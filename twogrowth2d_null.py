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

"""
dd = Misc.make_data_dict("5kA_dadi.data")
data = Spectrum.from_data_dict(dd, pop_ids=["O", "K"],projections=[40,30],polarized=True)
ns=data.sample_sizes
pts=[60,70,80]
"""

#-------------------
# model for two populations with no divergence and two growth periods

def two_growths2D((nu1,nu2, T1,T2), ns, pts):
      xx = Numerics.default_grid(pts) 
      phi = PhiManip.phi_1D(xx)        
      nu_func1 = lambda t: np.exp(np.log(nu1) * t/T1) 
      phi = Integration.one_pop(phi, xx, T1, nu_func1) 
      nu_func2 = lambda t: nu1*np.exp(np.log(nu2/nu1) * t/T2) 
      phi = Integration.one_pop(phi, xx, T2, nu_func2) 
      phi = PhiManip.phi_1D_to_2D(xx, phi)
      sfs = Spectrum.from_phi(phi, ns, (xx,xx)) 
      return sfs 

twog2d=Numerics.make_extrap_log_func(two_growths2D)

params=array([ 3,1,3,0.2])
upper_bound = [10,10,10,10]
lower_bound = [0.01,0.01,0.01,0.001]
poptg = dadi.Inference.optimize_log(params, data, twog2d, pts,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(params),
                                   maxiter=15)

model=twog2d(poptg,ns,pts)

print 'tg2d0 ',sys.argv[1],sys.argv[2],sys.argv[3],' params : ', poptg
ll_model = dadi.Inference.ll_multinom(model, data)
print 'tg2d0 ',sys.argv[1],sys.argv[2],sys.argv[3],' logLike :', ll_model
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print 'tg2d0 ',sys.argv[1],sys.argv[2],sys.argv[3],' theta:', theta

dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
                                    pop_ids = pop_ids)
plt.savefig('tg2d0_'+sys.argv[1]+sys.argv[2]+sys.argv[3]+'.pdf')


