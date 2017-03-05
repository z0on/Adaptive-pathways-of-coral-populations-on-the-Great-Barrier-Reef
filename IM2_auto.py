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

#infile="3kA_dadi.data"
#pop_ids=["skyblue","green4"]
#projections=[30,8]

dd = Misc.make_data_dict(infile)
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=False)
ns=data.sample_sizes
pts=[max(projections)+15,max(projections)+25,max(projections)+35]
np.set_printoptions(precision=3)     

#pylab.figure()
#dadi.Plotting.plot_single_2d_sfs(data, vmin=1)
#plt.show()

#-------------------
#  IM model with flexible initial sizes
def IM2(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu01,nu02,nu1,nu2,T,m12,m21)
    Isolation-with-migration model with exponential pop growth.
    Populations split into sizes nu01 and nu02.
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu01,nu02,nu1,nu2,T,m12,m21 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu01 * (nu1/nu01)**(t/T)  
    nu2_func = lambda t: nu02 * (nu2/nu02)**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

im2=Numerics.make_extrap_log_func(IM2)

params=array([0.078,  2.136,  2.174,  2.063,  1.174,  0.514,  0.379])
upper_bound = [100,100, 100, 100, 100, 200,200]
lower_bound = [1e-3,1e-3,1e-3,1e-3, 1e-3,0,0]
poptg = dadi.Inference.optimize_log(params, data, im2, pts,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(params),
                                   maxiter=20)

model=im2(poptg,ns,pts)
ll_model = dadi.Inference.ll_multinom(model, data)
theta = dadi.Inference.optimal_sfs_scaling(model, data)

print 'im2 ',sys.argv[1],sys.argv[2],sys.argv[3],' params : ', poptg
print 'im2 ',sys.argv[1],sys.argv[2],sys.argv[3],' logLike :', ll_model
print 'im2 ',sys.argv[1],sys.argv[2],sys.argv[3],' theta:', theta

dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
                                    pop_ids =pop_ids)
plt.savefig('IM2_'+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+"_"+sys.argv[4]+"_"+sys.argv[5]+'.pdf')

