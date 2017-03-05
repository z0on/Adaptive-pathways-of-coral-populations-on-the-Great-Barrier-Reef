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
#infile="5kA_dadi.data"
#pop_ids=["O1"]
#projections=[32]

import os

# replace this with your appropriate dir name
# os.chdir("/Users/c-monstr/Documents/allRAD_august2015/digitifera/dadi")

dd = Misc.make_data_dict(infile)
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=True)
ns=data.sample_sizes
pts=[65,80,95]
np.set_printoptions(precision=3)     

#-------------------
gr=Numerics.make_extrap_log_func(Demographics1D.growth)

params=array([ 3,1])
upper_bound = [100, 100]
lower_bound = [0.01, 0.01]
poptg = dadi.Inference.optimize_log(params, data, gr, pts,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(params),
                                   maxiter=15)

model=gr(poptg,ns,pts)

print 'gr ',sys.argv[1],sys.argv[2],' params : ', poptg
ll_model = dadi.Inference.ll_multinom(model, data)
print 'gr ',sys.argv[1],sys.argv[2],' logLike :', ll_model
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print 'gr ',sys.argv[1],sys.argv[2],' theta:', theta

#dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
#                                    pop_ids = pop_ids)
#plt.savefig('gr_'+sys.argv[1]+sys.argv[2]+sys.argv[3]+'.pdf')

