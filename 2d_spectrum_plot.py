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
dadi.Plotting.plot_single_2d_sfs(data, vmin=1)
#plt.show()
plt.savefig('2dAFS_'+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+"_"+sys.argv[4]+"_"+sys.argv[5]+'.pdf')
