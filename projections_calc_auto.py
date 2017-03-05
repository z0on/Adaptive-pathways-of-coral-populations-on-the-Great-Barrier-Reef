#!/usr/bin/env python
import dadi
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from dadi import Misc,Spectrum,Numerics,PhiManip,Integration,Demographics1D,Demographics2D
import sys
infile=sys.argv[1]
popid=[sys.argv[2]]
proj=range(int(sys.argv[3]),int(sys.argv[4]))
dd = Misc.make_data_dict(infile)
for p in range(len(proj)):
    data = Spectrum.from_data_dict(dd, pop_ids=popid,projections=[proj[p]],polarized=False)
    print proj[p],data.S()
    
