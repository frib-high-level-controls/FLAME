#!/usr/bin/env python
"""Pack TLM output files into an HDF5 file
with layout a sub-set of 'flame -F hdf5'

The input directory must contain the MCSModel*.txt files.

./tlm2h5.py <inputdir> <outfile.h5> [grpname]
"""

import sys
from os.path import join

import numpy as np
from scipy import loadtxt
from h5py import File

basedir=sys.argv[1]
H5=File(sys.argv[2])
HG=H5.create_group(sys.argv[3]) if len(sys.argv)>3 else H5

H5.attrs['sim_type']='MomentMatrix'

V=loadtxt(join(basedir, 'MCSModelCenVec.txt'))
S=V[:,0]

HG.create_dataset('pos', data=S)

D=HG.create_dataset('moment0_env', dtype=V.dtype, shape=(V.shape[0], 7))
D[:,:-1] = V[:,1:]

V=loadtxt(join(basedir, 'MCSModelRmsVec.txt'))
if not np.all(V[:,0]==S):
  print 'S positions inconsistent'

D=HG.create_dataset('moment0_rms', dtype=V.dtype, shape=(V.shape[0], 7))
D[:,:-1] = V[:,1:]

V=loadtxt(join(basedir, 'MCSModelRf.txt'))
if not np.all(V[:,0]==S):
  print 'S positions inconsistent'

HG.create_dataset('IonEk_env', data=V[:,1])
HG.create_dataset('phis_env', data=V[:,2])
