#!/usr/bin/env python
"""Plot simulation results
"""

from __future__ import print_function

import sys

from matplotlib.pylab import *

from h5py import File

fname, _junk, gname = sys.argv[1].partition(':')

data = File(fname)
grp  = data[gname or '/']

simtype = grp.attrs['sim_type']

print('sim_type', simtype)

def show_vector(grp, vector='state'):
    'Show envelope size and angle as a function of s position'
    pos   = grp['pos']
    state = grp[vector]

    subplot(2,1,1)
    plot(pos, state[:,0], '-b',
         pos, state[:,2], '-r',
         pos, state[:,4], '-g')
    xlabel('s')
    ylabel('size')
    legend(['x','y','z'])

    subplot(2,1,2)
    plot(pos, state[:,1], '-b',
         pos, state[:,3], '-r',
         pos, state[:,5], '-g')
    xlabel('s')
    ylabel('angle')

def show_moment2(grp):
    pos = grp['pos'][:]
    avg = grp['moment0'][:]
    rms = grp['moment0_rms'][:]

    rmsp, rmsn = avg+rms, avg-rms

    for i,L in zip(range(6), ('x','px','y','py','z','pz')):
        subplot(3,2,i+1)
        plot(pos, rmsp[:,i], '-b',
             pos, avg [:,i], '-r',
             pos, rmsn[:,i], '-b')
        xlabel('s')
        ylabel(L)


def show_generic(grp):
    print("Unknown sim_type")

showsim = {
  'Vector': show_vector,
  'MomentMatrix2': show_moment2,
}

showsim.get(simtype)(grp)

show()
