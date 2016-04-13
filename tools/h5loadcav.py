#!/usr/bin/env python
"""Load text based cavity data into an HDF5 file

./loadcav.py data.h5/cav/foo/field axisData_foo.txt

create dataset 'field' in group '/cav/foo' in file 'data.h5'.

File and group will be created if they do not exist.
Load will fail rather than overwrite an existing dataset unless
'--force' is given.
"""

import sys, argparse
import h5py

def getargs():
    from argparse import ArgumentParser
    P = ArgumentParser()
    P.add_argument('h5out', help='Output hdf5 file.h5/[group/]dataset (group is optional)"')
    P.add_argument('infile', help='Input cavity data text file')
    P.add_argument('-f','--force', action='store_true', default=False, help="Overwrite existing dataset")
    args = P.parse_args()
    hname = args.h5out.split('/')
    if len(hname)<2:
        P.error('h5out file must include file name and dataset name (eg. foo.h5/dset)')
    args.h5out = (hname[0], '/'.join(hname[1:-1]), hname[-1]) # (file, group, dataset)
    return args

def main(args):
    hfile, hgrp, hset = args.h5out

    from numpy import loadtxt

    D = loadtxt(args.infile)

    with h5py.File(hfile,'a') as F:
        G = F.require_group(hgrp) if hgrp else F

        if args.force and hset in G:
            del G[hset]
        S = G.create_dataset(hset, data=D)

if __name__=='__main__':
    main(getargs())
