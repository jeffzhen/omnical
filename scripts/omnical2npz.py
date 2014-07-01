#! /usr/bin/env python
import numpy as n
import sys

NCHAN = 203
NPRM = 353
NANT = 64

for omnical_file in sys.argv[1:]:
    pol = omnical_file.split('.')[-2].split('_')[-1][0]
    d = n.fromfile(omnical_file, dtype=n.float32)
    ntime = d.size / NCHAN / NPRM
    d.shape = (ntime, NCHAN, NPRM)
    d_npz = {}
    d_npz['iters'] = d[:,:,0]
    d_npz['chi2_log'] = d[:,:,1]
    d_npz['chi2_lin'] = d[:,:,2]
    g = 10**d[:,:,3:3+NANT]
    g = g * n.exp(-2j*n.pi*d[:,:,3+NANT:3+2*NANT]/360.) # minus to change conjugation convention between omnical and aipy
    for i in xrange(NANT):
        d_npz['%d,%s' % (i,pol)] = g[...,i]
    n.savez(omnical_file+'.npz', **d_npz)
