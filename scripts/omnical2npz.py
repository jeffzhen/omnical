#! /usr/bin/env python
import numpy as n
import sys

NCHAN = 203
NPRM = 353
NANT = 63 # because 37 was excluded

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
        if i >= 37: ant = i+1 # XXX cludge
        else: ant = i
        print ant,pol
        d_npz['%d,%s' % (ant,pol)] = g[...,i]
    for j in xrange((NPRM-(3+2*NANT))/2):
        print 'sep%d' % j
        d_npz['sep%d' % j] = d[:,:,3+2*NANT+2*j] + 1j*d[:,:,3+2*NANT+2*j+1]
    n.savez(omnical_file+'.npz', **d_npz)
