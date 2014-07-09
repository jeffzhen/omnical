#! /usr/bin/env python
import numpy as n, aipy as a
import sys,glob

NCHAN = 203
NPRM = 349 # 2 * nants + nubls
NANT = 61 # because 19,37,50 are excluded
calfile = 'psa6240_v003'
aa = a.cal.get_aa(calfile, n.array([.150]))

for omnical_file in sys.argv[1:]:
    jd = omnical_file.split('_')[2]
    timefile = glob.glob('miriadextract*%s*localtime.dat' % jd)[0]
    pol = omnical_file.split('.')[-2].split('_')[-2][0]
    times, lsts = [], []
    for t in open(timefile).readlines():
        aa.date = t
        times.append(aa.get_jultime())
        lsts.append(aa.sidereal_time())
    times, lsts = n.array(times), n.array(lsts)
    d = n.fromfile(omnical_file, dtype=n.float32)
    print omnical_file, d.size
    ntime = d.size / NCHAN / NPRM
    d.shape = (ntime, NCHAN, NPRM)
    
    d_npz = {}
    d_npz['times'] = times
    d_npz['lsts'] = lsts
    d_npz['iters'] = d[:,:,0]
    d_npz['chi2_log'] = d[:,:,1]
    d_npz['chi2_lin'] = d[:,:,2]
    g = 10**d[:,:,3:3+NANT]
    g = g * n.exp(-2j*n.pi*d[:,:,3+NANT:3+2*NANT]/360.) # minus to change conjugation convention between omnical and aipy
    for i in xrange(NANT):
        if i >= 19: ant = i+1 # XXX cludge
        if i >= 37: ant = i+2 # XXX cludge
        if i >= 50: ant = i+3 # XXX cludge
        else: ant = i
        d_npz['%d,%s' % (ant,pol)] = g[...,i]
    for j in xrange((NPRM-(3+2*NANT))/2):
        d_npz['sep%d' % j] = d[:,:,3+2*NANT+2*j] + 1j*d[:,:,3+2*NANT+2*j+1]

    n.savez(omnical_file+'.npz', **d_npz)
