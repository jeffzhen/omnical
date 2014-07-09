#! /usr/bin/env python
import numpy as n
import sys,glob

NCHAN = 203
NPRM = 349 # 2 * nants + nubls
NANT = 61 # because 19,37,50 are excluded
time_dir = '~/psa_live/forlstbinning_omnical_2/'

for omnical_file in sys.argv[1:]:
    #get times for the files.
    chunk = omnical_file.split('_')[2]
    time_files = glob.glob(time_dir + 'zen*%s*O'%chunk)
    pol = omnical_file.split('.')[-2].split('_')[-2][0]
    d = n.fromfile(omnical_file, dtype=n.float32)
    print omnical_file, d.size
    ntime = d.size / NCHAN / NPRM
    d.shape = (ntime, NCHAN, NPRM)
    
    st = 0 
    for tfile in time_files:
        sub_times = n.fromfile(tfile,sep='\n') #the times in a uvfile
        ntpuv = len(sub_times) #times per uv file
        d_npz = {}
        d_npz['times'] = sub_times
        d_npz['iters'] = d[st:st+ntpuv:,:,0]
        d_npz['chi2_log'] = d[st:st+ntpuv:,:,1]
        d_npz['chi2_lin'] = d[st:st+ntpuv:,:,2]
        g = 10**d[st:st+ntpuv:,:,3:3+NANT]
        g = g * n.exp(-2j*n.pi*d[st:st+ntpuv:,:,3+NANT:3+2*NANT]/360.) # minus to change conjugation convention between omnical and aipy
        for i in xrange(NANT):
            if i >= 19: ant = i+1 # XXX cludge
            if i >= 37: ant = i+2 # XXX cludge
            if i >= 50: ant = i+3 # XXX cludge
            else: ant = i
            print ant,pol
            d_npz['%d,%s' % (ant,pol)] = g[...,i]
        for j in xrange((NPRM-(3+2*NANT))/2):
            print 'sep%d' % j
            d_npz['sep%d' % j] = d[st:st+ntpuv:,:,3+2*NANT+2*j] + 1j*d[st:st+ntpuv:,:,3+2*NANT+2*j+1]
        st = st + ntpuv

        tfile = tfile.split('/')[-1]
        nametstamp = '_'.join( omnical_file.split('_')[0:2] + tfile.split('_')[0] + omnical_file.split('_')[3:] )
        n.savez(omnical_file+'.npz', **d_npz)
