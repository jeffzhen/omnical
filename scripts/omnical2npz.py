#! /usr/bin/env python
import numpy as n
import sys, optparse
import omnical.calibration_omni as omni


if __name__ == '__main__':
    ######################################################################
    ##############Config parameters###################################
    ######################################################################

    o = optparse.OptionParser()

    o.add_option('-i', '--infopath', action = 'store', default = '/data2/home/hz2ug/omnical/doc/redundantinfo_PSA128_17ba.bin', help = 'redundantinfo file to read')
    o.add_option('-a', '--nAntenna', action = 'store', type = 'int', default = 128, help = 'Number of antennas.')
    o.add_option('-f', '--nFrequency', action = 'store', type = 'int', default = 203, help = 'Number of frequency bins.')
    o.add_option('--overwrite', action = 'store_true', help = 'overwrite if file exists')
    opts,args = o.parse_args(sys.argv[1:])



    
    NCHAN = opts.nFrequency
    
    calibrator = omni.RedundantCalibrator(opts.nAntenna)
    print "Reading redundant info", opts.infopath, "...",
    sys.stdout.flush()
    calibrator.read_redundantinfo(opts.infopath)
    print "Done."
    sys.stdout.flush()

    NANT = calibrator.Info.nAntenna
    NPRM = 3 + 2 * (calibrator.Info.nAntenna + calibrator.Info.nUBL)
    
    for omnical_file in sys.argv[1:]:
        print "Processing", omnical_file, "..."
        sys.stdout.flush()
        pol = omnical_file.split('.')[-2].split('_')[-1][0]
        d = n.fromfile(omnical_file, dtype=n.float32)
        ntime = d.size / NCHAN / NPRM
        d.shape = (ntime, NCHAN, NPRM)
        d_npz = {}
        d_npz['iters'] = d[:,:,0]
        d_npz['chi2_log'] = d[:,:,1]
        d_npz['chi2_lin'] = d[:,:,2]
        g = 10**d[:,:,3:3+NANT] * n.exp(- 1.j * d[:,:,3+NANT:3+2*NANT]) # minus to change conjugation convention between omnical and aipy
        for i in xrange(NANT):
            d_npz['%d,%s' % (calibrator.Info.subsetant[i],pol)] = g[...,i]
        n.savez(omnical_file+'.npz', **d_npz)
