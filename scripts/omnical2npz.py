#! /usr/bin/env python
import numpy as n
import sys, optparse, os
import omnical.calibration_omni as omni


if __name__ == '__main__':
    ######################################################################
    ##############Config parameters###################################
    ######################################################################

    o = optparse.OptionParser()

    o.add_option('-i', '--infopath', action = 'store', default = '/data2/home/cacheng/calibrated128/psa6620/psa6620config3.bin', help = 'redundantinfo file to read')
    o.add_option('-a', '--nAntenna', action = 'store', type = 'int', default = 128, help = 'Number of antennas.')
    o.add_option('-f', '--nFrequency', action = 'store', type = 'int', default = 203, help = 'Number of frequency bins.')
    o.add_option('-o', '--outputpath', action = 'store', default = None, help = 'output folder. Same folder by default.')
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

        if opts.outputpath is None:
            if not (os.access(os.path.dirname(opts.outputpath), os.W_OK | os.X_OK)):
                raise IOError("Cannot write to %s."%opts.outputpath)
            else:
                oppath = omnical_file + '.npz'
        elif not (os.path.isdir(opts.outputpath) and os.access(opts.outputpath, os.W_OK | os.X_OK)):
            raise IOError("Cannot write to %s."%opts.outputpath)
        else:
            oppath = opts.outputpath + '/' + omnical_file.split('/')[-1] + '.npz'
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
        n.savez(oppath, **d_npz)
