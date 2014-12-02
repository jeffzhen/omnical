#!/usr/bin/env python

import aipy as ap
import numpy as np
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
import optparse, sys
FILENAME = "generate_redundantinfo.py"

if __name__ == '__main__':
    ######################################################################
    ##############Config parameters###################################
    ######################################################################

    o = optparse.OptionParser()

    o.add_option('-o', '--path', action = 'store', default = '', help = 'output name with path')
    o.add_option('-e', '--tol', action = 'store', type = 'float', default = 1e-2, help = 'tolerance of antenna location deviation when computing unique baselines.')
    o.add_option('-a', '--antenna', action = 'store', type = 'int', default = 64, help = 'Number of antennas.')
    o.add_option('--ba', action = 'store', default = '', help = 'bad antenna number indices seperated by commas')
    o.add_option('--bu', action = 'store', default = '', help = 'bad unique baseline indicated by ant pairs (seperated by .) seperated by commas: 1.2,3.4,10.11')
    o.add_option('--overwrite', action = 'store_true', help = 'overwrite if file exists')
    opts,args = o.parse_args(sys.argv[1:])

    if opts.path == '':
        raise Exception('Error: no output filename specified! Use -o to specify full name and path.')
    #if os.path.isfile(opts.path):
        #raise Exception('Error: output filename exists!')


    try:
        badAntenna = [int(i) for i in opts.ba.split(',')]
    except:
        badAntenna = []
    try:
        if opts.bu != '':
            badUBLpair = [[int(j) for j in i.split('.')] for i in opts.bu.split(',')]
        else:
            badUBLpair = []
    except:
        badUBLpair = []
    print 'Bad Antennas:', badAntenna
    print 'Bad unique baselines:', badUBLpair
    
    nant = opts.antenna
    
    calibrator = omni.RedundantCalibrator(nant)
    calibrator.badAntenna = badAntenna
    calibrator.badUBLpair = badUBLpair
    calibrator.antennaLocationTolerance = opts.tol
    timer = time.time()
    calibrator.compute_redundantinfo(verbose = True)
    print "Redundant info computed in %f minutes. %i good antenna, %i good UBL. Writing to %s..."%((time.time() - timer)/60., calibrator.Info.nAntenna, calibrator.Info.nUBL, opts.path),
    calibrator.write_redundantinfo(infoPath = opts.path, overwrite = opts.overwrite, verbose = True)
    print "Done."

