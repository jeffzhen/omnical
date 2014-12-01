#!/usr/bin/env python

import numpy as np
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
import optparse, sys
FILENAME = "generate_redundantinfo_X5.py"

##########################Sub-class#############################
class RedundantCalibrator_MITEoR(omni.RedundantCalibrator):
    def __init__(self, nTotalAnt = 64):
        omni.RedundantCalibrator.__init__(self, int(nTotalAnt))
        self.totalVisibilityId = np.concatenate([[[i,j] for j in range(i, self.nTotalAnt)] for i in range(self.nTotalAnt)])

    def compute_redundantinfo(self, badAntenna = [], badUBLpair = [], antennaLocationTolerance = 1e-6, verbose = False):
        self.antennaLocationTolerance = antennaLocationTolerance
        self.badAntenna = badAntenna + range(16) + range(56, 60)
        #print self.badAntenna
        self.badUBLpair = badUBLpair
        self.antennaLocation = np.zeros((self.nTotalAnt,3))
        for i in range(self.nTotalAnt):
            self.antennaLocation[i] = np.array([(i/16)*2 + (i%4)/2, ((i/4)%4)*2 + (i%4)%2, 0])# + np.random.randn(3)/10000.
        omni.RedundantCalibrator.compute_redundantinfo(self, verbose = verbose)




if __name__ == '__main__':
    ######################################################################
    ##############Config parameters###################################
    ######################################################################

    o = optparse.OptionParser()

    o.add_option('-o', '--path', action = 'store', default = '', help = 'output name with path')
    o.add_option('-e', '--tol', action = 'store', type = 'float', default = 1e-2, help = 'tolerance of antenna location deviation when computing unique baselines.')
    o.add_option('-a', '--antenna', action = 'store', type = 'int', default = 64, help = 'Number of antennas.')
    o.add_option('--ba', action = 'store', default = '', help = 'bad antenna number indices seperated by commas')
    o.add_option('--bu', action = 'store', default = '', help = 'bad unique baseline indices seperated by commas')
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
        badUBLpair = np.array([])
    nant = opts.antenna
    
    print '# Antennas:', nant
    print 'Bad Antennas:', badAntenna
    print 'Bad unique baselines:', badUBLpair

    calibrator = RedundantCalibrator_MITEoR(nant)
    timer = time.time()
    calibrator.compute_redundantinfo(badAntenna = badAntenna, badUBLpair = badUBLpair, antennaLocationTolerance = opts.tol, verbose = True)
    print "Redundant info computed in %f minutes. %i good antenna, %i good UBL."%((time.time() - timer)/60., calibrator.Info.nAntenna, calibrator.Info.nUBL)
    calibrator.write_redundantinfo(infoPath = opts.path, overwrite = opts.overwrite, verbose = False)


