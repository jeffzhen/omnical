#!/usr/bin/env python

import aipy as ap
import numpy as np
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
import optparse, sys
import scipy.signal as ss
import matplotlib.pyplot as plt
FILENAME = "omnical_X5.py"

##########################Sub-class#############################
class RedundantCalibrator_MITEoR(omni.RedundantCalibrator):
    def __init__(self):
        nTotalAnt = 64
        omni.RedundantCalibrator.__init__(self, nTotalAnt)
        self.totalVisibilityId = np.concatenate([[[i,j] for j in range(i, self.nTotalAnt)] for i in range(self.nTotalAnt)])

    def compute_redundantinfo(self, badAntenna = [], badUBL = [], antennaLocationTolerance = 1e-6):
        self.antennaLocationTolerance = antennaLocationTolerance
        self.badAntenna = badAntenna + range(16) + range(56, 60)
        #print self.badAntenna
        self.badUBL = badUBL
        self.antennaLocation = np.zeros((self.nTotalAnt,3))
        for i in range(self.nTotalAnt):
            self.antennaLocation[i] = np.array([(i/16)*2 + (i%4)/2, ((i/4)%4)*2 + (i%4)%2, 0])
        omni.RedundantCalibrator.compute_redundantinfo(self)




if __name__ == '__main__':
    ######################################################################
    ##############Config parameters###################################
    ######################################################################
    #o = optparse.OptionParser()

    #ap.scripting.add_standard_options(o, cal=True, pol=True)
    #o.add_option('-t', '--tag', action = 'store', default = 'PSA128', help = 'tag name of this calibration')
    #o.add_option('-d', '--datatag', action = 'store', default = 'PSA128', help = 'tag name of this data set')
    #o.add_option('-i', '--infopath', action = 'store', default = '/data2/home/hz2ug/omnical/doc/redundantinfo_PSA128_17ba.bin', help = 'redundantinfo file to read')
    #o.add_option('--add', action = 'store_true', help = 'whether to enable crosstalk removal')
    #o.add_option('--nadd', action = 'store', type = 'int', default = -1, help = 'time steps w to remove additive term with. for running average its 2w + 1 sliding window.')
    #o.add_option('--datapath', action = 'store', default = None, help = 'uv file or binary file folder')
    #o.add_option('--healthbar', action = 'store', default = '2', help = 'health threshold (0-100) over which an antenna is marked bad.')
    #o.add_option('-o', '--outputpath', action = 'store', default = ".", help = 'output folder')
    #o.add_option('-k', '--skip', action = 'store_true', help = 'whether to skip data importing from uv')
    #o.add_option('-u', '--newuv', action = 'store_true', help = 'whether to create new uv files with calibration applied')
    #o.add_option('-f', '--overwrite', action = 'store_true', help = 'whether to overwrite if the new uv files already exists')
    #o.add_option('--plot', action = 'store_true', help = 'whether to make plots in the end')
    #o.add_option('--crude', action = 'store_true', help = 'whether to apply crude calibration')

    make_plots = True
    calparano = '24.08.2014'##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
    dirpath = '/home/omniscope/fftt/calibration/fortran_code/'
    quat = 'q3'
    anos = [quat + letter for letter in ['A', 'B', 'C']]
    todofiles = [dirpath + '/logcal_todo' + ano + '.txt' for ano in anos]
    dataname = 'rawcal_visibilities'

    odfs = []
    for todofile in todofiles:
        with open(todofile) as f:
            odfs += [odf.replace('\n','') for odf in f.readlines()]


    healthbar = 2
    ubl_healthbar = 100

    infopaths = {}
    for pol in ['x','y']:
        infopaths[pol + pol] = '/home/omniscope/omnical/doc/redundantinfo_X5_%s%s.bin'%(quat, pol)
    wantpols = {'xx':0, 'yy':3}


    removedegen = True
    removeadditive = True

    keep_binary_data = False
    keep_binary_calpar = True


    converge_percent = 0.001
    max_iter = 50
    step_size = .3

    ######################################################################
    ######################################################################
    ######################################################################

    print FILENAME + " MSG:",  len(odfs), "odf files to be processed for ", anos
    sys.stdout.flush()

    print FILENAME + " MSG:", "Initializing Redundant Calibrators...",
    sys.stdout.flush()
    calibrators = {}
    for p, key in zip(range(len(wantpols)), wantpols.keys()):
        calibrators[key] = RedundantCalibrator_MITEoR()
        calibrators[key].read_redundantinfo(infopaths[key], verbose=False)
    print "Done."
    sys.stdout.flush()

    for odf in odfs:
        useful_int_keys = ['nChannels', 'nIntegrations', 'nAntennas']
        useful_float_keys = ['startFreq', 'endFreq', 'integrationTime']
        header = {}
        with open(odf + '/header.txt') as f:
            for l in f.readlines():
                if l.split()[0] in useful_float_keys:
                    header[l.split()[0]] = float(l.split()[1])
                if l.split()[0] in useful_int_keys:
                    header[l.split()[0]] = int(l.split()[1])
        rawdata = np.fromfile(odf + '/' + dataname, dtype = 'complex64').reshape((4, header['nIntegrations'], header['nChannels'], header['nAntennas'] * (header['nAntennas'] + 1) / 2))
        data = {}
        for pol in wantpols.keys():
            data[pol] = rawdata[wantpols[pol]]

        #sun = ephem.Sun()
        #sunpos  = np.zeros((len(timing), 2))
        #for nt,tm in zip(range(len(timing)),timing):
            #sa.date = tm

            #sun.compute(sa)
            #sunpos[nt] = sun.alt, sun.az
        #print FILENAME + " MSG: data time range UTC: %s to %s, sun altaz from (%f,%f) to (%f,%f)"%(timing[0], timing[-1], sunpos[0,0], sunpos[0,1], sunpos[-1,0], sunpos[-1,1])#, "CentaurusA altaz from (%f,%f) to (%f,%f)"%(cenApos[0,0], cenApos[0,1], cenApos[-1,0], cenApos[-1,1])
        #sys.stdout.flush()
        ####create redundant calibrators################
        for key in wantpols.keys():
            info = calibrators[key].Info.get_info()
            calibrators[key].nTime = header['nIntegrations']
            calibrators[key].nFrequency = header['nChannels']
            calibrators[key].rawCalpar = np.zeros((calibrators[key].nTime, calibrators[key].nFrequency, 3 + 2 * (calibrators[key].Info.nAntenna + calibrators[key].Info.nUBL)),dtype='float32')


            ####calibrate################
            calibrators[key].removeDegeneracy = removedegen
            calibrators[key].convergePercent = converge_percent
            calibrators[key].maxIteration = max_iter
            calibrators[key].stepSize = step_size

            ################first round of calibration  #########################
            print FILENAME + " MSG: starting calibration on %s %s. nTime = %i, nFrequency = %i ..."%(odf, key, calibrators[key].nTime, calibrators[key].nFrequency),
            sys.stdout.flush()
            timer = time.time()
            additivein = np.zeros_like(data[key])
            calibrators[key].logcal(data[key], additivein, verbose=True)
            additiveout = calibrators[key].lincal(data[key], additivein, verbose=True)
            #######################remove additive###############################
            if removeadditive:
                nadditiveloop = 1
                for i in range(nadditiveloop):
                    for f in range(header['nChannels']):
                        additivein[:,f,calibrators[key].Info.subsetbl] = np.average(additivein[:,f,calibrators[key].Info.subsetbl] + additiveout[:,f], axis = 0, weights = calibrators[key].rawCalpar[:,f,2])

                    calibrators[key].computeUBLFit = False
                    additiveout = calibrators[key].lincal(data[key], additivein, verbose=True)

            print "Done. %fmin"%(float(time.time()-timer)/60.)
            sys.stdout.flush()
            #######################save results###############################
            if keep_binary_calpar:
                print FILENAME + " MSG: saving calibration results on %s %s."%(odf, key),
                sys.stdout.flush()
                calibrators[key].rawCalpar.tofile(odf + 'a/' + calparano + "_%s.omnical"%key)
                if removeadditive:
                    additivein[0].tofile(odf + 'a/' + calparano + "_%s.omniadd"%key)
                print "Done"
                sys.stdout.flush()
            calibrators[key].diagnose(data = data[key], additiveout = additiveout, healthbar = healthbar, ubl_healthbar = ubl_healthbar)

        #if make_plots:
            #for p,pol in zip(range(len(wantpols)), wantpols.keys()):
                #plt.subplot(1,len(wantpols),p+1)
                #plot_data = (calibrators[pol].rawCalpar[:,:,2]/(len(calibrators[pol].Info.subsetbl)-calibrators[pol].Info.nAntenna - calibrators[pol].Info.nUBL))**.5
                #plt.imshow(plot_data, vmin = 0, vmax = (np.nanmax(calibrators[wantpols.keys()[0]].rawCalpar[:,30:-30:5,2])/(len(calibrators[pol].Info.subsetbl)-calibrators[pol].Info.nAntenna - calibrators[pol].Info.nUBL))**.5, interpolation='nearest')
            #plt.title('RMS fitting error per baseline')
            #plt.colorbar()
            #plt.show()
