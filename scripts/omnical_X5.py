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
        self.totalVisibilityId = np.concatenate([[[i,j] for j in range(i + 1)] for i in range(self.nTotalAnt)])

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
    calparano = '2014Aug23'##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
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

    infopaths = {'xx':'/home/omniscope/omnical/doc/redundantinfo_X5_5ba.bin', 'yy':'/home/omniscope/omnical/doc/redundantinfo_X5_5ba.bin'}
    wantpols = {'xx':0, 'yy':3}


    removedegen = True
    removeadditive = True

    keep_binary_data = False
    keep_binary_calpar = True


    converge_percent = 0.001
    max_iter = 20
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

    odf = '/home/omniscope/data/X5/day9/night/x5_2013-07-29_04_59_21_q4_167lo.odf'
    usefulkeys = ['startFreq', 'endFreq', 'nChannels', 'nIntegrations', 'integrationTime', 'nAntennas']
    header = {}
    with open(odf + '/header.txt') as f:
        for l in f.readlines():
            if l.split()[0] in usefulkeys:
                header[l.split()[0]] = float(l.split()[1])
    data = np.fromfile(odf + '/' + dataname, dtype = 'complex64').reshape((4, int(header['nIntegrations']), int(header['nChannels']), int(header['nAntennas'] * (header['nAntennas'] + 1) / 2)))
    print FILENAME + " MSG:",  header['nIntegrations'], "slices read."
    sys.stdout.flush()
    quit()

    #sun = ephem.Sun()
    #sunpos  = np.zeros((len(timing), 2))
    #for nt,tm in zip(range(len(timing)),timing):
        #sa.date = tm

        #sun.compute(sa)
        #sunpos[nt] = sun.alt, sun.az
    #print FILENAME + " MSG: data time range UTC: %s to %s, sun altaz from (%f,%f) to (%f,%f)"%(timing[0], timing[-1], sunpos[0,0], sunpos[0,1], sunpos[-1,0], sunpos[-1,1])#, "CentaurusA altaz from (%f,%f) to (%f,%f)"%(cenApos[0,0], cenApos[0,1], cenApos[-1,0], cenApos[-1,1])
    #sys.stdout.flush()
    ####create redundant calibrators################
    #calibrators = [omni.RedundantCalibrator(nant, info = infopaths[key]) for key in wantpols.keys()]
    omnigains = {}
    adds = {}
    for p, key in zip(range(len(data)), wantpols.keys()):

        calibrators[key] = RedundantCalibrator_PAPER(aa)
        calibrators[key].read_redundantinfo(infopaths[key], verbose=False)
        info = calibrators[key].Info.get_info()
        calibrators[key].nTime = len(timing)
        calibrators[key].nFrequency = nfreq

        ###prepare rawCalpar for each calibrator and consider, if needed, raw calibration################
        if need_crude_cal:
            initant, solution_path, additional_solution_path, degen, _ = omni.find_solution_path(info)

            crude_calpar = np.array([omni.raw_calibrate(data[p, 0, f], info, initant, solution_path, additional_solution_path, degen) for f in range(calibrators[key].nFrequency)])
            data[p] = omni.apply_calpar(data[p], crude_calpar, calibrators[key].totalVisibilityId)

        calibrators[key].rawCalpar = np.zeros((calibrators[key].nTime, calibrators[key].nFrequency, 3 + 2 * (calibrators[key].Info.nAntenna + calibrators[key].Info.nUBL)),dtype='float32')


        ####calibrate################
        calibrators[key].removeDegeneracy = removedegen
        calibrators[key].convergePercent = converge_percent
        calibrators[key].maxIteration = max_iter
        calibrators[key].stepSize = step_size

        ################first round of calibration  #########################
        print FILENAME + " MSG: starting calibration on %s %s. nTime = %i, nFrequency = %i ..."%(dataano, key, calibrators[key].nTime, calibrators[key].nFrequency),
        sys.stdout.flush()
        timer = time.time()
        additivein = np.zeros_like(data[p])
        calibrators[key].logcal(data[p], additivein, verbose=True)
        additiveout = calibrators[key].lincal(data[p], additivein, verbose=True)
        #######################remove additive###############################
        if removeadditive:
            nadditiveloop = 1
            for i in range(nadditiveloop):
                additivein[:,:,calibrators[key].Info.subsetbl] = additivein[:,:,calibrators[key].Info.subsetbl] + additiveout
                weight = ss.convolve(np.ones(additivein.shape[0]), np.ones(removeadditiveperiod * 2 + 1), mode='same')
                for f in range(additivein.shape[1]):#doing for loop to save memory usage at the expense of negligible time
                    additivein[:,f] = ss.convolve(additivein[:,f], np.ones(removeadditiveperiod * 2 + 1)[:, None], mode='same')/weight[:, None]
                calibrators[key].computeUBLFit = False
                additiveout = calibrators[key].lincal(data[p], additivein, verbose=True)

        print "Done. %fmin"%(float(time.time()-timer)/60.)
        sys.stdout.flush()
        #######################save results###############################
        calibrators[key].utctimes = timing
        omnigains[key] = calibrators[key].get_omnigain()
        adds[key] = additivein
        if keep_binary_calpar:
            print FILENAME + " MSG: saving calibration results on %s %s."%(dataano, key),
            sys.stdout.flush()
            #Zaki: catch these outputs and save them to wherever you like
            calibrators[key].rawCalpar.tofile(oppath + '/' + dataano + '_' + ano + "_xx.omnical")
            if removeadditive:
                adds[key].tofile(oppath + '/' + dataano + '_' + ano + "_xx.omniadd"+str(removeadditiveperiod))
            #calibrators[key].get_calibrated_data(data[p])
            #calibrators[key].get_omnichisq()
            #calibrators[key].get_omnifit()
            print "Done"
            sys.stdout.flush()
        calibrators[key].diagnose(data = data[p], additiveout = additiveout, healthbar = healthbar, ubl_healthbar = ubl_healthbar)
        ##print bad_ant_meter
        #nbad = 0
        #badstr = ''
        #for a in range(len(bad_ant_meter)):
            #if bad_ant_meter[a] > healthbar:
                #badstr += (str(a) + ',')
                #nbad += 1
        #if nbad > 0 :
            #print "BAD ANTENNA", badstr
        #print "%i NEW BAD ANTENNA(S)"%nbad
    if create_new_uvs:
        print FILENAME + " MSG: saving new uv files",
        sys.stdout.flush()
        infos = {}
        for key in wantpols.keys():
            infos[key] = omni.read_redundantinfo(infopaths[key])
        omni.apply_omnigain_uvs(uvfiles, omnigains, calibrators[wantpols.keys()[0]].totalVisibilityId, infos, wantpols, oppath, ano, adds= adds, verbose = True, comment = '_'.join(sys.argv), overwrite = overwrite_uvs)
        print "Done"
        sys.stdout.flush()
    if make_plots:
        for p,pol in zip(range(len(wantpols)), wantpols.keys()):
            plt.subplot(1,len(wantpols),p+1)
            plot_data = (calibrators[pol].rawCalpar[:,:,2]/(len(calibrators[pol].Info.subsetbl)-calibrators[pol].Info.nAntenna - calibrators[pol].Info.nUBL))**.5
            plt.imshow(plot_data, vmin = 0, vmax = (np.nanmax(calibrators[wantpols.keys()[0]].rawCalpar[:,30:-30:5,2])/(len(calibrators[pol].Info.subsetbl)-calibrators[pol].Info.nAntenna - calibrators[pol].Info.nUBL))**.5, interpolation='nearest')
        plt.title('RMS fitting error per baseline')
        plt.colorbar()
        plt.show()
