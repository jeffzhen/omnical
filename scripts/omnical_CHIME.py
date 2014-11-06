#!/usr/bin/env python

import numpy as np
import numpy.linalg as la
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
import optparse, sys
import scipy.signal as ss
import matplotlib.pyplot as plt
FILENAME = "omnical_CHIME.py"


##########################Sub-class#############################
class RedundantCalibrator_CHIME(omni.RedundantCalibrator):
    def __init__(self, nAnt):
        omni.RedundantCalibrator.__init__(self, nAnt)
        self.totalVisibilityId = np.concatenate([[[j,i] for j in range(i, self.nTotalAnt)] for i in range(self.nTotalAnt)])
        self.ps_k = None

    def compute_redundantinfo(self, npz, badAntenna = [], badUBLpair = [], antennaLocationTolerance = 1e-6):
        self.antennaLocationTolerance = antennaLocationTolerance
        self.badAntenna = badAntenna
        #print self.badAntenna
        self.badUBLpair = badUBLpair
        self.antennaLocation = np.ones((len(npz['feed_positions']),3))
        self.antennaLocation[:, :2] = npz['feed_positions']
        omni.RedundantCalibrator.compute_redundantinfo(self)
    
    def find_R(self, true_phase):#take the k vectors of cas (txfx3) a and figure out the R matrix
        if self.rawCalpar == None or len(self.rawCalpar.shape) != 3:
            raise TypeError("rawCalpar not properly initialized. The calibration is likely not performed yet.")
        
        
        max_ew = np.max([abs(blv[0]) for blv in self.Info.ubl])
        good_ubls = np.array([u for u in range(len(self.Info.ubl)) if abs(abs(self.Info.ubl[u][0]) - max_ew) < self.antennaLocationTolerance])

        calvis = self.rawCalpar[:,:,3+2*self.Info.nAntenna::2] + 1.j*self.rawCalpar[:,:,3+2*self.Info.nAntenna+1::2]
        true_amp = np.mean(np.abs(calvis)[:,:,good_ubls], axis = 2)
        true_vis = (true_amp[:,:,None] * np.exp(1.j * true_phase))
        return calvis - true_vis

def mod2p(arr):
    return np.mod(arr + np.pi, 2*np.pi) - np.pi

######################################################################
##############Config parameters###################################
######################################################################
o = optparse.OptionParser()
o.add_option('-a', '--nant', action = 'store', type = 'int', default = 8, help = 'n antenna. not necessary')
o.add_option('-f', '--nchan', action = 'store', type = 'int', default = 1, help = 'n frequency chan. not necessary')
o.add_option('-p', '--pol', action = 'store', default = 'x', help = 'polarization of the data, x for xx and y for yy.')
o.add_option('-t', '--tag', action = 'store', default = 'CHIME', help = 'tag name of this calibration')
o.add_option('-d', '--datatag', action = 'store', default = 'CHIME', help = 'tag name of this data set')
o.add_option('-i', '--infopath', action = 'store', default = '/data2/home/hz2ug/omnical/doc/redundantinfo_PSA128_17ba.bin', help = 'redundantinfo file to read')
o.add_option('--add', action = 'store_true', help = 'whether to enable crosstalk removal')
o.add_option('--nadd', action = 'store', type = 'int', default = -1, help = 'time steps w to remove additive term with. for running average its 2w + 1 sliding window.')
o.add_option('--datapath', action = 'store', default = None, help = 'uv file or binary file folder')
o.add_option('--healthbar', action = 'store', default = '2', help = 'health threshold (0-100) over which an antenna is marked bad.')
o.add_option('-o', '--outputpath', action = 'store', default = ".", help = 'output folder')
o.add_option('--plot', action = 'store_true', help = 'whether to make plots in the end')
o.add_option('--crude', action = 'store_true', help = 'whether to apply crude calibration')
o.add_option('--ba', action = 'store', default = '', help = 'bad antenna number indices seperated by commas')
o.add_option('--bu', action = 'store', default = '', help = 'bad unique baseline indicated by ant pairs (seperated by .) seperated by commas: 1.2,3.4,10.11')

opts,args = o.parse_args(sys.argv[1:])
make_plots = opts.plot
ano = opts.tag##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
dataano = opts.datatag#ano for existing data and lst.dat
sourcepath = opts.datapath
oppath = opts.outputpath
uvfiles = args
nfreq = opts.nchan
nant = opts.nant
try:
    badAntenna = [int(i) for i in opts.ba.split(',')]
except:
    badAntenna = []
try:
    badUBLpair = np.array([[int(j) for j in i.split('.')] for i in opts.bu.split(',')])
except:
    badUBLpair = np.array([])
print 'Bad Antennas:', badAntenna
print 'Bad unique baselines:', badUBLpair

#print opts.healthbar, opts.healthbar.split(), len(opts.healthbar.split())
if len(opts.healthbar.split(',')) == 1:
    healthbar = float(opts.healthbar)
    ubl_healthbar = 100
elif len(opts.healthbar.split(',')) == 2:
    healthbar = float(opts.healthbar.split(',')[0])
    ubl_healthbar = float(opts.healthbar.split(',')[1])
else:
    raise Exception("User input healthbar option (--healthbar %s) is not recognized."%opts.healthbar)
for uvf in uvfiles:
    if not os.path.isfile(uvf):
        uvfiles.remove(uvf)
        print "WARNING: file path %s does not exist!"%uvf
if len(uvfiles) == 0:
    raise Exception("ERROR: No valid uv files detected in input. Exiting!")

#trivial initialization of wantpols
wantpols = {}
for p in opts.pol.split(','):
    wantpols[p] = p

infopaths = {}
for key in wantpols.keys():
    infopaths[key]= opts.infopath


removedegen = True
know_sim_phase = True
fixed_ants = [0, 1, 4]
if opts.add and opts.nadd > 0:
    removeadditive = True
    removeadditiveperiod = opts.nadd
else:
    removeadditive = False
    removeadditiveperiod = -1

need_crude_cal = opts.crude #if true, (generally true for raw data) call raw_calibrate on first time slice of data set
#rawpaths = {'xx':"testrawphasecalparrad_xx", 'yy':"testrawphasecalparrad_yy"}



keep_binary_data = False
keep_binary_calpar = True


converge_percent = 0.001
max_iter = 200
step_size = .3

######################################################################
######################################################################
######################################################################

########Massage user parameters###################################
#sourcepath += '/'
oppath += '/'

###start reading data################
print FILENAME + " MSG:",  len(uvfiles), "data files to be processed for " + ano
sys.stdout.flush()
nt = 0
data = [None,None]
ra_list = []
for uv_file in uvfiles:
    data_pack = np.load(uv_file)
    nf = data_pack['corr_data'].shape[0]
    nt = nt + data_pack['corr_data'].shape[2]
    ra_list += np.array(data_pack['ra']).tolist()
    freqs = np.array(data_pack['f_MHz'])
    if freqs.shape == ():
        freqs = freqs.reshape(1)
    if data[0] == None:
        data[0] = data_pack['corr_data'].astype('complex64').transpose(2,0,1)
    else:
        data[0] = np.concatenate((data[0], data_pack['corr_data'].astype('complex64').transpose(2,0,1)))
print FILENAME + " MSG:",  nt, "slices read in %i channels"%nfreq
sys.stdout.flush()


####create redundant calibrators################
#calibrators = [omni.RedundantCalibrator(nant, info = infopaths[key]) for key in wantpols.keys()]
calibrators = {}
omnigains = {}
adds = {}
cdata = [None,None]
for p, key in zip(range(len(data)), wantpols.keys()):

    calibrators[key] = RedundantCalibrator_CHIME(len(data_pack['feed_positions']))
    calibrators[key].compute_redundantinfo(data_pack, badAntenna = badAntenna, badUBLpair = badUBLpair)
    info = calibrators[key].Info.get_info()
    calibrators[key].nTime = nt
    calibrators[key].nFrequency = nfreq

    ###prepare rawCalpar for each calibrator and consider, if needed, raw calibration################
    if need_crude_cal:
        initant, solution_path, additional_solution_path, degen, _ = omni.find_solution_path(info)
        crude_calpar[key] = np.array([omni.raw_calibrate(data[p, 0, f], info, initant, solution_path, additional_solution_path, degen) for f in range(calibrators[key].nFrequency)])
        data[p] = omni.apply_calpar(data[p], crude_calpar[key], calibrators[key].totalVisibilityId)

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

    #####absolute calibration for simulation only#########
    if know_sim_phase:
        degenmatrix = la.inv(info['antloc'][fixed_ants]).dot(calibrators[key].rawCalpar[:,:,3+info['nAntenna']:3+info['nAntenna']*2][:,:,fixed_ants].transpose(0,2,1)).transpose(1,0,2)
        calibrators[key].rawCalpar[:,:,3+info['nAntenna']:3+info['nAntenna']*2] = calibrators[key].rawCalpar[:,:,3+info['nAntenna']:3+info['nAntenna']*2] - info['antloc'].dot(degenmatrix).transpose(1,2,0)
        calvis = calibrators[key].rawCalpar[:,:,3+info['nAntenna']*2::2] + 1.j * calibrators[key].rawCalpar[:,:,3+info['nAntenna']*2+1::2]
        calvis = calvis * np.exp(1.j * (info['ubl'].dot(degenmatrix).transpose(1,2,0)))
        calibrators[key].rawCalpar[:,:,3+info['nAntenna']*2::2] = np.real(calvis)
        calibrators[key].rawCalpar[:,:,3+info['nAntenna']*2+1::2] = np.imag(calvis)
    
    #####calculate casa position and R#########
    chime = ephem.Observer()
    chime.lon = -119.618/180*np.pi
    chime.lat = 49.320/180*np.pi
    chime.pressure = 0
    casa = ephem.FixedBody()
    casa._ra = "23:23:26"#"23:23:26"
    casa._dec = "58:48:00"
    casa_altaz = np.zeros((len(ra_list), 2))
    for i, ra in enumerate(ra_list):
        chime.date = '2000/1/1'
        chime.date = chime.date - ephem.second*86164.0905 * (chime.sidereal_time() - ra / 180 * np.pi)/(2*np.pi)
        casa.compute(chime)
        casa_altaz[i] = [casa.alt, casa.az]
    t1 = np.argmax(casa_altaz[:,0]) - 10
    t2 = np.argmax(casa_altaz[:,0]) + 10
    casa_ks =np.array([np.cos(casa_altaz[:,0])*np.sin(casa_altaz[:,1]), np.cos(casa_altaz[:,0])*np.cos(casa_altaz[:,1]), np.sin(casa_altaz[:,0])]).transpose()
    casa_ks = np.outer(-freqs / 299.792458, casa_ks).reshape((nf,nt,3)).transpose(1,0,2)
    true_phase = 2 * np.pi * np.array(info['ubl']).dot(casa_ks.transpose(0,2,1)).transpose(1,2,0)
    #plt.scatter(casa_ks[:,0], casa_ks[:,1])
    #plt.scatter(casa_ks[1300:1400, 0, 0], casa_ks[1300:1400, 0, 1])
    #plotrange = freqs / 299.792458
    #plt.axis([-plotrange, plotrange, -plotrange, plotrange])
    #plt.show()
    
    
    for uu in range(len(info['ubl'])):
        plt.subplot(2, len(info['ubl'])/2, uu)
        plt.plot(np.angle(data_pack['true_corr'][0, info['subsetbl'][info['crossindex'][info['ublindex'][uu][0][2]]], t1:t2]))

        phase = mod2p(casa_ks.dot(info['ubl'][uu]) * 2 * np.pi)
        plt.plot(phase[t1:t2])
        plt.plot(mod2p(phase[t1:t2].flatten()-np.angle(data_pack['true_corr'][0, info['subsetbl'][info['crossindex'][info['ublindex'][uu][0][2]]], t1:t2])))
        plt.title(info['ubl'][uu])
        plt.ylim([-np.pi,np.pi])
    plt.show()



    ###temporarily use true_corr as true_phase
    #true_phase = np.angle([data_pack['true_corr'][:, info['subsetbl'][info['crossindex'][info['ublindex'][uu][0][2]]], :] for uu in range(info['nUBL'])]).transpose(2,1,0)
    R = calibrators[key].find_R(true_phase)
    
    #####store various useful quantities#########
    calpar = 10**(calibrators[key].rawCalpar[:,:,3:3+info['nAntenna']])*np.exp(1.j*calibrators[key].rawCalpar[:,:,3+info['nAntenna']:3+info['nAntenna']*2])
    cdata[p] = omni.apply_calpar(data[p],calpar,calibrators[key].totalVisibilityId)
    chisq = calibrators[key].rawCalpar[:,:,2]/(len(info['crossindex']) - info['nUBL'] - info['nAntenna'])/(53.**2/24200000)
    omni.omniview(np.array([data[p][np.argmax(casa_altaz[:,0]),0],cdata[p][np.argmax(casa_altaz[:,0]),0]]), info)
    #######################save results###############################
    if keep_binary_calpar:
        print FILENAME + " MSG: saving calibration results on %s %s."%(dataano, key),
        sys.stdout.flush()
        #Zaki: catch these outputs and save them to wherever you like
        calibrators[key].rawCalpar.tofile(oppath + '/' + dataano + '_' + ano + "_%s_%i_%i_%i.omnical"%(key, calibrators[key].rawCalpar.shape[0], calibrators[key].rawCalpar.shape[1], calibrators[key].rawCalpar.shape[2]))
        #if removeadditive:
            #adds[key].tofile(oppath + '/' + dataano + '_' + ano + "_%s.omniadd"%key + str(removeadditiveperiod))
        #calibrators[key].get_calibrated_data(data[p])
        #calibrators[key].get_omnichisq()
        #calibrators[key].get_omnifit()
        print "Done"
        sys.stdout.flush()
    #calibrators[key].diagnose(data = data[p], additiveout = additiveout, healthbar = healthbar, ubl_healthbar = ubl_healthbar)
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

if make_plots:
    for p,pol in zip(range(len(wantpols)), wantpols.keys()):
        plt.subplot(1,len(wantpols),p+1)
        plot_data = (calibrators[pol].rawCalpar[:,:,2]/(len(calibrators[pol].Info.subsetbl)-calibrators[pol].Info.nAntenna - calibrators[pol].Info.nUBL))**.5
        plt.imshow(plot_data, vmin = 0, vmax = (np.nanmax(calibrators[wantpols.keys()[0]].rawCalpar[:,30:-30:5,2])/(len(calibrators[pol].Info.subsetbl)-calibrators[pol].Info.nAntenna - calibrators[pol].Info.nUBL))**.5, interpolation='nearest')
    plt.title('RMS fitting error per baseline')
    plt.colorbar()
    plt.show()
