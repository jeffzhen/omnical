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

    def compute_redundantinfo(self, npz, badAntenna = [], badUBL = [], antennaLocationTolerance = 1e-6):
        self.antennaLocationTolerance = antennaLocationTolerance
        self.badAntenna = badAntenna
        #print self.badAntenna
        self.badUBL = badUBL
        self.antennaLocation = np.ones((len(npz['feed_positions']),3))
        self.antennaLocation[:, :2] = npz['feed_positions']
        omni.RedundantCalibrator.compute_redundantinfo(self)



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


opts,args = o.parse_args(sys.argv[1:])
make_plots = opts.plot
ano = opts.tag##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
dataano = opts.datatag#ano for existing data and lst.dat
sourcepath = opts.datapath
oppath = opts.outputpath
uvfiles = args
nfreq = opts.nchan
nant = opts.nant

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


removedegen = False
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
data = [None]
for uv_file in uvfiles:
	data_pack = np.load(uv_file)
	nf = data_pack['corr_data'].shape[0]
	nt = nt + data_pack['corr_data'].shape[2]
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
for p, key in zip(range(len(data)), wantpols.keys()):

	calibrators[key] = RedundantCalibrator_CHIME(len(data_pack['feed_positions']))
	calibrators[key].compute_redundantinfo(data_pack)
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

	################first round of calibration	#########################
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

	#####for simulation only#########
	if know_sim_phase:
		calibrators[key].rawCalpar[:,:,3+info['nAntenna']:3+info['nAntenna']*2] = calibrators[key].rawCalpar[:,:,3+info['nAntenna']:3+info['nAntenna']*2] - info['antloc'].dot(la.inv(info['antloc'][fixed_ants]).dot(calibrators[key].rawCalpar[:,:,3+info['nAntenna']:3+info['nAntenna']*2][:,:,fixed_ants].transpose(0,2,1)).transpose(1,0,2)).transpose(1,2,0)
	#######################save results###############################
	if keep_binary_calpar:
		print FILENAME + " MSG: saving calibration results on %s %s."%(dataano, key),
		sys.stdout.flush()
		#Zaki: catch these outputs and save them to wherever you like
		calibrators[key].rawCalpar.tofile(oppath + '/' + dataano + '_' + ano + "_%s.omnical"%key)
		if removeadditive:
			adds[key].tofile(oppath + '/' + dataano + '_' + ano + "_%s.omniadd"%key + str(removeadditiveperiod))
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

if make_plots:
	for p,pol in zip(range(len(wantpols)), wantpols.keys()):
		plt.subplot(1,len(wantpols),p+1)
		plot_data = (calibrators[pol].rawCalpar[:,:,2]/(len(calibrators[pol].Info.subsetbl)-calibrators[pol].Info.nAntenna - calibrators[pol].Info.nUBL))**.5
		plt.imshow(plot_data, vmin = 0, vmax = (np.nanmax(calibrators[wantpols.keys()[0]].rawCalpar[:,30:-30:5,2])/(len(calibrators[pol].Info.subsetbl)-calibrators[pol].Info.nAntenna - calibrators[pol].Info.nUBL))**.5, interpolation='nearest')
	plt.title('RMS fitting error per baseline')
	plt.colorbar()
	plt.show()
