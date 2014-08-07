#!/usr/bin/env python

import aipy as ap
import numpy as np
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
import optparse, sys
import scipy.signal as ss
import matplotlib.pyplot as plt
FILENAME = "omnical_PSA128.py"

##########################Sub-class#############################
class RedundantCalibrator_PAPER(omni.RedundantCalibrator):
	def __init__(self, aa):
		nTotalAnt = len(aa)
		omni.RedundantCalibrator.__init__(self, nTotalAnt)
		self.aa = aa

	def compute_redundantinfo(self, badAntenna = [], badUBL = [], antennaLocationTolerance = 1e-6):
		self.antennaLocationTolerance = antennaLocationTolerance
		self.badAntenna = badAntenna
		self.badUBL = badUBL
		self.antennaLocation = np.zeros((self.nTotalAnt,3))
		for i in range(len(self.aa.ant_layout)):
			for j in range(len(self.aa.ant_layout[0])):
				self.antennaLocation[self.aa.ant_layout[i][j]] = np.array([i, j, 0])
		self.preciseAntennaLocation = np.array([ant.pos for ant in self.aa])
		omni.RedundantCalibrator.compute_redundantinfo(self)





######################################################################
##############Config parameters###################################
######################################################################
o = optparse.OptionParser()

ap.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('-t', '--tag', action = 'store', default = 'PSA128', help = 'tag name of this calibration')
o.add_option('-d', '--datatag', action = 'store', default = 'PSA128', help = 'tag name of this data set')
o.add_option('-i', '--infopath', action = 'store', default = '/data2/home/hz2ug/omnical/doc/redundantinfo_PSA128_17ba.bin', help = 'redundantinfo file to read')
o.add_option('--add', action = 'store_true', help = 'whether to enable crosstalk removal')
o.add_option('--nadd', action = 'store', type = 'int', default = -1, help = 'time steps w to remove additive term with. for running average its 2w + 1 sliding window.')
o.add_option('--datapath', action = 'store', default = None, help = 'uv file or binary file folder')
o.add_option('-o', '--outputpath', action = 'store', default = None, help = 'output folder')
o.add_option('-k', '--skip', action = 'store_true', help = 'whether to skip data importing from uv')
o.add_option('-u', '--newuv', action = 'store_true', help = 'whether to create new uv files with calibration applied')
o.add_option('-f', '--overwrite', action = 'store_true', help = 'whether to overwrite if the new uv files already exists')
o.add_option('--plot', action = 'store_true', help = 'whether to make plots in the end')

opts,args = o.parse_args(sys.argv[1:])
skip = opts.skip
create_new_uvs = opts.newuv
overwrite_uvs = opts.overwrite
make_plots = opts.plot
ano = opts.tag##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
dataano = opts.datatag#ano for existing data and lst.dat
sourcepath = opts.datapath
oppath = opts.outputpath
uvfiles = args
for uvf in uvfiles:
	if not os.path.isdir(uvf):
		uvfiles.remove(uvf)
if len(uvfiles) == 0:
	print "ERROR: No valid uv files detected in input. Exiting!"
	quit()

wantpols = {}
for p in opts.pol.split(','): wantpols[p] = ap.miriad.str2pol[p]
#wantpols = {'xx':ap.miriad.str2pol['xx']}#, 'yy':-6}#todo:

print "Reading calfile %s"%opts.cal,
sys.stdout.flush()
aa = ap.cal.get_aa(opts.cal, np.array([.15]))
print "Done"
sys.stdout.flush()
infopathxx = opts.infopath
infopathyy = opts.infopath

#badAntenna = [37]
#badUBL = []



infopaths = {'xx': infopathxx, 'yy': infopathyy}


removedegen = True
if opts.add and opts.nadd > 0:
	removeadditive = True
	removeadditiveperiod = opts.nadd
else:
	removeadditive = False
	removeadditiveperiod = -1

needrawcal = False #if true, (generally true for raw data) you need to take care of having raw calibration parameters in float32 binary format freq x nant
#rawpaths = {'xx':"testrawphasecalparrad_xx", 'yy':"testrawphasecalparrad_yy"}



keep_binary_data = True
keep_binary_calpar = True


converge_percent = 0.001
max_iter = 20
step_size = .3

######################################################################
######################################################################
######################################################################

########Massage user parameters###################################
sourcepath += '/'
utcPath = sourcepath + 'miriadextract_' + dataano + "_localtime.dat"
lstPath = sourcepath + 'miriadextract_' + dataano + "_lsthour.dat"

####get some info from the first uvfile   ################
print "Getting some basic info from %s"%uvfiles[0],
sys.stdout.flush()
uv=ap.miriad.UV(uvfiles[0])
nfreq = uv.nchan;
nant = uv['nants']
#startfreq = uv['sfreq']
#dfreq = uv['sdf']
del(uv)
print "Done."
sys.stdout.flush()




###start reading miriads################
if skip:
	print FILENAME + " MSG: SKIPPED reading uvfiles. Reading binary data files directly...",
	sys.stdout.flush()
	with open(utcPath) as f:
		timing = f.readlines()
		timing = [t.replace('\n','') for t in timing]
	data = [np.fromfile(sourcepath + 'data_' + dataano + '_' + key, dtype = 'complex64').reshape((len(timing), nfreq, len(aa) * (len(aa) + 1) / 2)) for key in wantpols.keys()]
	print "Done."
	sys.stdout.flush()

else:
	print FILENAME + " MSG:",  len(uvfiles), "uv files to be processed for " + ano
	sys.stdout.flush()
	data, t, timing, lst = omni.importuvs(uvfiles, np.concatenate([[[i,j] for i in range(j + 1)] for j in range(len(aa))]), wantpols, timingTolerance=100)#, nTotalAntenna = len(aa))
	print FILENAME + " MSG:",  len(t), "slices read."
	sys.stdout.flush()

	if keep_binary_data:
		print FILENAME + " MSG: saving binary data to disk...",
		sys.stdout.flush()
		f = open(utcPath,'w')
		for time in timing:
			f.write("%s\n"%time)
		f.close()
		f = open(lstPath,'w')
		for l in lst:
			f.write("%s\n"%l)
		f.close()
		for p,key in zip(range(len(wantpols)), wantpols.keys()):
			data[p].tofile(sourcepath + 'data_' + dataano + '_' + key)
		print "Done."
		sys.stdout.flush()
print FILENAME + " MSG: data time range UTC: %s to %s"%(timing[0], timing[-1])
sys.stdout.flush()
####create redundant calibrators################
#calibrators = [omni.RedundantCalibrator(nant, info = infopaths[key]) for key in wantpols.keys()]
calibrators = {}
omnigains = {}
adds = {}
for p, key in zip(range(len(data)), wantpols.keys()):

	calibrators[key] = RedundantCalibrator_PAPER(aa)
	calibrators[key].read_redundantinfo(infopaths[key], verbose=False)
	calibrators[key].nTime = len(timing)
	calibrators[key].nFrequency = nfreq

	###prepare rawCalpar for each calibrator and consider, if needed, raw calibration################
	if needrawcal:
		raise Exception("Raw calpar not coded for the case where you need starting raw calibration!")
	else:
		calibrators[key].rawCalpar = np.zeros((calibrators[key].nTime, calibrators[key].nFrequency, 3 + 2 * (calibrators[key].Info.nAntenna + calibrators[key].Info.nUBL)),dtype='float32')


	####calibrate################
	calibrators[key].removeDegeneracy = removedegen
	calibrators[key].convergePercent = converge_percent
	calibrators[key].maxIteration = max_iter
	calibrators[key].stepSize = step_size

	################first round of calibration	#########################
	print FILENAME + " MSG: starting calibration on %s %s."%(dataano, key), calibrators[key].nTime, calibrators[key].nFrequency,
	sys.stdout.flush()
	additivein = np.zeros_like(data[p])
	calibrators[key].logcal(data[p], additivein, verbose=True)
	additivein[:,:,calibrators[key].Info.subsetbl] = calibrators[key].lincal(data[p], additivein, verbose=True)
	#######################remove additive###############################
	if removeadditive:
		nadditiveloop = 1
		for i in range(nadditiveloop):
			weight = ss.convolve(np.ones(additivein.shape[0]), np.ones(removeadditiveperiod * 2 + 1), mode='same')
			for f in range(additivein.shape[1]):#doing for loop to save memory usage at the expense of negligible time
				additivein[:,f] = ss.convolve(additivein[:,f], np.ones(removeadditiveperiod * 2 + 1)[:, None], mode='same')/weight[:, None]
			calibrators[key].computeUBLFit = False
			if i == nadditiveloop - 1:
				calibrators[key].lincal(data[p], additivein, verbose=True)
			else:
				additivein[:,:,calibrators[key].Info.subsetbl] = additivein[:,:,calibrators[key].Info.subsetbl] + calibrators[key].lincal(data[p], additivein, verbose=True)
	print "Done."
	sys.stdout.flush()
	#######################save results###############################
	calibrators[key].utctimes = timing
	omnigains[key] = calibrators[key].get_omnigain()
	adds[key] = additivein
	if keep_binary_calpar:
		print FILENAME + " MSG: saving calibration results on %s %s."%(dataano, key),
		sys.stdout.flush()
		#Zaki: catch these outputs and save them to wherever you like

		calibrators[key].get_calibrated_data(data[p])
		calibrators[key].get_omnichisq()

		calibrators[key].get_omnifit()
if create_new_uvs:
	print FILENAME + " MSG: saving new uv files",
	sys.stdout.flush()
	infos = {}
	for key in wantpols.keys():
		infos[key] = omni.read_redundantinfo(infopaths[key])
	omni.apply_omnigain_uvs(uvfiles, omnigains, calibrators[wantpols.keys()[0]].totalVisibilityId, infos, wantpols, oppath, ano, adds= adds, verbose = True, overwrite = overwrite_uvs)
	print "Done"
	sys.stdout.flush()
if make_plots:
	for p,pol in zip(range(len(wantpols)), wantpols.keys()):
		plt.subplot(1,len(wantpols),p+1)
		plt.imshow(calibrators[pol].rawCalpar[:,:,2], vmin = 0, vmax = np.nanmax(calibrators[wantpols.keys()[0]].rawCalpar[:,50:-50:5,2]), interpolation='nearest')
	plt.colorbar()
	plt.show()