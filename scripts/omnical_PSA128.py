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
o.add_option('--healthbar', action = 'store', type = 'float', default = 10, help = 'health threshold (0-100) over which an antenna is marked bad.')
o.add_option('-o', '--outputpath', action = 'store', default = ".", help = 'output folder')
o.add_option('-k', '--skip', action = 'store_true', help = 'whether to skip data importing from uv')
o.add_option('-u', '--newuv', action = 'store_true', help = 'whether to create new uv files with calibration applied')
o.add_option('-f', '--overwrite', action = 'store_true', help = 'whether to overwrite if the new uv files already exists')
o.add_option('--plot', action = 'store_true', help = 'whether to make plots in the end')
o.add_option('--crude', action = 'store_true', help = 'whether to apply crude calibration')

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
healthbar = opts.healthbar
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

need_crude_cal = opts.crude #if true, (generally true for raw data) call raw_calibrate on first time slice of data set
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
oppath += '/'
utcPath = sourcepath + 'miriadextract_' + dataano + "_localtime.dat"
lstPath = sourcepath + 'miriadextract_' + dataano + "_lsthour.dat"

####get some info from the first uvfile   ################
print "Getting some basic info from %s"%uvfiles[0],
sys.stdout.flush()
uv=ap.miriad.UV(uvfiles[0])
nfreq = uv.nchan;
nant = uv['nants']
sa = ephem.Observer()
sa.lon = uv['longitu']
sa.lat = uv['latitud']
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
	data = np.array([np.fromfile(sourcepath + 'data_' + dataano + '_' + key, dtype = 'complex64').reshape((len(timing), nfreq, len(aa) * (len(aa) + 1) / 2)) for key in wantpols.keys()])
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
		for qaz in timing:
			f.write("%s\n"%qaz)
		f.close()
		f = open(lstPath,'w')
		for l in lst:
			f.write("%s\n"%l)
		f.close()
		for p,key in zip(range(len(wantpols)), wantpols.keys()):
			data[p].tofile(sourcepath + 'data_' + dataano + '_' + key)
		print "Done."
		sys.stdout.flush()
sun = ephem.Sun()
sunpos  = np.zeros((len(timing), 2))
cenA = ephem.FixedBody()
cenA._ra = 3.5146
cenA._dec = -.75077
cenApos = np.zeros((len(timing), 2))
for nt,tm in zip(range(len(timing)),timing):
	sa.date = tm

	sun.compute(sa)
	sunpos[nt] = sun.alt, sun.az
	cenA.compute(sa)
	cenApos[nt] = cenA.alt, cenA.az
print FILENAME + " MSG: data time range UTC: %s to %s, sun altaz from (%f,%f) to (%f,%f)"%(timing[0], timing[-1], sunpos[0,0], sunpos[0,1], sunpos[-1,0], sunpos[-1,1])#, "CentaurusA altaz from (%f,%f) to (%f,%f)"%(cenApos[0,0], cenApos[0,1], cenApos[-1,0], cenApos[-1,1])
sys.stdout.flush()
####create redundant calibrators################
#calibrators = [omni.RedundantCalibrator(nant, info = infopaths[key]) for key in wantpols.keys()]
calibrators = {}
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
		print "Done"
		sys.stdout.flush()
	calibrators[key].diagnose(data = data[p], additiveout = additiveout, healthbar = healthbar)
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
	omni.apply_omnigain_uvs(uvfiles, omnigains, calibrators[wantpols.keys()[0]].totalVisibilityId, infos, wantpols, oppath, ano, adds= adds, verbose = True, overwrite = overwrite_uvs)
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
