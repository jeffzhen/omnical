import aipy as ap
import numpy as np
import commands, os, time, math, ephem
import calibration_omni as omni
FILENAME = "omnical.py"

##############Config parameters###################################
latP = -0.53619181096511903
lonP = 0.37399448506783717

for psa in range(918,920):
	ano = 'psa' + str(psa)##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical /data4/raw_data/psa903-930/psa903
	uvfiles = commands.getoutput('ls /data4/raw_data/psa903-930/' + ano + '/*.uvcRRE -d').split()
	wantpols = {'xx':-5, 'yy':-6}

	infopaths = {'xx':'./redundantinfo_PSA32.txt', 'yy':'./redundantinfo_PSA32.txt'}
	oppath = './results/'

	removedegen = 1

	needrawcal = True #if true, (generally true for raw data) you need to take care of having raw calibration parameters in float32 binary format freq x nant
	rawpaths = {'xx':"testrawphasecalparrad_xx", 'yy':"testrawphasecalparrad_yy"}
	############################################################
	sa = ephem.Observer()
	sa.lon = lonP
	sa.lat = latP
	sun = ephem.Sun()
	julDelta = 2415020 # =julian date - pyephem's Observer date
	print FILENAME + " MSG:",  len(uvfiles), "uv files to be processed for " + ano

	####read redundant info
	 
	info = [omni.read_redundantinfo(infopaths[key]) for key in wantpols.keys()]
	#print info[0]['bl1dmatrix']
	#exit(1)

	####get some info from the first uvfile
	uv=ap.miriad.UV(uvfiles[0])
	nfreq = uv.nchan;
	nant = uv['nants'] / 2 # 'nants' counting ant-pols, so divide 2
	startfreq = uv['sfreq']
	dfreq = uv['sdf']

	###read raw phase calibration prameters over frequencyfor each antenna, 203 by 32 in radiants; this part can be replaced
	if needrawcal:
		rawcalpar = np.asarray([np.fromfile(rawpaths[key], dtype="complex64").reshape(nfreq, nant) for key in wantpols.keys()])
		rawcorrection = np.zeros((len(wantpols), nfreq, nant*(nant+1)/2), dtype='complex64') + 1#to be dividing the data;  data/calpar = model
		for p in range(len(wantpols)):
			for i, bl in zip(range(len(info[p]['subsetbl'])), info[p]['subsetbl']):
				a1, a2 = info[p]['bl2d'][i]
				rawcorrection[p, :, bl] = np.conj(rawcalpar[p, :,a1]) *  rawcalpar[p, :,a2]

	####prepare processing
	deftime = 2000
	data = np.zeros((deftime, len(wantpols), nant * (nant + 1) / 2, nfreq), dtype = 'complex64')
	#sunpos = np.zeros((deftime, 2))
	t = []
	timing = []
	lst = []

	###start processing
	for uvfile in uvfiles:
		uv = ap.miriad.UV(uvfile)
		if len(timing) > 0:	
			print FILENAME + " MSG:",  timing[-1]#uv.nchan
		#print FILENAME + " MSG:",  uv['nants']
		currentpol = 0
		for preamble, rawd in uv.all():
			if len(t) < 1 or t[-1] != preamble[1]:#first bl of a timeslice
				t += [preamble[1]]
				sa.date = preamble[1] - julDelta
				#sun.compute(sa)
				timing += [sa.date.__str__()]
				lst += [(float(sa.sidereal_time()) * 24./2./math.pi)]
				if len(t) > len(data):
					print FILENAME + " MSG:",  "expanding number of time slices from", len(data), "to", len(data) + deftime
					data = np.concatenate((data, np.zeros((deftime, len(wantpols), nant * (nant + 1) / 2, nfreq), dtype = 'complex64'))) 
					#sunpos = np.concatenate((sunpos, np.zeros((deftime, 2))))
					#sunpos[len(t) - 1] = np.asarray([[sun.alt, sun.az]])
			for p, pol in zip(range(len(wantpols)), wantpols.keys()):
				if wantpols[pol] == uv['pol']:
					a1, a2 = preamble[2]
					bl = info[p]['bl1dmatrix'][a1, a2]
					if bl < info[p]['nbl']:
						#print info[p]['subsetbl'][info[p]['crossindex'][bl]],
						data[len(t) - 1, p, info[p]['subsetbl'][info[p]['crossindex'][bl]]] = rawd.data.astype('complex64')
		del(uv)

	print FILENAME + " MSG:",  len(t), "slices read.",
	#polstr = ['yx', 'xy', 'yy', 'xx', 'bla', 'bla', 'bla', 'bla']

	if not os.path.exists(oppath):
		os.makedirs(oppath)

	reorder = (0,2,1)
	for p, pol in zip(range(len(wantpols)), wantpols.keys()):
		print "Writing polarization: " + pol,  np.array(data[:len(t), p].shape)[np.array(reorder)]
		if needrawcal:
			(np.transpose(data[:len(t), p],reorder)/rawcorrection[p, np.newaxis,:,:]).tofile(oppath + 'miriadextract_' + pol + '_' + ano)
		else:
			np.transpose(data[:len(t), p],reorder).tofile(oppath + 'miriadextract_' + pol + '_' + ano)	

	#np.savetxt('miriadextract_' + ano + "_sunpos.dat", sunpos[:len(t)], fmt='%8.5f')	

	f = open(oppath + 'miriadextract_' + ano + "_localtime.dat",'w')
	for time in timing:
		f.write("%s\n"%time)
	f.close()
	f = open(oppath + 'miriadextract_' + ano + "_lsthour.dat",'w')
	for l in lst:
		f.write("%s\n"%l)
	f.close()

	del(data)

	for p, pol in zip(range(len(wantpols)), wantpols.keys()):
		command = "./omnical " + oppath + 'miriadextract_' + pol + '_' + ano + " " + infopaths[pol] + " " + str(len(t)) + " " + str(nfreq) + " "  + str(nant) + " " + str(removedegen)# + " " + oppath + 'miriadextract_' + pol + '_' + ano + ".omnical"
		print FILENAME + " MSG: System call: ",  command
		os.system(command)
		print np.fromfile(oppath + 'miriadextract_' + pol + '_' + ano + ".omnical", dtype = 'float32').reshape((len(t), nfreq, 3+2*(info[p]['nAntenna']+info[p]['nUBL'])))[:5,50,:3]






