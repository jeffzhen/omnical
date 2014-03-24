import aipy as ap
import numpy as np
import ephem
import commands, os
import time
import calibration_omni as omni
FILENAME = "omnical.py"

##############Config parameters###################################
latP = -0.53619181096511903
lonP = 0.37399448506783717

ano = 'test'##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
uvfiles = ['test.uv']
infopath = ['./redundantinfo_PSA32.txt', './redundantinfo_PSA32.txt']

needrawcal = True #if true, (generally true for raw data) you need to take care of having raw calibration parameters in float32 binary format freq x nant
rawpathxx = "testrawphasecalparrad_xx"
rawpathyy = "testrawphasecalparrad_yy"

############################################################
sa = ephem.Observer()
sa.lon = lonP
sa.lat = latP
sun = ephem.Sun()
julDelta = 2415020 # =julian date - pyephem's Observer date
print FILENAME + " MSG:",  len(uvfiles), "uv files to be processed"

####read redundant info
 
info = [omni.read_redundantinfo(infopath[-5 - pol]) for pol in [-5, -6]]
#print info['bl2d']
#exit(1)

####get some info from the first uvfile
uv=ap.miriad.UV(uvfiles[0])
nfreq = uv.nchan;
nant = uv['nants'] / 2 # 'nants' counting ant-pols, so divide 2
startfreq = uv['sfreq']
dfreq = uv['sdf']

###read raw phase calibration prameters over frequencyfor each antenna, 203 by 32 in radiants; this part can be replaced
if needrawcal:
	rawcalparxx = np.fromfile(rawpathxx, dtype="complex64").reshape(nfreq, nant)
	rawcalparyy = np.fromfile(rawpathyy, dtype="complex64").reshape(nfreq, nant)
	rawcalpar = np.asarray([rawcalparxx, rawcalparyy])
	rawcorrection = np.zeros((2, nfreq, nant*(nant+1)/2), dtype='complex64') + 1#to be dividing the data;  data/calpar = model
	for pol in [-5,-6]:
		for i, bl in zip(range(len(info[-5 - pol]['subsetbl'])), info[-5 - pol]['subsetbl']):
			a1, a2 = info[-5 - pol]['bl2d'][i]
			rawcorrection[-5 - pol, :, bl] = np.conj(rawcalpar[-5 - pol, :,a1]) *  rawcalpar[-5 - pol, :,a2]

####prepare processing
deftime = 2000
data = np.zeros((deftime, 2, nant * (nant + 1) / 2, nfreq), dtype = 'complex64')
#sunpos = np.zeros((deftime, 2))
t = []
timing = []

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
			if len(t) > len(data):
				print FILENAME + " MSG:",  "expanding number of time slices from", len(data), "to", len(data) + deftime
				data = np.concatenate((data, np.zeros((deftime, 2, nant * (nant + 1) / 2, nfreq), dtype = 'complex64'))) 
				#sunpos = np.concatenate((sunpos, np.zeros((deftime, 2))))
				#sunpos[len(t) - 1] = np.asarray([[sun.alt, sun.az]])
		if -5 == uv['pol'] or -6 == uv['pol']:
			if currentpol != uv['pol']:
				bl = 0
				currentpol = uv['pol']
			data[len(t) - 1, -5 - currentpol, bl] = rawd.data.astype('complex64')
			bl += 1

print FILENAME + " MSG:",  len(t), "slices read."
polstr = ['yx', 'xy', 'yy', 'xx', 'bla', 'bla', 'bla', 'bla']


for pol in [-5, -6]:
	print FILENAME + " MSG:",  data[:len(t), -5 - pol].shape
	if needrawcal:
		(np.transpose(data[:len(t), -5 - pol],(0,2,1))/rawcorrection[-5 - pol, np.newaxis,:,:]).tofile('miriadextract_' + polstr[pol] + '_' + ano)
	else:
		np.transpose(data[:len(t), -5 - pol],(0,2,1)).tofile('miriadextract_' + polstr[pol] + '_' + ano)	

#np.savetxt('miriadextract_' + ano + "_sunpos.dat", sunpos[:len(t)], fmt='%8.5f')	

f = open('miriadextract_' + ano + "_timing.dat",'w')
for time in timing:
	f.write("%s\n"%time)
f.close()

#data = 0

for pol in [-5, -6]:
	command = "./omnical " + 'miriadextract_' + polstr[pol] + '_' + ano + " " + infopath[-5 - pol] + " " + str(len(t)) + " " + str(nfreq) + " "  + str(nant)
	print FILENAME + " MSG: System call: ",  command
	os.system(command)
	print np.fromfile('miriadextract_' + polstr[pol] + '_' + ano + ".omnical", dtype = 'float32').reshape((len(t), nfreq, 3+2*(info[-5 - pol]['nAntenna']+info[-5 - pol]['nUBL'])))[:5,50,:3]






