import aipy as ap
import numpy as np
import ephem
import commands, os
import time
import calibration_omni as omni
latP = -0.53619181096511903
lonP = 0.37399448506783717
sa = ephem.Observer()
sa.lon = lonP
sa.lat = latP
sun = ephem.Sun()
julDelta = 2415020 # =julian date - pyephem's Observer date

ano = 'test'
uvfiles = ['zen.2455903.00836.uvcRRE']
print len(uvfiles), "uv files to be processed"

####read redundant info
print "Reading redundant info...", 
info = omni.read_redundantinfo('./redundantinfo_PSA32.txt')
print "done", len(info['subsetant']), info['nUBL'], info['nbl']
#print info['bl2d']
#exit(1)

####get some info from the first uvfile
uv=ap.miriad.UV(uvfiles[0])
nfreq = uv.nchan;
nant = uv['nants'] / 2 # 'nants' counting ant-pols, so divide 2
startfreq = uv['sfreq']
dfreq = uv['sdf']

###read raw phase calibration prameters over frequencyfor each antenna, 203 by 32 in radiants; this part can be replaced
rawcalpar = np.exp(1.j * np.fromfile("psa930rawphasecalparrad", dtype="float32").reshape(nfreq, nant))
rawcorrection = np.zeros((2, nfreq, nant*(nant+1)/2), dtype='complex64')#to be dividing the data;  data/calpar = model
bl = 0;
for a2 in range(nant):
	for a1 in range(a2 + 1):
		for pol in [-5,-6]:
			rawcorrection[pol + 6, :, bl] = np.conj(rawcalpar[:,a1]) *  rawcalpar[:,a2]
		bl += 1


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
		print timing[-1]#uv.nchan
	#print uv['nants']
	#cnter = 0
	currentpol = 0
	for preamble, rawd in uv.all():
		#if cnter >1000:
		#	exit(1)
		#print uv['pol'], preamble[1]
		#cnter = cnter +1
		if len(t) < 1 or t[-1] != preamble[1]:#first bl of a timeslice
			t += [preamble[1]]
			sa.date = preamble[1] - julDelta
			#sun.compute(sa)
			timing += [sa.date.__str__()]
			if len(t) > len(data):
				print "expanding number of time slices from", len(data), "to", len(data) + deftime
				data = np.concatenate((data, np.zeros((deftime, 2, nant * (nant + 1) / 2, nfreq), dtype = 'complex64'))) 
				#sunpos = np.concatenate((sunpos, np.zeros((deftime, 2))))
				#sunpos[len(t) - 1] = np.asarray([[sun.alt, sun.az]])
		if -5 == uv['pol'] or -6 == uv['pol']:
			if currentpol != uv['pol']:
				bl = 0
				currentpol = uv['pol']
			data[len(t) - 1, currentpol + 6, bl] = rawd.data.astype('complex64')
			bl += 1		

print len(t), "slices read."

#print (np.transpose(data[:len(t), -5 + 6],(0,2,1)))[5,50,:20]
#print (np.transpose(data[:len(t), -5 + 6],(0,2,1))/rawcorrection[-5+6,np.newaxis,:,:])[5,50,:20]
#exit(1)
polstr = ['yx', 'xy', 'yy', 'xx', 'bla', 'bla', 'bla', 'bla']
for pol in [-5, -6]:
	print data[:len(t), pol + 6].shape
	(np.transpose(data[:len(t), pol + 6],(0,2,1))/rawcorrection[pol + 6,np.newaxis,:,:]).tofile('miriadextract_' + polstr[pol] + '_' + ano)	

#np.savetxt('miriadextract_' + polstr[pol] + '_' + ano + "_sunpos.dat", sunpos[:len(t)], fmt='%8.5f')	

f = open('miriadextract_' + ano + "_timing.dat",'w')
for time in timing:
	f.write("%s\n"%time)
f.close()

for pol in [-5, -6]:
	command = "./omnical " + 'miriadextract_' + polstr[pol] + '_' + ano + " redundantinfo_PSA32.txt " + str(len(t)) + " " + str(nfreq) + " "  + str(nant)
	print command
	os.system(command)
	exit(1)




