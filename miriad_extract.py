import aipy as ap
import numpy as np
import ephem
import commands
import time
latP = -0.53619181096511903
lonP = 0.37399448506783717
sa = ephem.Observer()
sa.lon = lonP
sa.lat = latP
sun = ephem.Sun()
julDelta = 2415020 # =julian date - pyephem's Observer date

ano = 'psa903'
uvfiles = commands.getoutput('ls /data4/raw_data/psa903-930/' + ano + '/*.uvcRRE -d').split()
print len(uvfiles)

#pol = -5 -6 -7 -8 are xx yy xy yx

#print files
uv=ap.miriad.UV(uvfiles[0])
nfreq = uv.nchan;
nant = uv['nants'] / 2 # 'nants' counting ant-pols, so divide 2
startfreq = uv['sfreq']
dfreq = uv['sdf']

deftime = 2000
data = np.zeros((deftime, 2, nant * (nant + 1) / 2, nfreq), dtype = 'complex64')
sunpos = np.zeros((deftime, 2))
t = []
timing = []

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
			sun.compute(sa)
			timing += [sa.date.__str__()]
			if len(t) > len(data):
				print "expanding number of time slices from", len(data), "to", len(data) + deftime
				data = np.concatenate((data, np.zeros((deftime, 2, nant * (nant + 1) / 2, nfreq), dtype = 'complex64'))) 
				sunpos = np.concatenate((sunpos, np.zeros((deftime, 2))))
				sunpos[len(t) - 1] = np.asarray([[sun.alt, sun.az]])
		if -5 == uv['pol'] or -6 == uv['pol']:
			if currentpol != uv['pol']:
				bl = 0
				currentpol = uv['pol']
			data[len(t) - 1, currentpol + 6, bl] = rawd.data.astype('complex64')
			bl += 1		

print len(t), "slices read."

for pol in [-5, -6]:
	polstr = ['yx', 'xy', 'yy', 'xx', 'bla', 'bla', 'bla', 'bla']
	data[:len(t), pol + 6].tofile('miriadextract_' + polstr[pol] + '_' + ano)	
np.savetxt('miriadextract_' + polstr[pol] + '_' + ano + "_sunpos.dat", sunpos[:len(t)], fmt='%8.5f')	

f = open('miriadextract_' + polstr[pol] + '_' + ano + "_timing.dat",'w')
for time in timing:

	f.write("%s\n"%time)

f.close()

		





