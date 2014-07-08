import aipy as ap
import numpy as np
import commands, os, time, math, ephem, sys
import calibration_omni as omni
FILENAME = "omnical.py"

#########
#Simple script that apply a set of .omnical calibration parameters to a set of uv files to create a set of new uv files, ASSUMING they contain the same time stamps for the # of time stamps in input uv files. .omnical could have more stamps in the end. 

##############Config parameters###################################
latP = -0.53619181096511903
lonP = 0.37399448506783717

calparano = 'lst_v007_fg_crosstalkrm'
ano = 'lst_v006_I_w_v007par_crosstalkrm'##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
uvfiles = commands.getoutput('ls /data4/paper/arp/lst_v006_I/*.uv -d').split()[:15]
wantpols = {'SI':1}

infopaths = {'SI':'./redundantinfo_PSA32.txt'}
oppath = './lst_v006_I_results/'
calparpath = './lst_v007_results/'

removedegen = 1

needrawcal = False #if true, (generally true for raw data) you need to take care of having raw calibration parameters in float32 binary format freq x nant
rawpaths = {'xx':"testrawphasecalparrad_xx", 'yy':"testrawphasecalparrad_yy"}

keep_binary_data = False
########Massage user parameters################
oppath += '/' 

####read redundant info################
info = [omni.read_redundantinfo(infopaths[key]) for key in wantpols.keys()]
#print info[0]['bl1dmatrix']
#exit(1)

####get some info from the first uvfile   ################
uv=ap.miriad.UV(uvfiles[0])
nfreq = uv.nchan;
nant = uv['nants'] / 2 # 'nants' counting ant-pols, so divide 2
startfreq = uv['sfreq']
dfreq = uv['sdf']
del(uv)

###read raw phase calibration prameters over frequencyfor each antenna, 203 by 32 in radiants; this part can be replaced################
if needrawcal:
	rawcalpar = np.asarray([np.fromfile(rawpaths[key], dtype="complex64").reshape(nfreq, nant) for key in wantpols.keys()])
	rawcorrection = np.zeros((len(wantpols), nfreq, nant*(nant+1)/2), dtype='complex64') + 1#to be dividing the data;  data/calpar = model
	for p in range(len(wantpols)):
		for i, bl in zip(range(len(info[p]['subsetbl'])), info[p]['subsetbl']):
			a1, a2 = info[p]['bl2d'][i]
			rawcorrection[p, :, bl] = np.conj(rawcalpar[p, :,a1]) *  rawcalpar[p, :,a2]



###reorder and dump the binary data from miriad################
if not os.path.exists(oppath):
	os.makedirs(oppath)



#####apply calpar and create new uv##################################
calparfilenames = [calparpath + 'miriadextract_' + pol + '_' + calparano + '.omnical' for p, pol in zip(range(len(wantpols)), wantpols.keys())]
print FILENAME + " MSG: Applying", calparfilenames,
sys.stdout.flush()
omni.apply_omnical_uvs(uvfiles, calparfilenames, info, wantpols, oppath, ano)

print "Done!"



