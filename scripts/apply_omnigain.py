#!/usr/bin/env python

import aipy as ap
import numpy as np
import commands, os, time, math, ephem, optparse, sys, shutil
import omnical.calibration_omni as omni
import cPickle as pickle
import scipy.signal as ss
import scipy.ndimage.filters as sfil
FILENAME = "apply_omnigain.py"
print "#Omnical Version %s#"%omni.__version__
######################################################################
o = optparse.OptionParser()
o.add_option('-t', '--tag', action = 'store', default = 'PSA128', help = 'tag name of this calibration')

o.add_option('-o', '--outputpath', action = 'store', default = ".", help = 'output folder')
o.add_option('-x', '--xgain', action = 'store', default = "NOTGIVEN", help = 'output folder')
o.add_option('-y', '--ygain', action = 'store', default = "NOTGIVEN", help = 'output folder')

o.add_option('-f', '--overwrite', action = 'store_true', help = 'whether to overwrite if the new uv files already exists')
opts,args = o.parse_args(sys.argv[1:])
uvfiles = [os.path.expanduser(arg) for arg in args]

oppath = os.path.expanduser(opts.outputpath)
if not os.path.isdir(oppath):
	os.makedirs(oppath)

overwrite_uvs = opts.overwrite
ano = opts.tag

if opts.xgain == "NOTGIVEN" and opts.ygain == "NOTGIVEN":
	raise IOError("At least -x and -y  omnigain files are required.")

gainpaths = {}
omnigains = {}
infos = {}
flags = {}
if opts.xgain != "NOTGIVEN":
	gainpaths['x'] = os.path.expanduser(opts.xgain)

if opts.ygain != "NOTGIVEN":
	gainpaths['y'] = os.path.expanduser(opts.ygain)


for pol in gainpaths.keys():
	infos[pol] = omni.read_redundantinfo(gainpaths[pol].replace('.omnigain', '.binfo'))
	omnigains[pol] = omni.load_omnigain(gainpaths[pol], infos[pol])
	flags[pol] = np.fromfile(gainpaths[pol].replace('.omnigain', '.omniflag'), dtype='bool').reshape((len(omnigains[pol]), int(omnigains[pol][0,0,3])))
##########################################################################
omni.apply_omnigain_uvs(uvfiles, omnigains, [[0,0]], infos, oppath, ano, verbose = True, comment = '_'.join(sys.argv), flags = flags, overwrite = overwrite_uvs)
