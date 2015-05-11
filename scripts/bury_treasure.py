#!/usr/bin/env python

import numpy as np
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
import omnical._omnical as _O
import optparse, sys

FILENAME = "bury_treasure.py"
PI = np.pi
TPI = 2 * np.pi
print "#Omnical Version %s#"%omni.__version__

######################################################################
##############Config parameters###################################
######################################################################
o = optparse.OptionParser()

o.add_option('-i', '--infopath', action = 'store', default = None, help = 'redundantinfo file to read')
o.add_option('-p', '--pol', action = 'store', default = None, help = 'polarizations to store')

opts,args = o.parse_args(sys.argv[1:])

if len(args) != 3:
	raise IOError("Script requires 3 arguments: treasure path, number of lsts, and number of frequencies.")

treasure_path = args[0]
nLST = int(args[1])
nFrequency = int(args[2])

if os.path.isdir(os.path.expanduser(treasure_path)):
	raise IOError("Treasure path %s exists."%treasure_path)

pols = opts.pol
if pols is None:
	raise TypeError("polarization not set. Use -p to set polarizations like -p xx,yy")
else:
	pols = pols.split(',')

if np.max([len(pol) for pol in pols]) > 9:
	raise TypeError("polarization string cannot exceed 9 characters. Use -p to set polarizations like -p xx,yy")

treasure = omni.Treasure(treasure_path, nlst = nLST, nfreq = nFrequency)
try:
	if opts.infopath is not None:
		info = omni.read_redundantinfo(os.path.expanduser(opts.infopath))
		for p, pol in enumerate(pols):
			for i, ublvec in enumerate(info['ubl']):
				print pol, ublvec, "%i/%i"%(p * len(info['ubl']) + i, len(pols) * len(info['ubl']))
				sys.stdout.flush()
				treasure.add_coin((pol, ublvec))
		print "Done."
except KeyboardInterrupt:
	treasure.burn()
	print "Aborted due to KeyboardInterrupt."
