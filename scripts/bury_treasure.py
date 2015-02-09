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

opts,args = o.parse_args(sys.argv[1:])

if len(args) != 3:
	raise IOError("Script requires 3 arguments: treasure path, number of lsts, and number of frequencies.")

treasure_path = args[0]
nLST = int(args[1])
nFrequency = int(args[2])

if os.path.isdir(os.path.expanduser(treasure_path)):
	raise IOError("Treasure path %s exists."%treasure_path)

treasure = omni.Treasure(treasure_path, nlst = nLST, nfreq = nFrequency)

try:
	if opts.infopath is not None:
		info = omni.read_redundantinfo(os.path.expanduser(opts.infopath))
		for i, ublvec in enumerate(info['ubl']):
			print ublvec, "%i/%i"%(i, len(info['ubl']))
			sys.stdout.flush()
			treasure.add_coin(ublvec)
		print "Done."
except KeyboardInterrupt:
	treasure.burn()
	print "Aborted due to KeyboardInterrupt."
