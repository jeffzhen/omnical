#!/usr/bin/env python

import numpy as np
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
import omnical._omnical as _O
import optparse, sys
import matplotlib.pyplot as plt

print "#Omnical Version %s#"%omni.__version__
PI = np.pi
TPI = 2 * np.pi
######################################################################
##############Config parameters###################################
######################################################################
o = optparse.OptionParser()

o.add_option('--max', action = 'store', type = 'float', default = None, help = 'Max value of plotting.')
o.add_option('--min', action = 'store', type = 'float', default = None, help = 'Min value of plotting.')

opts,args = o.parse_args(sys.argv[1:])

treasure = omni.Treasure(args[0])
item = args[1]
if item not in ['count', 'variance_re', 'variance_im', 'weighted_mean', 'mean', 'weighted_variance']:
	print "%s not recognized."%item
	exit()
if item in ['mean', 'weighted_mean']:
	if len(args) < 3 or args[2] not in ['amp', 'phs']:
		print "Please specify amp or phs"
		exit()

for p, pol in enumerate(treasure.ubls.keys()):
	c = treasure.get_coin_now((pol,treasure.ubls[pol][np.argsort(np.linalg.norm(treasure.ubls[pol], axis=1))[0]]))

	plt.subplot('1%i%i'%(len(treasure.ubls.keys()), p+1))
	data = c.__getattr__(item)
	if item != 'count':
		flag = c.count==0
		data[flag] = np.nan
	if item in ['mean', 'weighted_mean']:
		if args[2] == 'amp':
			data = np.abs(data)
		else:
			data = np.angle(data)
	flag = np.isnan(data)|np.isinf(data)
	if opts.max is None:
		vmax = np.percentile(data[~flag], 95)
	else:
		vmax = opts.max
	if opts.min is None:
		vmin = np.percentile(data[~flag], 5)
	else:
		vmin = opts.min
	plt.imshow(data, aspect = 1/5., interpolation='none', vmin=vmin, vmax=vmax);plt.colorbar();plt.title(pol)

plt.show()
