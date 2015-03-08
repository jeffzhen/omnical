import numpy as np
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
import omnical._omnical as _O
import optparse, sys
import matplotlib.pyplot as plt

treasure = omni.Treasure(sys.argv[1])
item = sys.argv[2]
if item not in ['count', 'variance_re', 'variance_im', 'weighted_mean', 'mean', 'weighted_variance']:
	print "%s not recognized."%item
for p, pol in enumerate(treasure.ubls.keys()):
	c = treasure.get_coin_now((pol,treasure.ubls[pol][np.argsort(np.linalg.norm(treasure.ubls[pol], axis=1))[0]]))

	plt.subplot('1%i%i'%(len(treasure.ubls.keys()), p+1))
	plt.imshow(c.__getattr__(item), aspect = 1/5., interpolation='none');plt.colorbar();plt.title(pol)

plt.show()
