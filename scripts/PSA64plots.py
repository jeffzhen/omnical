import numpy as np
import omnical.calibration_omni as omni
import matplotlib.pyplot as plt
import aipy as a
import glob, os
from scipy import interpolate

D = 6266
uvfiles = glob.glob('/data4/paper/2012EoR/psa_live/psa%i/*uvcRREcAC'%D)

lst_file = '/data2/home/hz2ug/omnical/results/psa%i_lsts.txt'%D
if os.path.isfile(lst_file):
    lsts = np.loadtxt(lst_file)
else:
    aa = a.cal.get_aa('psa6240_v003', np.linspace(.1,.2,203))
    lsts = []
    tlast = None
    for uvfile in uvfiles:
        print "Processing", uvfile
        uv = a.miriad.UV(uvfile)
        for (uvw,t,bl), d in uv.all():
            if t != tlast:
                tlast = t
                aa.set_jultime(t)
                lsts.append(aa.sidereal_time())
    np.savetxt(lst_file, np.array(lsts))

info=omni.read_redundantinfo('/data2/home/hz2ug/omnical/doc/redundantinfo_PSA64_ba19_37_50.bin')


chis = [np.fromfile("/data4/paper/2012EoR/psa_live/forlstbinning_omnical_2/data_psa6266_2456266.%i_xx_add7.omnical"%i, dtype = 'float32')[2::3+2*(info['nAntenna'] + info['nUBL'])] for i in range(1,7)]
chis = [chi.reshape((len(chi)/203, 203)) for chi in chis]
chi2 = np.concatenate((chis), axis = 0)
dof = len(info['crossindex']) - (info['nAntenna'] + info['nUBL'] - 2)

im = plt.imshow(chi2 / dof)
im.set_clim(0,5e4)
plt.colorbar()
plt.show()

f = np.load('/data2/home/hz2ug/omnical/results/tsys_model_jy.npz')
noise = np.array(f['tsys_jy'] **2)

noise_model = interpolate.interp2d(np.arange(0, 203), f['lsts'], noise)
if min(lsts[1:] - lsts[:-1]) >= 0:
    modeled_noise = noise_model(range(0, 203), lsts)
else:
    turn_pt = np.argmin(lsts[1:] - lsts[:-1]) + 1
    modeled_noise = np.concatenate((noise_model(range(0, 203), lsts[:turn_pt]), noise_model(range(0, 203), lsts[turn_pt:])), axis = 1)

im = plt.imshow(modeled_noise)
im.set_clim(0,1e11)
plt.colorbar()
plt.show()
f.close()
