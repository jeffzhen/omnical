import numpy as np
import omnical.calibration_omni as omni
import matplotlib.pyplot as plt
import aipy as a
import glob


uvfiles = glob.glob('/data4/paper/2012EoR/psa_live/psa6266/*uvcRREcAC')
aa = a.cal.get_aa('psa6240_v003', np.linspace(.1,.2,203))
lsts = []

for uvfile in uvfiles:
    uv = a.miriad.UV(uvfile)
    for (uvw,t,bl), d in uv.all():
        aa.set_jultime(t)
        lsts.append(aa.sidereal_time())

print lsts

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
noise = f['tsys_jy'] **2
im = plt.imshow(noise)
im.set_clim(0,1e11)
plt.colorbar()
plt.show()
f.close()
