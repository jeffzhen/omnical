import numpy as np
import omnical.calibration_omni as omni
import matplotlib.pyplot as plt
import aipy as a
import glob, os
from scipy import interpolate

D = 6266#6301 1.9#6290  1.9#6300 1.63#6266 1.9#6242 1.9#
nadd = 7
nchan = 203
uvfiles = sorted(glob.glob('/data4/paper/2012EoR/psa_live/psa%i/*uvcRREcAC'%D))

lst_file = '/data2/home/hz2ug/omnical/results/psa%i_lsts.txt'%D
if os.path.isfile(lst_file):
    lsts = np.loadtxt(lst_file)
else:
    aa = a.cal.get_aa('psa6240_v003', np.linspace(.1,.2,nchan))
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
    lsts = np.array(lsts)

#info=omni.read_redundantinfo('/data2/home/hz2ug/omnical/doc/redundantinfo_PSA64_ba19_37_50.bin')
#chis = [np.fromfile("/data4/paper/2012EoR/psa_live/forlstbinning_omnical_2/data_psa%i_245%i.%i_xx_add7.omnical"%(D, D, i), dtype = 'float32')[2::3+2*(info['nAntenna'] + info['nUBL'])] for i in range(1,7)]

info=omni.read_redundantinfo('/data2/home/hz2ug/omnical/doc/redundantinfo_PSA64_7ba_7bu_08-15-2014.bin')
chis = [np.fromfile("/data4/paper/2012EoR/psa_live/PSA64_omnical_results_Aug_2014/psa%i_245%i.%i_aug2014_xx.omnical"%(D, D, i), dtype = 'float32')[2::3+2*(info['nAntenna'] + info['nUBL'])] for i in range(1,7)]


chis = [chi.reshape((len(chi)/nchan, nchan)) for chi in chis]
chi2 = np.concatenate((chis), axis = 0)
dof = len(info['crossindex']) * 2*nadd / (2*nadd + 1) - (info['nAntenna'] + info['nUBL'] - 2)


print "WARNING: THE TSYS MODEL IS INCORRECT SINCE ITS (XX+YY)/2"
f = np.load('/data2/home/hz2ug/PSA64_plot_for_zaki/tsys_model_jy_corrected.npz')
noise = np.array(f['tsys_jy'] **2) / (43*100./nchan*1e6)

noise_model = interpolate.interp2d(np.arange(0, nchan), f['lsts']/12.*np.pi, noise)
#if min(lsts[1:] - lsts[:-1]) >= 0:
#    modeled_noise = noise_model(range(0, nchan), lsts)
#else:
#    turn_pt = np.argmin(lsts[1:] - lsts[:-1]) + 1
#    modeled_noise = np.concatenate((noise_model(range(0, nchan), lsts[:turn_pt]), noise_model(range(0, nchan), lsts[turn_pt:])), axis = 0)
modeled_noise = np.zeros((len(lsts), nchan))
modeled_noise[np.argsort(lsts)] = noise_model(range(0, nchan), np.sort(lsts))
model_range = (lsts < max(f['lsts']/12.*np.pi)) & (lsts > min(f['lsts']/12.*np.pi))

onedchi = (chi2[model_range, 25:160] / dof / modeled_noise[model_range,25:160]).flatten()
nanmask = (~(np.isnan(onedchi)))&np.isfinite(onedchi)

plt.subplot('141')
im = plt.imshow(chi2[model_range] / dof)
im.set_clim(0,1.5e4)
plt.colorbar()
#plt.show()

plt.subplot('142')
im = plt.imshow(modeled_noise[model_range])
im.set_clim(0,1.5e4/np.median(onedchi[nanmask]))
plt.colorbar()

plt.subplot('143')
im = plt.imshow(chi2[model_range] / dof / modeled_noise[model_range])
np.savez('/data2/home/zakiali/jeff_psa64/chisqdof_%i.npz'%D, chi2[model_range] / dof / modeled_noise[model_range])
im.set_clim(0,3)
plt.colorbar()


plt.subplot('144')
plt.hist(onedchi[nanmask], np.arange(0,3*np.median(onedchi[nanmask]),3*np.median(onedchi[nanmask])/100.))
plt.title('%.2f, %.2f'%(np.mean(onedchi[nanmask]), np.median(onedchi[nanmask])))
plt.show()
f.close()
