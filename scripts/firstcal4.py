#!/usr/bin/env python

import aipy as ap
import numpy as np
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
import omnical._omnical as _O
import optparse, sys
import scipy.signal as ss
import scipy.linalg as la
from scipy.stats import nanmedian
import matplotlib.pyplot as plt
import cPickle as pickle

FILENAME = "firstcal4.py"
PI = np.pi
TPI = 2 * np.pi
print "#Omnical Version %s#"%omni.__version__
def gain_curve(mhz):
    ghz = np.array(mhz) / 1000.
    return np.inner([-167333390752.98276, 198581623581.65594, -102487141227.4993, 30027423590.548084, -5459067124.669095, 630132740.98792362, -45056600.848056234, 1822654.0034047314, -31892.9279846797], [ghz**i for i in range(8,-1,-1)])

######################################################################
##############Config parameters###################################
######################################################################
o = optparse.OptionParser()

#ap.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('-C', '--cal', action = 'store', default = None, help = 'calfile for processing uv files')
o.add_option('--max', action = 'store', type = 'int', default = 5, help = 'Max number of iterations when removing bad antennas.')
o.add_option('--healthbar', action = 'store', type = 'float', default = 2, help = 'Health threshold (0-100) over which an antenna is marked bad. 2 by default.')
o.add_option('--suppress', action = 'store', type = 'float', default = 1, help = 'Amplitude of the gains for the bad antennas. Larger means more suppressed. Only useful if you dont want to use flagging.')
o.add_option('-f', '--freq_range', action = 'store', default = '0_0', help = 'Frequency bin number range to use for fitting amp and delay seperated by underscore. 0_0 by default and will process all frequencies.')
o.add_option('-o', '--outputpath', action = 'store', default = "DONT_WRITE", help = 'Output folder. No output by default.')
o.add_option('-t', '--info_tag', action = 'store', default = "DEFAULT", help = 'Name tag for output redundantinfo file.')
o.add_option('--overwrite', action = 'store_true', help = 'whether to overwrite if the new uv files already exists')
o.add_option('--ampdelay', action = 'store_true', help = 'whether to print out amplitude and delay for Calfile usage. Amplitude print out disabled because theres no justifiable reason to ever use it.')
o.add_option('--plot', action = 'store_true', help = 'Whether to make plots in the end.')
o.add_option('-e', '--tol', action = 'store', type = 'float', default = 1e-2, help = 'tolerance of antenna location deviation when computing unique baselines.')
o.add_option('--ba', action = 'store', default = '', help = 'bad antenna number indices seperated by commas')
o.add_option('--bu', action = 'store', default = '', help = 'bad unique baseline indicated by ant pairs (seperated by .) seperated by commas: 1.2,3.4,10.11')
o.add_option('--flagsigma', action = 'store', type = 'float', default = 4, help = 'Number of sigmas to flag on chi^2 distribution. 4 sigma by default. For full cadence data, 3 is recommended.')



opts,args = o.parse_args(sys.argv[1:])
overwrite = opts.overwrite
make_plots = opts.plot
print_ampdelay = opts.ampdelay
oppath = os.path.expanduser(opts.outputpath)
info_tag = opts.info_tag
uvfiles = [os.path.expanduser(arg) for arg in args]
healthbar = opts.healthbar
bad_ant_suppress = opts.suppress
max_try = opts.max
nsigma = opts.flagsigma

try:
    badAntenna = [int(i) for i in opts.ba.split(',')]
except:
    badAntenna = []
try:
    if opts.bu != '':
        badUBLpair = [[int(j) for j in i.split('.')] for i in opts.bu.split(',')]
    else:
        badUBLpair = []
except:
    badUBLpair = []
redundancy_tol = opts.tol

[fstart,fend] = [int(x) for x in opts.freq_range.split('_')]

wantpols = {}
for p in ['xx', 'yy']:
    wantpols[p] = ap.miriad.str2pol[p]

for uvf in uvfiles:
    if not os.path.isdir(uvf):
        uvfiles.remove(uvf)
        print "WARNING: uv file path %s does not exist!"%uvf

if '.odf' in uvfiles[0]:
    input_type = 'odf'
else:
    input_type = 'uv'

if len(uvfiles) == 0:
    raise IOError("ERROR: No valid uv files detected in input. Exiting!")
elif len(uvfiles) != 4 and len(uvfiles) != 1 and input_type == 'uv':
    raise IOError("ERROR: %i uvfiles are inputed and assumed to be corresponding to different polarizations, but we expect %s polarizations. Exiting!"%(len(uvfiles), 4))
elif len(uvfiles) != 1 and input_type == 'odf':
    raise IOError("ERROR: %i odfs are inputed and only one at a time is supported."%len(uvfiles))
else:
    if input_type == 'uv' and len(uvfiles) == 4:
        uvfiles_dic = {}
        for i, p in enumerate(['xx', 'xy', 'yx', 'yy']):
            uvfiles_dic[p] = uvfiles[i]


if input_type == 'uv':
    print "Reading calfile %s..."%opts.cal,
    sys.stdout.flush()
    aa = ap.cal.get_aa(opts.cal, np.array([.15]))
    print "Done. Antenna layout:"
    print aa.ant_layout
    sys.stdout.flush()


removedegen = True #this is only for amplitude, because for phase there's a special degenaracy removal anyways
removeadditive = False
removeadditiveperiod = -1

need_crude_cal = True

converge_percent = 0.0001
max_iter = 50
step_size = .3

######################################################################
######################################################################
######################################################################

########Massage user parameters###################################
oppath += '/'

####get some info from the first uvfile   ################
print "Getting some basic info from %s"%uvfiles[0],
sys.stdout.flush()

sa = ephem.Observer()
sa.pressure = 0

if input_type == 'odf':
    with open(uvfiles[0]+'/header.txt') as f:
        for l in f.readlines():
            if l.split()[0] == 'nChannels':
                nfreq = int(l.split()[1])
            elif l.split()[0] == 'nAntennas':
                nant = int(l.split()[1])
            elif l.split()[0] == 'longitude':#odf header has bug that switched lat and lon
                sa.lat = float(l.split()[1]) * PI/180.
            elif l.split()[0] == 'latitude':
                sa.lon = float(l.split()[1]) * PI/180.
            elif l.split()[0] == 'startFreq':
                startfreq = float(l.split()[1]) / 1000.
            elif l.split()[0] == 'endFreq':
                endfreq = float(l.split()[1]) / 1000.
        dfreq = (endfreq - startfreq) / nfreq

    header_antloc = np.loadtxt(uvfiles[0].replace('.odf','.odfa') + '/antlocx5.dat')
else:
    uv=ap.miriad.UV(uvfiles[0])
    nfreq = uv.nchan;
    nant = uv['nants']
    sa.lon = aa.lon
    sa.lat = aa.lat
    startfreq = uv['sfreq']#GHz
    dfreq = uv['sdf']#GHz
    del(uv)
print "Done."
sys.stdout.flush()




###start reading miriads################
print FILENAME + " MSG:",  len(uvfiles), "%s files to be processed"%input_type
sys.stdout.flush()

if input_type == 'odf':
    nTimes = []
    for uvfile in uvfiles:
        with open(uvfile+'/header.txt') as f:
            for l in f.readlines():
                if l.split()[0] == 'nIntegrations':
                    nTimes = nTimes + [int(l.split()[1])]
    pol_select = []
    for key in wantpols.keys():
        if 'xx' == key:
            pol_select = pol_select + [0]
        elif 'xy' == key:
            pol_select = pol_select + [1]
        elif 'yx' == key:
            pol_select = pol_select + [2]
        elif 'yy' == key:
            pol_select = pol_select + [3]
    rawdata = np.zeros((len(pol_select), np.sum(nTimes), nfreq, nant * (nant + 1) / 2), dtype='complex64')
    timing = []
    lst = []
    for i, uvfile in enumerate(uvfiles):
        rawdata[:, len(timing):len(timing)+nTimes[i]] = np.fromfile(uvfile+'/visibilities', dtype='complex64').reshape((4, nTimes[i], nfreq, nant * (nant + 1) / 2))[pol_select]
        with open(uvfile.replace('.odf', '.odfa') + '/timing.dat') as f:
            for l in f.readlines():
                sa.date = l.replace('\n','')
                sa.date += ephem.second * 3600 * 4
                lst += [sa.sidereal_time()]
                timing += [sa.date.__str__()]
else:
    if len(uvfiles) == 1:
        rawdata, t, timing, lst, rawflag = omni.importuvs(uvfiles, wantpols, totalVisibilityId = np.concatenate([[[i,j] for i in range(j + 1)] for j in range(len(aa))]), timingTolerance=100)#, nTotalAntenna = len(aa))
    else:

        for p, pol in enumerate(wantpols.keys()):
            if p == 0:
                rawdata, t, timing, lst, rawflag = omni.importuvs([uvfiles_dic[pol]], {pol:wantpols[pol]}, totalVisibilityId = np.concatenate([[[i,j] for i in range(j + 1)] for j in range(len(aa))]), timingTolerance=100)
            else:
                tmpdata, t, timing, lst, tmpflag = omni.importuvs([uvfiles_dic[pol]], {pol:wantpols[pol]}, totalVisibilityId = np.concatenate([[[i,j] for i in range(j + 1)] for j in range(len(aa))]), timingTolerance=100)
                rawdata = np.concatenate((rawdata, tmpdata))
    print FILENAME + " MSG:",  len(t), "slices read."
    sys.stdout.flush()

sun = ephem.Sun()
sunpos  = np.zeros((len(timing), 2))
southern_points = {'hyd':{'ra': '09:18:05.7', 'dec': '-12:05:44'},
'cen':{'ra': '13:25:27.6', 'dec': '-43:01:09'},
'cyg':{'ra': '19:59:28.3', 'dec': '40:44:02'},
'pic':{'ra': '05:19:49.7', 'dec': '-45:46:44'},
'vir':{'ra': '12:30:49.4', 'dec': '12:23:28'},
'for':{'ra': '03:22:41.7', 'dec': '-37:12:30'},
'sag':{'ra': '17:45:40.045', 'dec': '-29:0:27.9'},
'cas':{'ra': '23:23:26', 'dec': '58:48:00'},
'crab':{'ra': '5:34:31.97', 'dec': '22:00:52.1'}}
for source in southern_points.keys():
    southern_points[source]['body'] = ephem.FixedBody()
    southern_points[source]['body']._ra = southern_points[source]['ra']
    southern_points[source]['body']._dec = southern_points[source]['dec']
    southern_points[source]['pos'] = np.zeros((len(timing), 2))
for nt,tm in zip(range(len(timing)),timing):
    sa.date = tm
    sun.compute(sa)
    sunpos[nt] = sun.alt, sun.az
    for source in southern_points.keys():
        southern_points[source]['body'].compute(sa)
        southern_points[source]['pos'][nt] = southern_points[source]['body'].alt, southern_points[source]['body'].az
print FILENAME + " MSG:"
print "data time range UTC: %s to %s"%(timing[0], timing[-1])
print "sun altaz from (%f,%f) to (%f,%f)"%(sunpos[0,0], sunpos[0,1], sunpos[-1,0], sunpos[-1,1])
for source in southern_points.keys():
    print "%s altaz from (%f,%f) to (%f,%f)"%(source, southern_points[source]['pos'][0,0], southern_points[source]['pos'][0,1], southern_points[source]['pos'][-1,0], southern_points[source]['pos'][-1,1])
sys.stdout.flush()
####create redundant calibrators################
new_bad_ant = ["Just to get while loop started"]
trials = 0
calibrators = {}
data = {}
flags = {}
while new_bad_ant != [] and trials < max_try and len(badAntenna + new_bad_ant) < nant / 2:
    trials = trials + 1
    if trials > 1:
        print "##########################################################################"
        print FILENAME + " trial #%i: Recalculating redundant info removing new bad antennas..."%trials, new_bad_ant
        sys.stdout.flush()


    ant_bad_meter = {}
    ubl_bad_meter = {}
    crude_calpar = {}
    if trials > 1:
        #figure out old and new bad antennas
        badAntenna = list(np.sort(badAntenna + new_bad_ant))
        print 'Current bad Antennas (%i):'%len(badAntenna), badAntenna
        print 'Bad unique baselines (%i):'%len(badUBLpair), badUBLpair

    for p, pol in enumerate(wantpols.keys()):


        if input_type == 'uv':
            calibrators[pol] = omni.RedundantCalibrator_PAPER(aa)
        elif input_type == 'odf':
            calibrators[pol] = omni.RedundantCalibrator_X5(header_antloc)

        calibrators[pol].nTime = len(timing)
        calibrators[pol].nFrequency = nfreq
        calibrators[pol].removeDegeneracy = removedegen
        calibrators[pol].convergePercent = converge_percent
        calibrators[pol].maxIteration = max_iter
        calibrators[pol].stepSize = step_size


        timer = time.time()
        calibrators[pol].compute_redundantinfo(badAntenna = badAntenna, badUBLpair = badUBLpair, antennaLocationTolerance = redundancy_tol)
        print "Redundant info on %s computed in %f minutes."%(pol, (time.time() - timer)/60.)
        info = calibrators[pol].Info.get_info()

        ###prepare rawCalpar for each calibrator and consider, if needed, raw calibration################
        if need_crude_cal:
            initant, solution_path, additional_solution_path, degen, _ = omni.find_solution_path(info, tol = calibrators[pol].antennaLocationTolerance, verbose = False)
            crude_calpar[pol] = np.array([omni.raw_calibrate(rawdata[p, 0, f], info, initant, solution_path, additional_solution_path, degen) for f in range(calibrators[pol].nFrequency)])
            data[p] = omni.apply_calpar(rawdata[p], crude_calpar[pol], calibrators[pol].totalVisibilityId)
        else:
            data[p] = rawdata[p]
        #calibrators[pol].rawCalpar = np.zeros((calibrators[pol].nTime, calibrators[pol].nFrequency, 3 + 2 * (calibrators[pol].Info.nAntenna + calibrators[pol].Info.nUBL)),dtype='float32')

        ################################
        ########calibrate################
        ################################

        ################first round of calibration  #########################
        print FILENAME + " MSG: starting calibration on %s. nTime = %i, nFrequency = %i ..."%(pol, calibrators[pol].nTime, calibrators[pol].nFrequency),
        sys.stdout.flush()
        timer = time.time()
        additivein = np.zeros_like(data[p])
        calibrators[pol].logcal(data[p], additivein, verbose=True)
        additiveout = calibrators[pol].lincal(data[p], additivein, verbose=True)
        print "Done. %fmin"%(float(time.time()-timer)/60.)
        sys.stdout.flush()

        flags[pol] = calibrators[pol].flag(nsigma = nsigma, fwindow=3, twindow=99999)
        if input_type == 'odf' and .137 > startfreq and .137 < endfreq:
            flags[pol][:, int((.137-startfreq)/dfreq):int((.138-startfreq)/dfreq)] = True
        #################degeneracy removal on 3 antennas [initant, degen[0], degen[1]]######################################
        print FILENAME + " MSG: Performing additional degeneracy removal on %s"%pol,
        sys.stdout.flush()
        timer = time.time()
        A = np.zeros((calibrators[pol].Info.nAntenna, 4))
        masker = np.zeros((calibrators[pol].Info.nAntenna, calibrators[pol].Info.nAntenna))
        for a in [initant, degen[0], degen[1]]:
            A[a] = list(calibrators[pol].Info.antloc[a]) + [1.]
            masker[a,a] = 1.
        matrix = np.identity(calibrators[pol].Info.nAntenna) - np.array([list(vec) + [1.] for vec in info['antloc']]).dot(la.pinv(A.transpose().dot(A)).dot(A.transpose()).dot(masker))

        calibrators[pol].rawCalpar[:, :, (3 + calibrators[pol].Info.nAntenna):(3 + 2 * calibrators[pol].Info.nAntenna)] = ((matrix.dot(calibrators[pol].rawCalpar[:, :, (3 + calibrators[pol].Info.nAntenna):(3 + 2 * calibrators[pol].Info.nAntenna)].transpose(0,2,1)).transpose(1,2,0)) + PI)%TPI - PI
        print "Done. %fmin"%(float(time.time()-timer)/60.)
        sys.stdout.flush()

        #######################diagnose###############################
        ant_bad_meter[pol], ubl_bad_meter[pol] = calibrators[pol].diagnose(data = data[p], additiveout = additiveout, flag = flags[pol], healthbar = healthbar, verbose = False)
        nbad = 0
        for ab in ant_bad_meter[pol]:
            if ab > healthbar:
                nbad += 1
        if nbad > 0:
            print FILENAME + " MSG: %i bad antennas found on %s:"%(nbad, pol),
            for i, ab in enumerate(ant_bad_meter[pol]):
                if ab > healthbar:
                    print "%i:%i"%(calibrators[pol].Info.subsetant[i],ab),
            print ""
        else:
            print FILENAME + " MSG: %i bad antennas found on %s"%(nbad, pol)
        sys.stdout.flush()


    new_bad_ant = []
    for a in range(calibrators[wantpols.keys()[0]].Info.nAntenna):
        for pol in wantpols.keys():
            if ant_bad_meter[pol][a] > healthbar:
                new_bad_ant.append(calibrators[wantpols.keys()[0]].Info.subsetant[a])
                break

    if new_bad_ant == [] or trials == max_try:#start delay fitting
        print "--------------------------------------------------------------------------"
        print FILENAME + " MSG: Starting delay fit"
        worst_delay_ant = {}
        n_bad_delay = {}
        linearcalpar = {}#combine lincal results averaged over time and the initial crude_calpar results
        linearcalpar_model = {}
        solution = {}
        delay = {}
        delay_error = {}
        delay_fit_count = {}
        freq_flag = np.zeros(nfreq, dtype='bool')
        for pol in wantpols.keys():
            linearcalpar[pol] = 10.**(calibrators[pol].rawCalpar[:,:,3:3+calibrators[pol].Info.nAntenna]) * np.exp(1j * calibrators[pol].rawCalpar[:,:,3+calibrators[pol].Info.nAntenna:3+2*calibrators[pol].Info.nAntenna])
            if input_type != 'odf':
                linearcalpar[pol][flags[pol]] = np.nan
            linearcalpar[pol] = nanmedian(np.real(linearcalpar[pol]), axis = 0) + 1j * nanmedian(np.imag(linearcalpar[pol]), axis = 0)
            linearcalpar[pol][np.isnan(linearcalpar[pol])] = 1
            linearcalpar[pol] = linearcalpar[pol] * crude_calpar[pol][:, calibrators[pol].Info.subsetant]

            freq_flag = freq_flag|(np.sum(flags[pol], axis = 0) > .4 * calibrators[pol].nTime)|np.isnan(np.sum(linearcalpar[pol], axis=1))
            freq_flag = freq_flag|(ss.convolve(freq_flag,[.2,.4,1,.4,.2],mode='same')>=1)
            if freq_flag.all():
                raise ValueError("Entire %s data is flagged. Try to make --flagsigma larger to have less aggressive flagging."%pol)
        for pol in wantpols.keys():
            antlocs = np.array(calibrators[pol].antloc)
            antlocs = antlocs - np.mean(antlocs, axis = 0)#move origin to center of array
            A = np.array([list(a) + [1] for a in antlocs])
            AAA = A.dot(np.linalg.pinv(A.transpose().dot(A)).dot(A.transpose()))

            ################################################
            ####delay fitting and fancy degeneracy removal####
            ################################################

            delay_fit_count[pol] = 0
            MAX_DELAY_ITER = 10 #rarely need more than 1 iteration;
            DELAY_ERROR_THRESH = .1 #RMS error in radians for delay fitting
            #solution1_history = []
            #solution0_history = []
            #error_history = []
            #rough_slope_history = []
            #subant49_history = []
            delay[pol] = np.zeros(calibrators[pol].nTotalAnt, dtype='float')
            delay_error[pol] = np.zeros(calibrators[pol].nTotalAnt, dtype='float')+np.inf
            while (delay_fit_count[pol] == 0) or (np.max(delay_error[pol][calibrators[pol].Info.subsetant]) > DELAY_ERROR_THRESH and delay_fit_count[pol] < MAX_DELAY_ITER):
                delay_fit_count[pol] = delay_fit_count[pol] + 1
                if fend <= fstart:
                    fend = nfreq
                sub_freq_flag = freq_flag[fstart:fend]
                avg_angle = np.angle(linearcalpar[pol][fstart:fend][~sub_freq_flag].transpose()).astype('float32')#2D nant x freq

                ####unwrap phase
                valid_freq_index = np.arange(len(freq_flag))[fstart:fend][~sub_freq_flag]
                valid_freq_index = valid_freq_index - valid_freq_index[int(len(valid_freq_index)/2)]#important for stability
                #step = 2
                rough_slopes = np.median(np.concatenate([((avg_angle[:, step:] - avg_angle[:, :-step] + PI)%TPI - PI) / (valid_freq_index[None, step:]-valid_freq_index[None, :-step]) for step in [1,2,3]], axis = 1), axis = 1)
                line_model = np.outer(rough_slopes, valid_freq_index)
                avg_angle = avg_angle + TPI * np.round((line_model - avg_angle) / TPI)

                #make sure no line model is close to PI away from the avg_angle, which can cause instability of fits to jump from -PI to PI
                unstable_ants = np.std((line_model - avg_angle + PI) % TPI - PI, axis = 1) > 0.5 #this threshold is somewhat arbitrary
                line_model[unstable_ants] = line_model[unstable_ants] + PI
                avg_angle = avg_angle + TPI * np.round((line_model - avg_angle) / TPI)

                ##I am allowing a phase offset per antenna. iterating on degeneracy removal does not seem to help
                A = np.ones((fend-fstart, 2),dtype='float32')
                A[:, 0] = np.arange(fstart, fend) * dfreq + startfreq
                A = A[~sub_freq_flag]
                #######
                #matrix_f0 = (la.pinv(A.transpose().dot(A)).dot(A.transpose()))[1]
                #offsets = matrix_f0.dot(avg_angle.transpose())

                #intersect = np.round(offsets/TPI) * TPI #find closest multiple of 2pi for the intersect
                #overall_angle = 0#omni.meanAngle(offsets, weights = 1/(intersect+1))
                #avg_angle = avg_angle - intersect[:,None] - overall_angle

                ###now fit
                #A = np.arange(fstart, fend) * dfreq + startfreq
                #A = A[~sub_freq_flag].reshape((np.sum(~sub_freq_flag), 1))
                #######
                matrix = (la.pinv(A.transpose().dot(A)).dot(A.transpose()))

                solution[pol] = matrix.dot(avg_angle.transpose())


                ##use delay fit to further correct rephasing degeneracies in linearcalpar
                avg_angle_error = np.angle(linearcalpar[pol]) - np.outer(np.arange(nfreq) * dfreq + startfreq, solution[pol][0]) - solution[pol][1]
                avg_angle_error = (avg_angle_error + PI)%TPI - PI
                fit_degen = avg_angle_error.dot(AAA.transpose())
                fit_degen[np.isnan(fit_degen)] = 0
                linearcalpar[pol] = linearcalpar[pol] / np.exp(1.j*fit_degen)
                avg_angle_error = (avg_angle_error - fit_degen + PI)%TPI - PI

                #remove degeneracy in solution
                solution_degen = np.zeros_like(solution[pol])
                solution_degen[1] = omni.medianAngle(solution[pol][1])
                solution[pol][1] = (solution[pol][1] - solution_degen[1] + PI) % TPI - PI
                for sdi in range(5):
                    solution_degen = solution_degen + solution[pol].dot(AAA.transpose())
                    solution[pol] = solution[pol] - solution[pol].dot(AAA.transpose())
                    solution[pol][1] = (solution[pol][1] + PI) % TPI - PI
                    delay[pol][calibrators[pol].Info.subsetant] = solution[pol][0] / TPI
                    #print la.norm(solution[pol][1])
                linearcalpar[pol] = linearcalpar[pol] / np.exp(1.j* (np.outer(np.arange(nfreq) * dfreq + startfreq, solution_degen[0]) + solution_degen[1]))

                #for flagged frequencies, use NaN#the phase model
                model_phase = np.array([(np.arange(nfreq)*dfreq + startfreq) * solution[pol][0, a] + solution[pol][1, a] for a in range(linearcalpar[pol].shape[1])]).transpose()
                if input_type != 'odf':
                    linearcalpar[pol][freq_flag] = np.nan#np.exp(1.j * model_phase[freq_flag])
                #linearcalpar[pol][np.isnan(linearcalpar[pol])] = np.exp(1.j * model_phase[np.isnan(linearcalpar[pol])])

                #also create linearcalpar counter part that is purely the model phase and unity amp
                linearcalpar_model[pol] = np.exp(1.j * model_phase)

                #error terms
                delay_error[pol][calibrators[pol].Info.subsetant] = np.linalg.norm(avg_angle_error[fstart:fend][~sub_freq_flag,:], axis = 0) / (np.sum(~sub_freq_flag))**.5
                #error_history = error_history + [delay_error[pol][calibrators[pol].Info.subsetant]]
                #solution1_history = solution1_history + [solution[pol][1]]
                #solution0_history = solution0_history + [solution[pol][0]]
                #rough_slope_history = rough_slope_history + [rough_slopes]
                #subant49_history = subant49_history + [avg_angle[49]]
            #print "Used %i iterations to find the delays"%delay_fit_count[pol]
            if delay_fit_count[pol] == MAX_DELAY_ITER:
                n_bad_delay[pol] = np.sum(delay_error[pol][calibrators[pol].subsetant] >= DELAY_ERROR_THRESH)
                worst_delay_ant[pol] = calibrators[pol].subsetant[np.argsort(delay_error[pol][calibrators[pol].subsetant])[-n_bad_delay[pol]:]][::-1]
                new_bad_ant = new_bad_ant + list(calibrators[pol].subsetant[delay_error[pol][calibrators[pol].subsetant] >= DELAY_ERROR_THRESH])


        print FILENAME + " MSG: %i bad antennas found during delay fitting."%(len(new_bad_ant))
        if new_bad_ant == [] or len(new_bad_ant) > calibrators[wantpols.keys()[0]].nAntenna / 20.:#end case all good or too many delay fitting failures
            new_bad_ant = []
            print "##########################################################################"
            print "Identified possible bad baselines (not automatically excluded):"
            for u in range(calibrators[wantpols.keys()[0]].nUBL):#assuming calibrators on all pols are the same ubl config
                for p, pol in enumerate(wantpols.keys()):
                    c = calibrators[pol]
                    if ubl_bad_meter[pol][u] > 1:
                        print np.round(c.ubl[u]), ubl_bad_meter[pol][u], "ant pair %i.%i"%(c.subsetant[int(c.ublindex[u][0][0])], c.subsetant[int(c.ublindex[u][0][1])])
                        break

            for pol in wantpols.keys():
                ######amplitude: a single amplitude over all frequency is only useful for calfile and not really accurate
                ##if print_ampdelay:
                    ##amp = np.ones(calibrators[pol].nTotalAnt, dtype='float') * bad_ant_suppress
                    ##amp[calibrators[pol].Info.subsetant] = nanmedian(np.abs(linearcalpar[pol]), axis = 0)
                    ##print FILENAME + " MSG: amplitude factor on %s as |g|:"%pol
                    ##print '{'
                    ##for a1, a2 in zip(range(len(amp)), amp):
                        ##print "%i: %f, "%(a1,a2)
                    ##print '}'
                    ##sys.stdout.flush()

                #print stuff for calfile
                if print_ampdelay:
                    print FILENAME + " MSG: delay on %s in nanoseconds:"%pol
                    print '{'
                    for a1, a2 in enumerate(delay[pol]):
                        print "%i: %f, "%(a1,a2)
                    print '}'
                    sys.stdout.flush()

                if make_plots:
                    nplot = 8
                    plot_a = np.argsort(delay_error[pol][calibrators[pol].subsetant])[range(0, calibrators[pol].nAntenna - 1, calibrators[pol].nAntenna/min(nplot,calibrators[pol].nAntenna)) + [-1]]#[range(nplot/2)+range(-nplot/2,0)]
                    for i, a in enumerate(plot_a):
                        plt.subplot(1, len(plot_a), i+1)
                        pltdata = (np.angle(linearcalpar[pol][:, a]) + PI)%TPI - PI
                        pltdata[freq_flag] = np.nan
                        pltdata_bad = (np.angle(linearcalpar[pol][:, a]) + PI)%TPI - PI
                        pltdata_bad[~freq_flag] = np.nan
                        plt.plot(np.arange(nfreq), pltdata)
                        plt.plot(np.arange(nfreq), ((np.arange(nfreq)*dfreq + startfreq) * solution[pol][0, a] + solution[pol][1, a] + PI)%TPI - PI)
                        #plt.plot(np.arange(nfreq), pltdata_bad, 'r--')
                        plt.title("Ant#%i %.2f"%(calibrators[pol].subsetant[a], delay_error[pol][calibrators[pol].Info.subsetant[a]]))
                        plt.axis([0, nfreq, -PI, PI])
                        #plt.axes().set_aspect('equal')
                    plt.show()
                    plt.hist(delay_error[pol][calibrators[pol].subsetant], 20)
                    if delay_fit_count[pol] == MAX_DELAY_ITER:
                        plt.title("The %i worst fitting antennas on %s are (worst first):\n %s"%(n_bad_delay[pol], pol, worst_delay_ant[pol]))
                    else:
                        plt.title(pol)
                    plt.show()

####################################################################################
#####################cross-pol calibration############################
####################################################################################
###start reading miriads################
print FILENAME + " MSG: Processing cross pol:"
sys.stdout.flush()

del(rawdata)
crosspols = {}
for p in ['xy', 'yx']:
    crosspols[p] = ap.miriad.str2pol[p]

if input_type == 'odf':
    nTimes = []
    for uvfile in uvfiles:
        with open(uvfile+'/header.txt') as f:
            for l in f.readlines():
                if l.split()[0] == 'nIntegrations':
                    nTimes = nTimes + [int(l.split()[1])]
    pol_select = []
    for key in crosspols.keys():
        if 'xx' == key:
            pol_select = pol_select + [0]
        elif 'xy' == key:
            pol_select = pol_select + [1]
        elif 'yx' == key:
            pol_select = pol_select + [2]
        elif 'yy' == key:
            pol_select = pol_select + [3]
    rawdata = np.zeros((len(pol_select), np.sum(nTimes), nfreq, nant * (nant + 1) / 2), dtype='complex64')
    timing = []
    lst = []

    for i, uvfile in enumerate(uvfiles):
        rawdata[:, len(timing):len(timing)+nTimes[i]] = np.fromfile(uvfile+'/visibilities', dtype='complex64').reshape((4, nTimes[i], nfreq, nant * (nant + 1) / 2))[pol_select]
        with open(uvfile.replace('.odf', '.odfa') + '/timing.dat') as f:
            for l in f.readlines():
                sa.date = l.replace('\n','')
                sa.date += ephem.second * 3600 * 4
                lst += [sa.sidereal_time()]
                timing += [sa.date.__str__()]
else:
    if len(uvfiles) == 1:
        rawdata, t, timing, lst, rawflag = omni.importuvs(uvfiles, crosspols, totalVisibilityId = np.concatenate([[[i,j] for i in range(j + 1)] for j in range(len(aa))]), timingTolerance=100)#, nTotalAntenna = len(aa))
    else:

        for p, pol in enumerate(crosspols.keys()):
            if p == 0:
                rawdata, t, timing, lst, rawflag = omni.importuvs([uvfiles_dic[pol]], {pol:crosspols[pol]}, totalVisibilityId = np.concatenate([[[i,j] for i in range(j + 1)] for j in range(len(aa))]), timingTolerance=100)
            else:
                tmpdata, t, timing, lst, tmpflag = omni.importuvs([uvfiles_dic[pol]], {pol:crosspols[pol]}, totalVisibilityId = np.concatenate([[[i,j] for i in range(j + 1)] for j in range(len(aa))]), timingTolerance=100)
                rawdata = np.concatenate((rawdata, tmpdata))
    print FILENAME + " MSG:",  len(t), "slices read."
    sys.stdout.flush()


c = calibrators['xx']
calpar = {}
for pol in ['xx','yy']:
    calpar[pol[0]] = np.ones((linearcalpar[pol].shape[0], c.nTotalAnt), dtype='complex64')
    calpar[pol[0]][:, c.subsetant] = linearcalpar[pol]

for p, pol in enumerate(crosspols.keys()):
    rawdata[p] = omni.apply_calpar2(rawdata[p], calpar[pol[0]], calpar[pol[1]], c.totalVisibilityId)

#construct A and b for cross-pol
A = []
bindex = []
ublcount_sort = np.argsort(c.ublcount)
for u in range(c.nUBL):
    if c.ublcount[u] >= c.ublcount[ublcount_sort[-min(c.nUBL,9999)]]:
        cbl = c.ublindex[u][:,2].astype(int)
        mask = c.reversed[cbl] == 1
        if mask.any():
            A = A + (c.antloc[c.ublindex[u][mask, 0].astype(int)] - c.antloc[int(c.ublindex[u][mask, 0][0])]).tolist()
            bindex = bindex + [[a, int(c.ublindex[u][mask, 2][0])] for a in c.ublindex[u][mask, 2].astype(int)]#inside crossindex
        if not mask.all():
            A = A + (c.antloc[c.ublindex[u][~mask, 0].astype(int)] - c.antloc[int(c.ublindex[u][~mask, 0][0])]).tolist()
            bindex = bindex + [[a, int(c.ublindex[u][~mask, 2][0])] for a in c.ublindex[u][~mask, 2].astype(int)]#inside crossindex
A = np.array(A)[:, :2]
A = np.concatenate((A, -A), axis=0)
bindex = c.subsetbl[c.crossindex[np.array(bindex).astype(int)]]#inside all bl

b = np.empty((len(A), rawdata.shape[1], rawdata.shape[2]), dtype='float32')
for p, pol in enumerate(crosspols.keys()):
    b[p*len(bindex):(p+1)*len(bindex)] = (np.angle(rawdata[p][..., bindex[:,0]]) - np.angle(rawdata[p][..., bindex[:,1]])).transpose((2,0,1))

crosstimer = omni.Timer()
print "Solving cross calibration....",
sys.stdout.flush()
sol = omni.solve_slope(A, b, 1)
sol[:,flags['xx']] = np.nan
sol[:,flags['yy']] = np.nan
median_sol = np.array([nanmedian(sol[i]) for i in range(len(sol))])
print "Done.",
sys.stdout.flush()
crosstimer.tick()
if make_plots:
    for i in range(len(sol)):
        plt.subplot('21%i'%(i+1))
        min_period = TPI / min(abs(A[:,i][abs(A[:,i]) > 1]))
        plt.imshow((sol[i]+min_period/2.)%min_period-min_period/2.);plt.colorbar()
    plt.show()
#apply to data
crosscalpar = {}
#for pol in ['xx','yy']:
    #crosscalpar[pol[0]] = np.ones(c.nTotalAnt, dtype='complex64')
#crosscalpar[crosspols.keys()[0][1]][c.subsetant] = np.exp(1.j * c.antloc[:,:2].dot(median_sol))

for pol in ['xx','yy']:
    crosscalpar[pol[0]] = np.ones((rawdata.shape[1], rawdata.shape[2], c.nTotalAnt), dtype='complex64')
crosscalpar[crosspols.keys()[0][1]][..., c.subsetant] = np.exp(1.j * c.antloc[:,:2].dot(sol.transpose(1,0,2))).transpose(1,2,0)

cdata = np.empty_like(rawdata)
for p, pol in enumerate(crosspols.keys()):
    cdata[p] = omni.apply_calpar2(rawdata[p], crosscalpar[pol[0]], crosscalpar[pol[1]], c.totalVisibilityId)

#check application
#A = []
#bindex = []
#ublcount_sort = np.argsort(c.ublcount)
#for u in range(c.nUBL):
    #if c.ublcount[u] >= c.ublcount[ublcount_sort[-min(c.nUBL,99999)]]:
        #cbl = c.ublindex[u][:,2].astype(int)
        #mask = c.reversed[cbl] == 1
        #A = A + (c.antloc[c.ublindex[u][mask, 0].astype(int)] - c.antloc[int(c.ublindex[u][mask, 0][0])]).tolist()
        #bindex = bindex + [[a, int(c.ublindex[u][mask, 2][0])] for a in c.ublindex[u][mask, 2].astype(int)]#inside crossindex
#A = np.array(A)[:, :2]
#A = np.concatenate((A, -A), axis=0)
#bindex = c.subsetbl[c.crossindex[np.array(bindex).astype(int)]]#inside all bl

#b = np.empty((len(A), rawdata.shape[1], rawdata.shape[2]), dtype='float32')
#for p, pol in enumerate(crosspols.keys()):
    #b[p*len(bindex):(p+1)*len(bindex)] = (np.angle(cdata[p][..., bindex[:,0]]) - np.angle(cdata[p][..., bindex[:,1]])).transpose((2,0,1))

#crosstimer = omni.Timer()
#print "Solving cross calibration again....",
#sys.stdout.flush()
#sol2 = omni.solve_slope(A, b, 1, verbose=True)
#sol2[:,flags['xx']] = np.nan
#sol2[:,flags['yy']] = np.nan
#median_sol2 = np.array([np.nanmedian(sol2[i]) for i in range(len(sol2))])
#print "Done."
#sys.stdout.flush()

#################################################################################
############################solve for overall constant###########################
#################################################################################
crossextract, cross_chisq = omni.extract_crosspol_ubl(cdata, c.Info.get_info())
red_enough = c.ublcount > (max(c.ublcount) / 2) #treat all ubls above half max ublcount equally when computing the x-y angle
xy_candidates = np.empty(list(cdata.shape)[1:3] + [np.sum(red_enough)], dtype='float32')

for i,u in enumerate(np.arange(c.nUBL)[red_enough].astype(int)):
    xy_candidates[..., i] = np.angle(crossextract[0,...,u]/crossextract[1,...,u])

overall_angle = omni.medianAngle(omni.medianAngle(xy_candidates, axis = -1), axis = 0) / 2
for i in range(1, len(overall_angle)):
    offset = overall_angle[i-1]-PI/2
    overall_angle[i] = (overall_angle[i] + offset)%(PI) - offset
if make_plots:
    plt.plot(overall_angle)
    plt.ylim(-PI, PI)
    plt.show()

crosscalpar2={}
for pol in ['xx','yy']:
    crosscalpar2[pol[0]] = np.ones((rawdata.shape[1], rawdata.shape[2], c.nTotalAnt), dtype='complex64')
crosscalpar2[crosspols.keys()[0][1]][..., c.subsetant] = np.exp(1.j * overall_angle[None, :, None])

cdata2 = np.empty_like(rawdata)
for p, pol in enumerate(crosspols.keys()):
    cdata2[p] = omni.apply_calpar2(cdata[p], crosscalpar2[pol[0]], crosscalpar2[pol[1]], c.totalVisibilityId)

####apply cross pol solution to 'xx'############
apply_pol = crosspols.keys()[0][1] + crosspols.keys()[0][1]
#model sol as linear over freq
solf = nanmedian(sol, axis=-2)
Af = np.ones((len(solf[0]), 2), dtype='float32')
Af[:,0] = np.arange(len(solf[0]))
nan_freq = np.isnan(solf[0])|np.isnan(solf[1])
lin_solf = Af.dot(la.inv(Af[~nan_freq].transpose().dot(Af[~nan_freq])).dot(Af[~nan_freq].transpose().dot(solf.transpose()[~nan_freq]))).transpose()
##apply
linearcalpar[apply_pol] =  linearcalpar[apply_pol] * (np.exp(1.j * c.antloc[:,:2].dot(solf)).transpose()) * np.exp(1.j * overall_angle[:, None])
linearcalpar_model[apply_pol] = linearcalpar_model[apply_pol] * (np.exp(1.j * c.antloc[:,:2].dot(lin_solf)).transpose()) * np.exp(1.j * overall_angle[:, None])
crosstimer.tick("cross calibration")
####################################################################################
######################output results ##################################
####################################################################################
if oppath != "DONT_WRITE/":
    if info_tag == "DEFAULT":
        op_info_path = oppath + 'redundantinfo_first_cal_' + time.strftime("%Y_%m_%d_%H_%M_%S") + ".bin"
        op_calpar_path = oppath + 'calpar_first_cal_' + time.strftime("%Y_%m_%d_%H_%M_%S") + ".p"
        op_calpar_path_model = oppath + 'calpar_first_cal_' + time.strftime("%Y_%m_%d_%H_%M_%S") + ".pmodel"
    else:
        op_info_path = oppath + 'redundantinfo_first_cal_' + info_tag + ".bin"
        op_calpar_path = oppath + 'calpar_first_cal_' + info_tag + ".p"
        op_calpar_path_model = oppath + 'calpar_first_cal_' + info_tag + ".pmodel"

    print FILENAME + " MSG: Writing redundant info to %s"%op_info_path,
    sys.stdout.flush()
    calibrators['xx'].write_redundantinfo(infoPath = op_info_path, verbose = False, overwrite = overwrite)
    print "Done."
    sys.stdout.flush()

    #for pol in crosspols.keys():
        #print FILENAME + " MSG: Writing %s raw calpar to %s"%(pol, op_calpar_path.replace('.bin', pol+'.bin')),
        #sys.stdout.flush()
        #linearcalpar[pol].astype('float32').tofile(op_calpar_path.replace('.bin', pol+'.bin'))
        #print "Done."
        #sys.stdout.flush()


    for linpar, linpath in zip([linearcalpar, linearcalpar_model], [op_calpar_path, op_calpar_path_model]):
        linearcalpar_out = {}#include all antennas, not just good ones
        for pol in ['xx', 'yy']:
            linearcalpar_out[pol] = np.ones((linpar[pol].shape[0], calibrators[pol].nTotalAnt), dtype='complex64')
            linearcalpar_out[pol][:, calibrators[pol].subsetant] = linpar[pol]
        with open(linpath, 'wb') as outfile:
            print FILENAME + " MSG: Writing crude calpar to %s"%linpath,
            sys.stdout.flush()
            pickle.dump(linearcalpar_out, outfile, protocol=pickle.HIGHEST_PROTOCOL)
    print "Done."
    sys.stdout.flush()

else:
    print FILENAME + " MSG: Not outputting redundantinfo or rawcalpar by default."
    sys.stdout.flush()





print "%i Bad antennas found: "%len(badAntenna),
bad_str = ""
for bant in badAntenna:
    bad_str += (str(int(bant)) + ",")
if bad_str[-1] == ",":
    bad_str = bad_str[:-1]
print bad_str

for pol in n_bad_delay.keys():
    print "Possible bad antennas found from %s delay fitting: "%pol,
    bad_str = ""
    for bant in worst_delay_ant[pol]:
        bad_str += (str(int(bant)) + ",")
    if bad_str[-1] == ",":
        bad_str = bad_str[:-1]
    print bad_str
