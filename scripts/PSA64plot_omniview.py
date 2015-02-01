import numpy as np
import omnical.calibration_omni as omni
import matplotlib.pyplot as plt
import aipy as ap


paper_uvfiles = ['/data4/paper/2012EoR/psa_live/psa6266/zen.2456266.64714.uvcRREcAC']
omni_uvfiles = ["/data4/paper/2012EoR/psa_live/forlstbinning_omnical_2/zen.2456266.64714.uvcRREcACO"]
info=omni.read_redundantinfo('/data2/home/hz2ug/omnical/doc/redundantinfo_PSA64_ba19_37_50.bin')
nant = 64
wantpols = {}
wantpols['xx'] = ap.miriad.str2pol['xx']

papercal_data, _, _, _, rawflag = omni.importuvs(paper_uvfiles, wantpols, timingTolerance=100)
omnical_data, _, _, _, rawflag = omni.importuvs(omni_uvfiles, wantpols, timingTolerance=100)

omni.omniview(np.array([papercal_data[0,0,120], omnical_data[0,0,120]]), info)
