import datetime
import socket, multiprocessing, math, random, traceback, ephem, string, commands, datetime
import time
from time import ctime
import aipy as ap
import struct
import numpy as np
import os, sys
import datetime
from optparse import OptionParser
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import scipy as sp
    import scipy.sparse as sps
    import scipy.linalg as la
    
import calibration_omni as omni

correctinfo = omni.read_redundantinfo('redundantinfo_PSA32.txt')
calibrator = omni.RedundantCalibrator(64)

calibrator.antennaLocationTolerance = .1
calibrator.badAntenna = []
calibrator.badUBL = [0, 1, 2]

#1234
#calibrator.antennaLocation = 

calibrator.compute_info()


