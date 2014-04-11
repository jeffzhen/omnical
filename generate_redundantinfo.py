import datetime
import socket, multiprocessing, math, random, traceback, ephem, string, commands, datetime
import time
from time import ctime
import struct
import numpy as np
import os, sys
import datetime
from optparse import OptionParser
import calibration_omni as omni
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import scipy as sp
    import scipy.sparse as sps
    import scipy.linalg as la

FILENAME = "*generate_redundantinfo.py*"

def match_vec(vec, item, tol):##vec is 2d array, this method finds the first index i at which norm(vec[i] - item) < tol
	for i in range(len(vec)):
		if sum((vec[i] - item)**2.)**0.5 < tol:
			return i, 1
		if sum((vec[i] + item)**2.)**0.5 < tol:
			return i, -1
	return -1, -1

def generate_redundantinfo(arrayinfopath, outputpath):
	with open(arrayinfopath) as f:
		rawinfo = [[float(x) for x in line.split()] for line in f]
	METHODNAME = "generate_redundantinfo"
	print FILENAME + METHODNAME + " MSG:",  "Reading", arrayinfopath, "...",
	badant = np.array(rawinfo[0]).astype(int)
	if badant[0] < 0:
		badant = []
	badubl = np.array(rawinfo[1]).astype(int)
	if badubl[0] < 0:
		badubl = []
	
	tol = rawinfo[2][0]
	
	cnter = 3
	antloc = []
	if len(rawinfo[cnter]) != 3:
		print FILENAME + METHODNAME + " MSG:",  "Format error in", arrayinfopath, ": The antenna locations should start on the third line, with 3 numbers in each line!"
		return
	while len(rawinfo[cnter]) == 3:
		antloc += [rawinfo[cnter]]
		cnter += 1
	nantT = len(antloc)
	antloc = np.array(antloc)

	bl2dT = []
	if len(rawinfo[cnter]) != 2:
		print FILENAME + METHODNAME + " MSG:",  "Format error in", arrayinfopath, ": The baseline to antenna mapping should start after antenna locations, with 2 numbers (conj index, index) in each line!"
		return
	while len(rawinfo[cnter]) == 2:
		bl2dT += [rawinfo[cnter]]
		cnter += 1
	bl2dT = np.array(bl2dT).astype(int)
		
	
	print badant
	print badubl

	
info0 = omni.read_redundantinfo('./redundantinfo_PSA32.txt')
generate_redundantinfo('./arrayinfo_PAPER32.txt', '')

