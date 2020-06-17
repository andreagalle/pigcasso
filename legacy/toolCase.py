# -*- coding: utf-8 -*-
"""
toolCase.py
"""

import draft as fct
import numpy as np
import string
import time
import os, sys, shutil
import matplotlib.pyplot as plt
from xml.dom.minidom import parse
import xml.dom.minidom
import scipy.signal
from scipy import stats
import scipy.integrate




def fct_read_file(file):
    rfile = open(file.strip(),'r')
    LINELIST = rfile.read().split('\n') # read all lines splitting them each time there is \n 
    # remove empty lines or lines with blanks only
    RMV_IDX = []
    for cnt in range(len(LINELIST)):
        if len(string.split(LINELIST[cnt].strip())) ==  0: RMV_IDX.append(cnt)
    for cnt in range(len(RMV_IDX)): del LINELIST[RMV_IDX[-cnt-1]]
    rfile.close()
    
    return LINELIST

def read_probes_pos(workdir,probe):
	    filex = workdir+probe+'/line_probe.X_CV'
	    filey = workdir+probe+'/line_probe.Y_CV'
	    filez = workdir+probe+'/line_probe.Z_CV'
	    X_POS = fct_read_file(filex)
	    Y_POS = fct_read_file(filey)
	    Z_POS = fct_read_file(filez)
	    
	    TMPX = X_POS[0].strip().split()
	    TMPY = Y_POS[0].strip().split()
	    TMPZ = Z_POS[0].strip().split()
	    nb_probes = int(TMPX[2])
	    data = [(TMPX[i],TMPY[i],TMPZ[i]) for i in range(3,nb_probes+3)]
	    data = np.array(data,dtype=[('x',float),('y',float),('z',float)])
	    
	    return data,nb_probes
    
def read_probes_data(workdir,probe,var,start=0):
    fileData = workdir+probe+'/line_probe.%s'%var
    
    LINELIST = fct_read_file(fileData)
    
    Data = [] ; Time = [] ; cnt = 0
    for line in LINELIST:
        TMP = line.strip().split()
        
        nb = int(TMP[2])
        
        if float(TMP[0])>=start:
			# if cnt%10==0:
            Time.append(TMP[1])
            Data.append(TMP[3:nb+3])
        
            cnt += 1
        
    Data = np.array(Data,dtype=float)
    Time = np.array(Time,dtype=float)
    return Data,Time

	# def calcAvgBulk(var)
	# Pinf_VRE = np.sum(P_VRE, axis=[1,2])/np.size(P_VRE)

def getVarDataAtWall(data_outVTK, ori, norm , var, alt=None):
    # function display 
    print '---- DAEPy::getVarDataAtPosition ----'

    # test if the field var is present
    if data_outVTK.GetPointData().HasArray(var)!=1:
        raise ValueError("Error : field %s not present"%var)

    # slice the 2D field data using "ori" as the origin and
    # "norm" as the normal vector 
    data1D = fct.getSlice(data_outVTK, ori, norm)

    # extract the tangential velocity
    [Vcoords, data] = fct.getArrayFromPointData(data1D, [var])

	# print Vcoords

    # define a new scalar coordinates along the line orthogonal to the wall
	# x_coo = np.array(fct.getScalarCoord(Vcoords, 1))
    x_coo = np.array(Vcoords[:,0])
    id_sort = x_coo.argsort()

    # sort cordinates along the line orthogonal to the wall
    x_coo = x_coo[id_sort] - x_coo[id_sort[0]]
    data = data[id_sort]

    if alt == None:
        return [x_coo, data]
    else:
        return [ori, data[alt]]

# this function calculates the span-wise and stream-wise spectra
def computeSpectrum(nx,nz,grid_x,grid_z,grid_data):

    KZ = [] ; KX = []
    FFTZ = [] ; FFTX = []

    for i in xrange(nx):
        signal = scipy.signal.detrend(grid_data.T[i], type='constant') # the mean of data is subtracted (detrending).
        fftz = np.fft.rfft(signal, norm='ortho') 
        n = signal.size # number of array's elements = product of the array’s dimensions.
        dz = 3./63 # np.ediff1d(grid_z.T[i]) # spacing between the elemtns of the i-th column of grid_z.
        fftz = fftz*np.conj(fftz) # to calculate the module of the complex fftz? is what we want to plot!
        kz = 2.*np.pi*np.fft.rfftfreq(n, d=dz)#[0]) # d: Sample spacing (inv. sampling rate: number samples per second), defaults=1. 

        KZ.append(kz)
        FFTZ.append(fftz)

    for i in xrange(nz):
        signal = scipy.signal.detrend(grid_data[i], type='constant')
        fftx = np.fft.rfft(signal, norm='ortho')
        n = signal.size
        dx = 6./63 # np.ediff1d(grid_x[i]) MODIFIED TO AVOID SOME ISSUES DUE TO POINT (AND NOT CELL) DATAS...
        fftx = fftx*np.conj(fftx)
        kx = 2.*np.pi*np.fft.rfftfreq(n, d=dx)#[0]) BEEING NOT AN ARRAY, BUT A SCALAR VALUE!

        KX.append(kx)
        FFTX.append(fftx)

    return [KX, FFTX, KZ, FFTZ]

def readDNSdb_Umean(workdir,database,var):
	fileData = workdir + '/DNS_profiles/' + database
	LINELIST = fct_read_file(fileData)

	coord = [] ; coordP = [] ; profile = [] ; RMV_IDX = []

	for cnt in range(len(LINELIST)):
		if LINELIST[cnt].startswith("%") | LINELIST[cnt].startswith("--") | LINELIST[cnt].startswith("#") is True: RMV_IDX.append(cnt) 
	for cnt in range(len(RMV_IDX)): del LINELIST[RMV_IDX[-cnt-1]]

	for line in LINELIST: 
		TMP = line.strip().split()
		# if (var == 'u+' | var == 'up' | var == 'Up'): # is True:
		coord.append(TMP[0])
		coordP.append(TMP[1])
		profile.append(TMP[2])
		# else: raise ValueError("Error : Input variable not known")
	coord = np.array(coord,dtype=float)
	coordP = np.array(coordP,dtype=float)
	profile = np.array(profile,dtype=float)
	return coord,coordP,profile

def readDNSdb_ReStress(workdir,database):
	fileData = workdir + '/DNS_profiles/' + database
	LINELIST = fct_read_file(fileData)

	profileUU = [] ; profileVV = [] ; profileWW = [] ; profileUV = [] ; RMV_IDX = []

	for cnt in range(len(LINELIST)):
		if LINELIST[cnt].startswith("%") | LINELIST[cnt].startswith("--") | LINELIST[cnt].startswith("#") is True: RMV_IDX.append(cnt) 
	for cnt in range(len(RMV_IDX)): del LINELIST[RMV_IDX[-cnt-1]]

	for line in LINELIST: 
		TMP = line.strip().split()
		# if (var == 'u+' | var == 'up' | var == 'Up'): # is True:
		# coord.append(TMP[0])
		# coordP.append(TMP[1])
		profileUU.append(TMP[2])
		profileVV.append(TMP[3])
		profileWW.append(TMP[4])
		profileUV.append(TMP[5])
		# else: raise ValueError("Error : Input variable not known")
	# coordP = np.array(coordP,dtype=float)
	profileUU = np.array(profileUU,dtype=float)
	profileVV = np.array(profileVV,dtype=float)
	profileWW = np.array(profileWW,dtype=float)
	profileUV = np.array(profileUV,dtype=float)
	return profileUU,profileVV,profileWW,profileUV



# this function calculate the mean velocity slope near the wall, using a linear regression 

# def calc_tau_wall(u,y,Muinf):
# 	slope, intercept, r_value, p_value, std_err = stats.linregress(y[:10],u[:10])
# 	tau_wall = Muinf * slope
# 	print 'JUST TO CHECK: the interpolated du/dy|0 value is:', slope
# 	return tau_wall



# this function should read the mean profiles calculated, but is not completed

# def readChanDatas(workdir,model)
# 	fileData = resdir + '/Profiles_%s.dat'%model
# 	LINELIST = fct_read_file(fileData)

# 	yp = [] ; Up = [] ; Mu_sgsp =[] ; y_coo = [] ; Ut = [] ; Urms = [] ; Urey = [] ;  RMV_IDX = []

# 	for cnt in range(len(LINELIST)):
# 		if LINELIST[cnt].startswith("%") | LINELIST[cnt].startswith("--") | LINELIST[cnt].startswith("#") is True: RMV_IDX.append(cnt) 
# 	for cnt in range(len(RMV_IDX)): del LINELIST[RMV_IDX[-cnt-1]]

# 	for line in LINELIST: 
# 		TMP = line.strip().split()
# 		yp.append(TMP[1])
# 		Up.append(TMP[2])
# 		Mu_sgsp.append(TMP[2])
# 		y_coo.append(TMP[2])
# 		Ut.append(TMP[2])
# 		Urms.append(TMP[2])
# 		Urey.append(TMP[2])

# 	coord = np.array(coord,dtype=float)
# 	profile = np.array(profile,dtype=float)
# 	return coord,profile


# this function read the spectra once they have been calculated from computeSpectrum, the first time, 
# and the results have been saved in the result directory 
def readSpectraDatas(workdir,model):

	coordX = [] ; coordZ = []  
	spectraX = [] ; spectraZ = []

	ID = ['1', '2', '3', '4']
	for plane in range(4):
		idp = ID[plane]

		modZ1 = [] ; modZ2 = [] ; modZ3 = [] ; specZ1 = [] ; specZ2 = [] ; specZ3 = []
		SpanData = workdir + '/PP_results/spectra_plot/plane_%s'%idp + '/SpanSpectrum_plane%s_%s.dat'%(idp,model)
		LINELIST = fct_read_file(SpanData)

		RMV_IDX = []
		for cnt in range(len(LINELIST)):
			if LINELIST[cnt].startswith("%") | LINELIST[cnt].startswith("--") | LINELIST[cnt].startswith("#") is True: RMV_IDX.append(cnt) 
		for cnt in range(len(RMV_IDX)): del LINELIST[RMV_IDX[-cnt-1]]

		for line in LINELIST: 
			TMP = line.strip().split()
			modZ1.append(TMP[0]) 
			modZ2.append(TMP[0]) 
			modZ3.append(TMP[0])
			specZ1.append(TMP[1])
			specZ2.append(TMP[2])
			specZ3.append(TMP[3])

		modZ = [np.array(modZ1, dtype=float), np.array(modZ2, dtype=float), np.array(modZ3, dtype=float)]
		modZ = np.array(modZ).real 
		specZ = [np.array(specZ1, dtype=float), np.array(specZ2, dtype=float), np.array(specZ3, dtype=float)]
		specZ = np.array(specZ).real 
		# print modZ
		# print 'shape of modZ e specZ : ', modZ.shape, specZ.shape

		#########################################################

		modX1 = [] ; modX2 = [] ; modX3 = [] ; specX1 = [] ; specX2 = [] ; specX3 = []
		StreamData = workdir + '/PP_results/spectra_plot/plane_%s'%idp + '/StreamSpectrum_plane%s_%s.dat'%(idp,model)
		LINELIST = fct_read_file(StreamData)

		RMV_IDX = []
		for cnt in range(len(LINELIST)):
			if LINELIST[cnt].startswith("%") | LINELIST[cnt].startswith("--") | LINELIST[cnt].startswith("#") is True: RMV_IDX.append(cnt) 
		for cnt in range(len(RMV_IDX)): del LINELIST[RMV_IDX[-cnt-1]]

		for line in LINELIST: 
			TMP = line.strip().split()
			modX1.append(TMP[0]) 
			modX2.append(TMP[0]) 
			modX3.append(TMP[0])
			specX1.append(TMP[1])
			specX2.append(TMP[2])
			specX3.append(TMP[3])

		modX = [np.array(modX1, dtype=float), np.array(modX2, dtype=float), np.array(modX3, dtype=float)]
		modX = np.array(modX).real 
		specX = [np.array(specX1, dtype=float), np.array(specX2, dtype=float), np.array(specX3, dtype=float)]
		specX = np.array(specX).real 
		# print 'shape of modX e specX : ', modX.shape, specX.shape

		#########################################################

		coordX.append(modX) ; coordZ.append(modZ)
		spectraX.append(specX) ; spectraZ.append(specZ)

	return [coordX, spectraX, coordZ, spectraZ]

# this function is needed to calculate the grid dimension in the Z direction
def getVarDataAtWall_GridZ(data_outVTK, ori, norm , var, alt=None):
	# function display 
	print '---- DAEPy::getVarDataAtPosition ----'

	# test if the field var is present
	if data_outVTK.GetPointData().HasArray(var)!=1:
		raise ValueError("Error : field %s not present"%var)

	# slice the 2D field data using "ori" as the origin and
	# "norm" as the normal vector 
	data1D = fct.getSlice(data_outVTK, ori, norm)
	
	# extract the tangential velocity
	[Vcoords, data] = fct.getArrayFromPointData(data1D, [var])

	# print Vcoords

	# define a new scalar coordinates along the line orthogonal to the wall
	# z_coo = np.array(fct.getScalarCoord(Vcoords, 1))
	z_coo = np.array(Vcoords[:,2])
	id_sort = z_coo.argsort()

	# sort cordinates along the line orthogonal to the wall
	z_coo = z_coo[id_sort] - z_coo[id_sort[0]]
	data = data[id_sort]

	if alt == None:
		return [z_coo, data]
	else:
		return [ori, data[alt]]

def readDNSdb_Spectra(workdir,database):
	kDNS = [] ; SpectraDNS = []
	planes = [5, 10, 98, 392]
	for idp in planes:
		k = [] ; E_uu = [] ; E_vv = [] ; E_ww = []
		fileData = workdir + '/DNS_profiles/spectra' + database + '%d'%idp
		LINELIST = fct_read_file(fileData)

		RMV_IDX = []
		for cnt in range(len(LINELIST)):
			if LINELIST[cnt].startswith("%") | LINELIST[cnt].startswith("--") | LINELIST[cnt].startswith("#") is True: RMV_IDX.append(cnt) 
		for cnt in range(len(RMV_IDX)): del LINELIST[RMV_IDX[-cnt-1]]

		for line in LINELIST: 
			TMP = line.strip().split()
			k.append(TMP[0])
			E_uu.append(TMP[1])
			E_vv.append(TMP[2])
			E_ww.append(TMP[3])
			# else: raise ValueError("Error : Input variable not known")
			# coordP = np.array(coordP,dtype=float)

		k = np.array(k,dtype=float)
		E_uu = np.array(E_uu,dtype=float)
		E_vv = np.array(E_vv,dtype=float)
		E_ww = np.array(E_ww,dtype=float)
		E = np.array([E_uu, E_vv, E_ww],dtype=float)

		kDNS.append(k)
		SpectraDNS.append(E)

	return kDNS, SpectraDNS

def readTimePerIter(workdir):
	models = ['VRE', 'DSM', 'WAL', 'NON']
	minTimeXiter_MOD = [] ; meanTimeXiter_MOD = [] ; StDtimeXiter_MOD = [] ; timeXiter_MOD = []
	for mod in models:
		fileData = workdir + '/%s_cc/'%mod + '/timeperiter.dat'

		LINELIST = fct_read_file(fileData)

		timeXiter = [] ; # RMV_IDX = []

		# for cnt in range(len(LINELIST)):
		# 	# if LINELIST[cnt].startswith("%") is True: RMV_IDX.append(cnt)
		# 	# if LINELIST[cnt]>=50 is True: RMV_IDX.append(cnt)
		# 	if LINELIST[cnt]>=50 is True: LINELIST[cnt]= 0
		# # for cnt in range(len(RMV_IDX)): del LINELIST[RMV_IDX[-cnt-1]]

		for line in LINELIST: 
			TMP = line.strip().split()
			# if TMP[0]>50: 
			# 	TMP[0]=50
			timeXiter.append(TMP[0])
		timeXiter = np.array(timeXiter,dtype=float)
		# for cnt in np.arange(len(timeXiter)):
		# 	if timeXiter[cnt]>50: 
		# 		timeXiter[cnt]=0
		minTime = min(timeXiter)
		minTimeXiter_MOD.append(minTime)

		meanTime = np.mean(timeXiter)
		meanTimeXiter_MOD.append(meanTime)

		StDtimeXiter = np.std(timeXiter)
		StDtimeXiter_MOD.append(StDtimeXiter)

		timeXiter_MOD.append(timeXiter)


	return minTimeXiter_MOD, meanTimeXiter_MOD, StDtimeXiter_MOD, timeXiter_MOD

















