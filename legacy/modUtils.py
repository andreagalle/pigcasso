# -*- coding: utf-8 -*-
"""
@author: a.grebert
"""

import draft as fct
import numpy as np
import string
import time
import os, sys, shutil
import matplotlib.pyplot as plt
from xml.dom.minidom import parse
import xml.dom.minidom

def time_stat(tstart):
    minutes, seconds= divmod(time.time()-tstart, 60)

    if minutes != 0:
        txt_min = '%3d minute(s)'%minutes
    else:
        txt_min = ''
        
    return txt_min,seconds

def fct_progress_bar(cmd,iter,width):
    bar_strg = '|%s|' %('.'*width)
    fillsign = (u'\N{BLACK SQUARE}')
    bar_strg = bar_strg.replace('.',fillsign,iter)
    print '\r%s %s' %(cmd,bar_strg) ,
    sys.stdout.flush()

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def fct_dir_exist(dir):
    check = True
    if os.path.isdir(dir) == False:
        check = False
    return check

def fct_rawinp(strg):
    out = raw_input(strg)
    return out

def fct_query(question,selection):
    # Function to get user input (selection of string combinations like 'y','n','2'...)
    query = '---'
    while not query in selection:
        query = string.lower(fct_rawinp('>>>>>>>> '+question))
        if len(query) == 0: query = '---'
        elif query == 'x': sys.exit()
    return query

def fct_read_file(file):
    rfile = open(file.strip(),'r')
    LINELIST = rfile.read().split('\n')                 # read all lines
    # remove empty lines or lines with blanks only
    RMV_IDX = []
    for cnt in range(len(LINELIST)):
        if len(string.split(LINELIST[cnt].strip())) ==  0: RMV_IDX.append(cnt)
    for cnt in range(len(RMV_IDX)): del LINELIST[RMV_IDX[-cnt-1]]
    rfile.close()
    
    return LINELIST
    
def fct_read_exp(file):
    
    LINELIST = fct_read_file(file)
    
    data = []
    for line in LINELIST:
        TMP = line.strip().split()
        if TMP[0][0] <> '%' and TMP[0][0] <> '#':
            TMP = [float(i) for i in TMP]
            # l = np.array((TMP[0],TMP[1],TMP[2],TMP[3],TMP[4],TMP[5],TMP[6]))
            l = np.array((TMP))
            data.append(l)
        
    data = np.asarray(data).T    
    
    return data

def fct_read(file):
    
    LINELIST = fct_read_file(file)
    
    data = [] ; header = []
    for line in LINELIST:
        TMP = line.strip().split()
        if TMP[0][0] == '%' or TMP[0][0] == '#':
            header.append(TMP[1:])
        # if TMP[0][0] <> '%' and TMP[0][0] <> '#':
        else:
            TMP = [float(i) for i in TMP]
            # l = np.array((TMP[0],TMP[1],TMP[2],TMP[3],TMP[4],TMP[5],TMP[6]))
            l = np.array((TMP))
            data.append(l)
        
    data = np.asarray(data).T
    header = np.asarray(header)
    
    return [header, data]

def fct_read_Cf(file):
    
    LINELIST = fct_read_file(file)
    
    data = []
    for line in LINELIST:
        TMP = line.strip().split()
        if TMP[0] <> '#' and TMP[0] <> '%':
            TMP = [float(i) for i in TMP]
            l = np.array((TMP[0],TMP[1]))
            data.append(l)
        
    data = np.asarray(data).T    
    
    return data

def fct_read_Schaltter_LES(workdir):
    
    data_re = fct_read_exp(workdir + '/re.data')
    data_cf = fct_read_exp(workdir + '/cf.data')

    data = np.array((data_re[2,:5000],data_cf[1,:5000])) 
    
    return data

def read_probes_pos(workdir,probe):
    filex = workdir+probe+'_line.X_CV'
    filey = workdir+probe+'_line.Y_CV'
    filez = workdir+probe+'_line.Z_CV'
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
    fileData = workdir+probe+'_line.%s'%var
    
    LINELIST = fct_read_file(fileData)
    
    Data = [] ; Time = [] ; cnt = 0
    for line in LINELIST:
        TMP = line.strip().split()
        
        nb = int(TMP[2])
        
        if float(TMP[0])>=start:
#            if cnt%10==0:
            Time.append(TMP[1])
            Data.append(TMP[3:nb+3])
        
            cnt += 1
        
    Data = np.array(Data,dtype=float)
    Time = np.array(Time,dtype=float)
    return Data,Time

def getInputs(workdir,inp='/inputPy.in',inf=True):

    field = ['data', 'wall', 'norm', 'ori', 'length', 'x_pos', 'delta_i',\
              'x_start', 'Pinf', 'Tinf', 'Rhoinf', 'Uinf', 'Muref', 'Mupower', 'Tref']

    inFile = workdir + '/' + inp
    if os.path.isfile(inFile) == False: 
      raise ValueError("Error : Input file does not exist")
    linelist =  fct_read_file(inFile)

    data = []
    for line in linelist:
        TMP = line.strip().split()
        if TMP[0] <> '#':
            data.append(TMP)

    data = np.asarray(data)

    file_d = os.path.join(workdir,data[np.where(data==field[0])[0][0]][-1])
    file_w = os.path.join(workdir,data[np.where(data==field[1])[0][0]][-1])
    nor1 = np.array(data[np.where(data==field[2])[0][0]][-1].strip("[").strip("]").split(","),dtype=float)
    ori1 = np.array(data[np.where(data==field[3])[0][0]][-1].strip("[").strip("]").split(","),dtype=float)
    length = np.float(data[np.where(data==field[4])[0][0]][-1])
    x_pos =  np.float(data[np.where(data==field[5])[0][0]][-1])
    delta_i = np.float(data[np.where(data==field[6])[0][0]][-1])
    x_start = np.float(data[np.where(data==field[7])[0][0]][-1])
    if inf == True:
        Pinf = np.float(data[np.where(data==field[8])[0][0]][-1])
        Tinf = np.float(data[np.where(data==field[9])[0][0]][-1])
        Rhoinf = np.float(data[np.where(data==field[10])[0][0]][-1])
        Uinf = np.float(data[np.where(data==field[11])[0][0]][-1])
        Muref = np.float(data[np.where(data==field[12])[0][0]][-1])
        Mupower = np.float(data[np.where(data==field[13])[0][0]][-1])
        Tref = np.float(data[np.where(data==field[14])[0][0]][-1])
        Muinf = Muref*np.power(Tinf/Tref,Mupower)
    elif inf == False:
        Pinf = None
        Tinf = None
        Rhoinf = None
        Uinf = None
        Muinf = None

    x_pos = x_pos/length

    return [file_d, file_w, nor1, ori1, length, x_pos, delta_i, x_start, Pinf, Tinf, Rhoinf, Uinf, Muinf]

def getVarDataAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, var, alt=None):
    # function display 
    print '---- DAEPy::getVarDataAtPosition ----'

    # test if the field var is present
    if data_outVTK.GetPointData().HasArray(var)!=1:
        raise ValueError("Error : field %s not present"%var)

    # get the vector used to slice the field data
    [ori, vec_tan] = fct.getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    data1D = fct.getSlice(data_outVTK, ori, vec_tan)

    # extract the tangential velocity
    [Vcoords, data] = fct.getArrayFromPointData(data1D, [var])

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(fct.getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()

    # sort cordinates along the line orthogonal to the wall
    y_coo = y_coo[id_sort]
    data = data[id_sort]

    if alt == None:
        return [y_coo, data]
    else:
        return [ori, data[alt]]
        
def getVarDataAtPositionNoWall(data_outVTK, ori, norm , var, alt=None):
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

#    print Vcoords

    # define a new scalar coordinates along the line orthogonal to the wall
#    y_coo = np.array(fct.getScalarCoord(Vcoords, 1))
    y_coo = np.array(Vcoords[:,1])
    id_sort = y_coo.argsort()

    # sort cordinates along the line orthogonal to the wall
    y_coo = y_coo[id_sort] - y_coo[id_sort[0]]
    data = data[id_sort]

    if alt == None:
        return [y_coo, data]
    else:
        return [ori, data[alt]]

def getVarDataBetweenPosition(data_outVTK, wall_outVTK, cord_choice, x_p0, x_p1, Npts, var, alt):
    """
    
    """

    # function display 
    print '---- DAEPy::getVarDataBetweenPosition ----'

    data = np.zeros(Npts)
    pos = np.zeros((Npts, 3))

    for i in range(Npts):
        
        x_p_temp = x_p0 + (x_p1-x_p0)/(Npts-1)*i

        [_, data[i]] = getVarDataAtPosition(data_outVTK, wall_outVTK, cord_choice, x_p_temp, var, alt)

    return [_, data]

def readCBC(cbcdir, inputXMLFileName):

    # Open XML document using minidom parser
    dom = xml.dom.minidom.parse(cbcdir+inputXMLFileName)

    #       1) Type of profile
    profilesType = dom.getElementsByTagName("Profiles")[0]. \
      attributes["Type"].value

    #       2) Wall-normal coordinates
    coordinatesX2 = \
      dom.getElementsByTagName("Coordinates_2")[0].firstChild.data
    scalingX2 = \
      dom.getElementsByTagName("Coordinates_2")[0]. \
      attributes["ScalingFactor"].value


    coordinatesX2 = coordinatesX2.strip().replace("\n"," ").split(" ")
    coordinatesX2 = np.array([float(coordinatesX2[idx]) for idx in range(len(coordinatesX2))])
    coordinatesX2 = coordinatesX2*float(scalingX2)

    nbPoints = len(coordinatesX2)

    #       3) Spanwise coordinates
    coordinatesX3 = \
      dom.getElementsByTagName("Coordinates_3")[0].firstChild.data
    scalingX3 = \
      dom.getElementsByTagName("Coordinates_3")[0]. \
      attributes["ScalingFactor"].value


    coordinatesX3 = coordinatesX3.strip().replace("\n"," ").split(" ")
    coordinatesX3 = np.array([float(coordinatesX3[idx]) for idx in range(len(coordinatesX3))])
    coordinatesX3 = coordinatesX3*float(scalingX3)

    nbProfiles = len(coordinatesX3)


    #       4) Fields
    fields = dom.getElementsByTagName("Field")
    data = np.zeros((len(fields),len(coordinatesX2)))
    dataName = []

    for i,field in enumerate(fields):
        fieldName          = field.attributes["Name"].value
        fieldScalingFactor = field.attributes["ScalingFactor"].value
        # fieldValues        = dom.getElementsByTagName(fieldName)[0].firstChild.data

        fieldValues = field.firstChild.nodeValue.strip().split(" ")
        fieldValues = np.array([float(fieldValues[idx]) for idx in range(len(fieldValues)) if fieldValues[idx]]).\
                      reshape((nbProfiles,nbPoints))
        data[i] = fieldValues[0]*float(fieldScalingFactor)
        dataName.append(fieldName)

    return [coordinatesX2, data, dataName]

def getPositionExp(workdir,target,data='theta',comp=False):

    if comp == True:
        LINELIST = fct_read_file(workdir + '/Reynolds.dat')
        theta = 'ReTheta'
    elif comp == False:
        LINELIST = fct_read_file(workdir + '/Cfi.dat')
        theta = 'ReTheta_i'

    ReTheta = [] ; ReTau = [] ; x_coo = []
    for line in LINELIST:
        TMP = line.strip().split()
        if TMP[0] == '#':
            idx = TMP.index('x') - 1 
            idTheta = TMP.index(theta) - 1
            idTau = TMP.index('ReTau') - 1
        if TMP[0] <> '#':
            TMP = [float(i) for i in TMP]
            x_coo.append(TMP[idx])
            ReTheta.append(TMP[idTheta])
            ReTau.append(TMP[idTau])
    
    if data == 'theta':
        Re = np.asarray(ReTheta)
    elif data == 'tau':
        Re = np.asarray(ReTau)

    if min(Re) > target or max(Re) < target:
        raise ValueError("Error : Experimental value of Reynolds theta is out of bound of dataset")

    ordre = abs(Re-target).argsort()[0:2]
    ordre = np.sort(ordre)

    # compute the position using a simple linear model
    xa = x_coo[ordre[0]]
    Rea = Re[ordre[0]]
    xb = x_coo[ordre[1]] 
    Reb = Re[ordre[1]]
    tau = (Reb-Rea)/(xb-xa)
    pos = xa + (target-Rea)/tau

    return pos

