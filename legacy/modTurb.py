# -*- coding: utf-8 -*-
"""
Created on Fri Jan 7 08:08:12 2016

@author: a.grebert
"""

import vtk
import draft as fct
import numpy as np
import string
import os, sys, shutil
import matplotlib.pyplot as plt
import modUtils as util
import numpy.linalg as la

def get1DSpectrum(myArray, Dz):

    Ruu = np.correlate(myArray,myArray,mode='full')#[-myArray.shape[0]:]
    Ruu = Ruu[Ruu.size/2:]

    # fig = plt.figure(figsize=(10, 6))
    # ax = fig.add_subplot(111)
    # ax.plot(Ruu/Ruu[0])
    # # ax.set_xscale('log')
    # # ax.set_yscale('log')
    # plt.show()

    # Euu = np.zeros(((Nz-1)/2))
    # kz = np.zeros(((Nz-1)/2))
    # for n in range((Nz-1)/2):
    #     su = 0
    #     for k in range(1,(Nz-1)/2):
    #         su += Ruu[k]*np.cos(2.0*np.pi*n*k/(Nz-1))
    #     Euu[n] = 1 + 2*np.abs(su)
    #     kz[n] = n/Dz*1e-3

    # print Euu.shape, kz.shape
    # print Euu
    # print kz

    Euu = np.abs(np.fft.rfft(Ruu))**2
    kz = np.fft.rfftfreq(Ruu.size, Dz)
    idx = np.argsort(kz)

    Euu = Euu[idx] ; kz = kz[idx]

    return [kz, Euu, Ruu]

def getTurbulenceTriangle(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf=None, Rhoinf=None):

    """
    

    Parameters
    ----------

    wall_outVTK : VTK output object from a VTK reader
        This field is supposed to describe a curve with 1D cell type, like in getLineTangentialVector().

    cord_choice : integer, 0, 1 or 2
        Gives the axis that is going to be used to define the point where the velocity
        profile is extracted. Convention :
            - 0 = x axis
            - 1 = y axis
            - 2 = z axis

    x_perc : float between 0 and 1
        Gives the percentage of the relative cord along the axis defined with cord_d.

    Uinf : float 
        (optional) Free stream velocity

    Rhoinf : float 
        (optional) Free stream density

    Returns
    -------
    ori : tuple(3)
        The point coordinates of the point at the wall where we extract the data.

    yp : vector of floats
        dimensionless wall normal coordinates at x_perc.

    Urms : vector of floats (3 components)
        normalised (using Utau) Reynolds stress tensor components : Rxx, Ryy, Rzz

    Urey : vector of floats (3 components) 
        normalised (using Utau) Reynolds stress tensor components : Ryz, Rzx, Rxy

    scaling : vector of floats
        Density correction factor sqrt(Rho/Rho[0])

    """

    # function display 
    print '---- DAEPy::getRMSAtPosition ----'

    # test if the field RHO_AVG is present
    if data_outVTK.GetPointData().HasArray('RHO_AVG')!=1:
        raise ValueError("Error : field RHO_AVG not present")

    # test if the field MU_LAM_AVG is present
    if data_outVTK.GetPointData().HasArray('MU_LAM_AVG')!=1:
        raise ValueError("Error : field MU_LAM_AVG not present")

    # test if the field U_RMS is present
    if data_outVTK.GetPointData().HasArray('U_RMS')!=1:
        raise ValueError("Error : field U_RMS not present")

    # test if the field U_REY is present
    if data_outVTK.GetPointData().HasArray('U_REY')!=1:
        raise ValueError("Error : field U_REY not present")

    # get the vector used to slice the field data
    [ori, vec_tan] = fct.getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    data1D = fct.getSlice(data_outVTK, ori, vec_tan)

    # extract the tangential velocity
    [Vcoords, Urms, Urey, Rho, Mu] = fct.getArrayFromPointData(data1D, ['U_RMS','U_REY','RHO_AVG','MU_LAM_AVG'])

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(fct.getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()
    
    # sort cordinates along the line orthogonal to the wall
    y_coo = y_coo[id_sort]
    Rho = Rho[id_sort]
    Mu = Mu[id_sort]

    [_, Cf, Utau, _, _] = fct.getSkinFrictionAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf, Rhoinf)

    # sort the RMS and REY quantities according to the coordinates y_coo. Normalisation using Utau is also performed

    Urms = (Urms[id_sort]**2).T
    Urey = (Urey[id_sort]).T

    yp = y_coo*Utau*Rho[0]/Mu[0]

    IIb = np.zeros((len(yp)))
    IIIb = np.zeros((len(yp)))
    Ni = np.zeros((len(yp)))
    Ei = np.zeros((len(yp)))

    for idy in range(len(yp)):

        k = (Urms[0,idy]+Urms[1,idy]+Urms[2,idy])

        bij = np.matrix(([[Urms[0,idy]/k-1/3., Urey[2,idy]/k,Urey[1,idy]/k],\
                         [Urey[2,idy]/k, Urms[1,idy]/k-1/3., Urey[0,idy]/k],\
                         [Urey[1,idy]/k, Urey[0,idy]/k, Urms[2,idy]/k-1/3.]]))

        # ubji = np.triu(bij**2,0)
        # ubji = np.triu(np.multiply(bij,bij),0)

        # IIb[idy] = -0.5*np.sum(ubji)
        # IIIb[idy] = la.det(bij)
        IIb[idy] = -(bij[0,0]**2+bij[1,1]**2+bij[2,2]**2)/2.

        IIIb[idy] = (bij[0,0]**3+bij[1,1]**3+bij[2,2]**3)/3.

        Ni[idy] = np.sqrt(-1./3.*IIb[idy])
        Ei[idy] = np.power(np.abs(IIIb[idy])*0.5,1./3.)

    # function output
    print ''
    return [ori, y_coo, yp, Ei, Ni, IIIb, IIb]

def getStressAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf=None, Rhoinf=None):

    # function display 
    print '---- DAEPy::getRMSAtPosition ----'

    # test if the field RHO_AVG is present
    if data_outVTK.GetPointData().HasArray('RHO_AVG')!=1:
        raise ValueError("Error : field RHO_AVG not present")

    # test if the field MU_LAM_AVG is present
    if data_outVTK.GetPointData().HasArray('MU_LAM_AVG')!=1:
        raise ValueError("Error : field MU_LAM_AVG not present")

    # test if the field U_RMS is present
    if data_outVTK.GetPointData().HasArray('U_RMS')!=1:
        raise ValueError("Error : field U_RMS not present")

    # test if the field U_REY is present
    if data_outVTK.GetPointData().HasArray('U_REY')!=1:
        raise ValueError("Error : field U_REY not present")

    # test if the field TauWallAvg is present
    if wall_outVTK.GetPointData().HasArray('TauWallAvg')!=1:
        raise ValueError("Error : field TauWallAvg not present")

    # get the vector used to slice the field data
    [ori, vec_tan] = fct.getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    data1D = fct.getSlice(data_outVTK, ori, vec_tan)
    wall1D = fct.getSlice(wall_outVTK, ori, vec_tan)

    # extract the shear stress
    [_, Tw] = fct.getArrayFromPointData(wall1D, ['TauWallAvg'])

    # extract the tangential velocity
    [Vcoords, U, Urms, Urey, Rho, Mu, MuSGS] = fct.getArrayFromPointData(data1D, ['U_AVG','U_RMS','U_REY','RHO_AVG','MU_LAM_AVG','MU_SGS_AVG'])
    Ut = np.sum(U*vec_tan, axis=1)

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(fct.getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()
    
    # sort cordinates along the line orthogonal to the wall
    y_coo = y_coo[id_sort]
    Rho = Rho[id_sort]
    Mu = Mu[id_sort]
    MuSGS = MuSGS[id_sort]
    Ut = Ut[id_sort]
    U = U[id_sort]
    Urms = Urms[id_sort]**2
    Urey = Urey[id_sort]
    Utau = np.sqrt(np.sqrt(Tw[0,0]**2+Tw[0,1]**2+Tw[0,2]**2)/Rho[0])
    yp = y_coo*Utau*Rho[0]/Mu[0]

    Rey = -Rho*Urey.T[2]

    dy = np.gradient(y_coo,edge_order=2)
    dudy = np.gradient(Ut,dy,edge_order=2)
    dudy = (MuSGS+Mu)*dudy

    # print Mu

    # data_outVTK = computeVarGradient(data_outVTK,'U_AVG')
    # data1D = fct.getSlice(data_outVTK, ori, vec_tan)

    # [Vcoords, GradU] = fct.getArrayFromPointData(data1D, ['GRAD_U_AVG'])

    # y_coo = np.array(fct.getScalarCoord(Vcoords, 1))
    # id_sort = y_coo.argsort()
    # GradU = GradU[id_sort].T
    # dudy = Mu*GradU[1]

    # function output
    print ''
    return [ori, yp, Rey, dudy, Tw]

def computeVarGradient(data_outVTK,var):

    # test if the field "var" is present
    if data_outVTK.GetPointData().HasArray(var)!=1:
        raise ValueError("Error : field %s not present"%var)

    gradientFilter = vtk.vtkGradientFilter()
    gradientFilter.SetInputData(data_outVTK)
    gradientFilter.SetInputArrayToProcess(0,0,0,0,var)
    gradientFilter.SetResultArrayName('GRAD_%s'%var)
    gradientFilter.Update()
    data_outVTK = gradientFilter.GetOutput()

    return data_outVTK

def getPressureGradientAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, y_perc, var, Uinf=None, Rhoinf=None):

    # test if the field U_REY is present
    if data_outVTK.GetPointData().HasArray('GRAD_%s'%var)!=1:
        raise ValueError("Error : field U_REY not present")

    # test if the field TauWallAvg is present
    if wall_outVTK.GetPointData().HasArray('TauWallAvg')!=1:
        raise ValueError("Error : field TauWallAvg not present")

    # get the vector used to slice the field data
    [ori, vec_tan] = fct.getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    wall1D = fct.getSlice(wall_outVTK, ori, vec_tan)

    # extract the shear stress
    [_, Tw] = fct.getArrayFromPointData(wall1D, ['TauWallAvg'])

    data1D = fct.getSlice(data_outVTK, ori, vec_tan)

    [Vcoords, GradP, P] = fct.getArrayFromPointData(data1D, ['GRAD_%s'%var,var])

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(fct.getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()
    y_coo = y_coo[id_sort]
    GradP = GradP[id_sort].T
    P = P[id_sort]

    [_, _, deltaS, _] = fct.getDeltaAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf, Rhoinf)

    magTw = np.sqrt(Tw[0,0]**2+Tw[0,1]**2+Tw[0,2]**2)

    y_pos = util.find_nearest(y_coo, y_perc*(y_coo[-1]-y_coo[0]))

    scaling = (deltaS/magTw)

    BetaX = (deltaS/magTw)*GradP[0][y_pos]
    BetaY = (deltaS/magTw)*GradP[1][y_pos]

    ori[1] = y_coo[y_pos]

    return [ori, BetaX, BetaY, P[y_pos], scaling]

def getPressureGradientBetweenPosition(data_outVTK, wall_outVTK, cord_choice, x_p0, x_p1, Npts, y_perc, var, Uinf=None, Rhoinf=None):

    # function display 
    print '---- DAEPy::getPressureGradientBetweenPosition ----'

    data_outVTK = computeVarGradient(data_outVTK,var)

    Beta = np.zeros((Npts,2))
    P = np.zeros(Npts)
    scaling = np.zeros(Npts)
    pos = np.zeros((Npts, 3))

    for i in range(Npts):
        
        x_p_temp = x_p0 + (x_p1-x_p0)/(Npts-1)*i

        [pos[i,:], Beta[i,0], Beta[i,1], P[i], scaling[i]] = getPressureGradientAtPosition(data_outVTK, wall_outVTK, cord_choice, x_p_temp, y_perc, var, Uinf, Rhoinf)

    return [pos, Beta.T, P, scaling]

def getViscosity(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf=None, Rhoinf=None):
    # function display 
    print '---- DAEPy::getViscosity ----'

    # test if the field RHO_AVG is present
    if data_outVTK.GetPointData().HasArray('MU_LAM_AVG')!=1:
        raise ValueError("Error : field MU_LAM_AVG not present")

    # test if the field RHO_AVG is present
    if data_outVTK.GetPointData().HasArray('MU_SGS_AVG')!=1:
        raise ValueError("Error : field MU_SGS_AVG not present")

    # get the vector used to slice the field data
    [ori, vec_tan] = fct.getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    # data1D = fct.getSlice(data_outVTK, ori, vec_tan)
    data1D = fct.getSlice(data_outVTK, ori, vec_tan)

    [Vcoords, Mu, MuSGS, Rho] = fct.getArrayFromPointData(data1D, ['MU_LAM_AVG','MU_SGS_AVG','RHO_AVG'])

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(fct.getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()
    
    # sort cordinates along the line orthogonal to the wall
    y_coo = y_coo[id_sort]

    # sort the other quantities according to the coordinates y_coo. Velocity at the wall Ut[0] is
    # set to 0 to respect no slip condition
    Mu = Mu[id_sort]
    MuSGS = MuSGS[id_sort]
    Rho = Rho[id_sort]

    [_, Cf, Utau, _, _] = fct.getSkinFrictionAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf, Rhoinf)

    yp = y_coo*Utau*Rho[0]/Mu[0]

    return [ori, y_coo, yp, Mu, MuSGS]




def getTurbBudget(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf=None, Rhoinf=None):

    # function display 
    print '---- DAEPy::getTurbBudget ----'

    # test if the field RHO_AVG is present
    if data_outVTK.GetPointData().HasArray('RHO_AVG')!=1:
        raise ValueError("Error : field RHO_AVG not present")

    # test if the field RHO_AVG is present
    if data_outVTK.GetPointData().HasArray('RHO_RMS')!=1:
        raise ValueError("Error : field RHO_RMS not present")

    # test if the field U_AVG is present
    if data_outVTK.GetPointData().HasArray('U_AVG')!=1:
        raise ValueError("Error : field U_AVG not present")

    # test if the field U_RMS is present
    if data_outVTK.GetPointData().HasArray('U_RMS')!=1:
        raise ValueError("Error : field U_RMS not present")

    # test if the field U_REY is present
    if data_outVTK.GetPointData().HasArray('U_REY')!=1:
        raise ValueError("Error : field U_REY not present")

    # test if the field TauWallAvg is present
    if wall_outVTK.GetPointData().HasArray('TauWallAvg')!=1:
        raise ValueError("Error : field TauWallAvg not present")

    # get the vector used to slice the field data
    [ori, vec_tan] = fct.getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    # data1D = fct.getSlice(data_outVTK, ori, vec_tan)
    wall1D = fct.getSlice(wall_outVTK, ori, vec_tan)

    # extract the shear stress
    [_, Tw] = fct.getArrayFromPointData(wall1D, ['TauWallAvg'])

    [Vcoords, U, Urms, Urey, Rho, Rhorms, Mu] = fct.getArrayFromPointData(data_outVTK, ['U_AVG','U_RMS','U_REY','RHO_AVG','RHO_RMS','MU_LAM_AVG'])

    end = 200 ; nstep = 3

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(fct.getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()
    y_coo = y_coo[id_sort]
    Urms = (Urms[id_sort]**2).T
    Urey = Urey[id_sort].T
    U = U[id_sort].T
    Rho = Rho[id_sort]
    Rhorms = Rhorms[id_sort]
    Mu = Mu[id_sort]

    Utau = np.sqrt(np.sqrt(Tw[0,0]**2+Tw[0,1]**2+Tw[0,2]**2)/Rho[0])
    yp = y_coo*Utau*Rho[0]/Mu[0]

    # vec = (Rho*U[0]+Rhorms*Urms[0])/Rho
    vec = Urms[0]/Rho

    dy = np.gradient(y_coo[0:end:nstep],edge_order=2)
    # dy = np.ediff1d(yp)
    dudy = np.gradient(vec[0:end:nstep],dy,edge_order=2)

    Pr = Rho[0:end:nstep]*Urey[2][0:end:nstep]*dudy
    scaling = Rho[0]**2*Utau**4/Mu[0]

    return [ori, yp[0:end:nstep], Pr, scaling]



def getTurbulenceSpectrum(data_outVTK, cord_choice, x_perc):

    # function display 
    print '---- DAEPy::getTurbulenceSpectrum ----'

    # test if the field RHO_AVG is present
    if data_outVTK.GetPointData().HasArray('U_AVG')!=1:
        raise ValueError("Error : field U_AVG not present")

    # test if the field U_RMS is present
    if data_outVTK.GetPointData().HasArray('U_RMS')!=1:
        raise ValueError("Error : field U_RMS not present")

    # test if the field U_REY is present
    if data_outVTK.GetPointData().HasArray('U_REY')!=1:
        raise ValueError("Error : field U_REY not present")

    # # get the vector used to slice the field data
    # [ori, vec_tan] = fct.getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # print ori, vec_tan

    # # slice the 2D field data using the tangential vector to the wall
    # data1D = fct.getSlice(data_outVTK, ori, vec_tan)

    # extract the tangential velocity
    [Vcoords, U, Uavg, Urms, Urey, Rho, Rhoavg, P] = fct.getArrayFromPointData(data_outVTK, ['U','U_AVG','U_RMS','U_REY','RHO','RHO_AVG','T'])

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(fct.getScalarCoord(Vcoords, 2))
    id_sort = y_coo.argsort()
    y_coo = y_coo[id_sort]
    Urms = (Urms[id_sort]**2).T
    Urey = Urey[id_sort].T
    U = U[id_sort].T
    Uavg = Uavg[id_sort].T
    Rho = Rho[id_sort]
    Rhoavg = Rhoavg[id_sort]
    P = P[id_sort]

    UArray = (U[0]-U[0].mean())
    VArray = (U[1]-U[1].mean())
    WArray = (U[2]-U[2].mean())
    RhoArray = Rho-Rho.mean()
    PArray = P-P.mean()

    Dz = np.ediff1d(y_coo).mean()

    [kp, Epp, Rpp] = get1DSpectrum(PArray, Dz)
    [kr, Err, Rrr] = get1DSpectrum(RhoArray, Dz)
    [ku, Euu, Ruu] = get1DSpectrum(UArray, Dz)
    # [kv, Evv, Rvv] = get1DSpectrum(VArray, Dz)
    # [kw, Eww, Rww] = get1DSpectrum(WArray, Dz)


    idmax = np.argmax(Euu)

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.plot(ku,Euu/Euu[1],'b-')
    # ax.plot(kv,Evv/Evv[1],'r-')
    # ax.plot(kw,Eww/Eww[1],'k-')
    ax.plot(kp,Epp/Epp[1],'r-')
    ax.plot(kr,Err/Err[1],'k-')
    # ax.plot(kz,Euu[idmax]/Euu[1]*(kz/kz[idmax])**(-5./3.),'k--')
    ax.plot(ku,10*(ku/1000)**(-5./3.),'k--')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.show()

    # Ruu = np.correlate(Urms[0],Urms[0],mode='same')
    # plt.plot(Ruu)
    # plt.show()
    # nz = Ruu.shape[0]
    # kz = np.asarray([0]+range(0, nz/2))
    # kz = np.fft.fftfreq(nz,d=1.42857316e-05)
    # FFT = np.fft.rfft(0.5*Ruu)

    # nz = U[0].shape[0]
    # kz = np.asarray(range(0, nz/2)+[0] + range(-nz/2+1, 0))
    # FFT = np.fft.rfftn(U[0])

    # nz = myArray.shape[0]

    # kz = np.asarray(range(0, nz/2)+[0] + range(-nz/2+1, 0))

    # FFT = np.fft.rfftn(myArray)

    # nbins = 128

    # size = FFT.shape[0]/2
    # if nbins is None:
    #     nbins = size
    # dk = float(size)/float(nbins)
    # spectrum_energy = np.zeros(nbins)
    # samples = np.zeros(nbins)

    # for k in range(nz):
    #     kmod = float(kz[k]**2)**0.5
    #     if kmod < size:
    #         index = int(kmod/dk)
    #         k_approx = (index+0.5)*dk
    #         # print index,size,kz,ky_tmp,kx_tmp,myArray.shape
    #         samples[index] += 1
    #         spectrum_energy[index] += \
    #             np.abs(FFT[k])**2 * \
    #             (np.pi*4*k_approx**2*dk)

    # npts = (2.*nz)

    # for k in range(nbins):
    #     spectrum_energy[k] /= samples[k]*npts**2

    # kmod_list = np.arange(0, size, dk)

    return

def getTurbulenceTrianglePiroz(yp,Urms):

    IIb = np.zeros((len(yp)))
    IIIb = np.zeros((len(yp)))
    Ni = np.zeros((len(yp)))
    Ei = np.zeros((len(yp)))

    for idy in range(len(yp)):

        k = (Urms[0,idy]+Urms[1,idy]+Urms[2,idy])

        bij = np.matrix(([[Urms[0,idy]/k-1/3., 0., 0.],\
                         [0., Urms[1,idy]/k-1/3., 0.],\
                         [0., 0., Urms[2,idy]/k-1/3.]]))

        # ubji = np.triu(bij**2,0)
        # ubji = np.triu(np.multiply(bij,bij),0)

        # IIb[idy] = -0.5*np.sum(ubji)
        # IIIb[idy] = la.det(bij)
        IIb[idy] = -(bij[0,0]**2+bij[1,1]**2+bij[2,2]**2)/2.

        IIIb[idy] = (bij[0,0]**3+bij[1,1]**3+bij[2,2]**3)/3.

        Ni[idy] = np.sqrt(-1./3.*IIb[idy])
        Ei[idy] = np.power(np.abs(IIIb[idy])*0.5,1./3.)

    # function output
    print ''
    return [Ei, Ni, IIIb, IIb]

