# -*- coding: utf-8 -*-

import vtk
import os
import numpy as np
import matplotlib.pyplot as plt
import modTurb as mod

from vtk.util.numpy_support import vtk_to_numpy

def getOutputVTKwithPointDataFromFile(fileName):

    """
    Uses the correct VTK reader to obtain a VTK output object so that one can
    work with the data in *fileName*.
    NOTE : the cell data are transfered to the nodes.  

    Parameters
    ----------
    fileName : string
        The .vt* file to open

    Returns
    -------
    data_outVTK : VTK output object from a VTK reader
        VTK output object associated to the file type of *fileName*

    """

    # function display 
    print '---- DAEPy::getOutputVTKwithPointDataFromFile ----'
    
    # test if the file exists
    print '--> Reading', fileName
    if not os.path.isfile(fileName):
        raise ValueError("Error : file does not exists")

    extension = os.path.splitext(fileName)[-1]
    if extension == '.vtu':
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif extension == '.pvtu':
        reader = vtk.vtkXMLPUnstructuredGridReader()
    elif extension == '.vtp':
        reader = vtk.vtkXMLPolyDataReader()
    elif extension == '.vtm':
        # TODO : To check
        reader = vtk.vtkXMLMultiBlockDataReader()
        reader = vtk.MergeBlocks(reader)
    else:
        raise ValueError("Error: unknown extension of file "+fileName)

    reader.SetFileName(fileName)
    reader.Update()
    data_outVTK = reader.GetOutput()

    # All the data are transfered to the nodes
    c2p = vtk.vtkCellDataToPointData()
    c2p.SetInputData(data_outVTK)
    c2p.Update()
    data_outVTK = c2p.GetOutput()

    # list the fields available
    n_fields = data_outVTK.GetPointData().GetNumberOfArrays()
    print '--> Available:', n_fields, 'fields'
    for i in range(n_fields):
        print '        -', data_outVTK.GetPointData().GetArrayName(i)


    print ''
    return data_outVTK


def writeDataToFile(data, fileName):

    if not hasattr(data, 'GetClassName'):
        raise ValueError('data is not a vtk data obj')
    if data.GetClassName() == 'vtkUnstructuredGrid':
        writer = vtk.vtkXMLUnstructuredGridWriter()
    else:
        raise ValueError('data format not implemented')

    writer.SetInputData(data)
    writer.SetFileName(fileName)
    writer.Write()


def getArrayFromPointData(data_outVTK, field_name):
    
    """
    The string list *field_name* command the data contained in the VTK output object *data_outVTK* 
    that are going to be stored the returned array.

    WARNING : USE ONLY FOR POINTS DATA, for cell data, use getArrayFromCellData()

    Parameters
    ----------
    field_name : list of string
        The fields names, which would be for this BE :

    Returns
    -------
    coord : numpy array

    data_arr : numpy array

    """

    # function display 
    print '---- DAEPy::getArrayFromPointData ----'

    coord = np.array(
        [data_outVTK.GetPoint(i) for i in range(data_outVTK.GetNumberOfPoints())],
        dtype=np.float32)

    print '--> extract fields', [f for f in field_name]
    data_arr = [vtk_to_numpy(data_outVTK.GetPointData().GetArray(f)) for f in field_name]


    print ''
    return [coord] + data_arr


def getArrayFromCellData(data_outVTK, field_name):
    
    """
    The string list *field_name* command the data contained in the VTK output object *data_outVTK* 
    that are going to be stored the returned array.

    WARNING : USE ONLY FOR CELL DATA, for point data, use getArrayFromPointData()

    Parameters
    ----------
    field_name : list of string
        The fields names, which would be for this BE :

    Returns
    -------
    coord : numpy array

    data_arr : numpy array

    """

    # function display 
    print '---- DAEPy::getArrayFromCellData ----'
    print 'FUNCTION NOT TESTED YET ---- TO DEBUG'

    coord = np.array(
        [data_outVTK.GetCell(i) for i in range(data_outVTK.GetNumberOfCells())],
        dtype=np.float32)

    print '--> extract fields', [f for f in field_name]
    data_arr = [vtk_to_numpy(data_outVTK.GetCellData().GetArray(f)) for f in field_name]


    print ''
    return [coord] + data_arr


def getSlice(data_outVTK, orig, nor):
    
    """
    Defines the plane normal to the vector *nor* that contains the point *orig* 
    and then extract the data from the reader *reader*. The dimension of the data
    returned is reduced by one (a 2D data field gives a 1D data filed).

    Parameters
    ----------
    data_outVTK : VTK output object from a VTK reader
        VTK output object associated to the file we want to extract a slice of data.

    orig : tuple(3)
        A point element of the plane used to slice the data.

    normal : tuple(3)
        A vector normal to the plane used to slice the data.

    Returns
    -------
    slice_outVTK : VTK output object from a VTK reader
        Returned data are still not in an array form, but a VTK output object 
        is returned where the data are stored.
        
    """

    #function display
    print '---- DAEPy::getSlice ----'

    # Plane of the Slice
    plane = vtk.vtkPlane()
    plane.SetOrigin(orig)
    plane.SetNormal(nor)
    print '--> normal used to slice: ', nor

    # Cutter plane
    planeCut = vtk.vtkCutter()
    planeCut.SetInputData(data_outVTK)
    planeCut.SetCutFunction(plane)

    # Make the slice
    planeCut.Update()
    slice_outVTK = planeCut.GetOutput()

    print ''
    return slice_outVTK


def getLine(data_outVTK, orig, n1, n2):

    """
    Uses twice the getSlice() function to extract the line defined with the point *orig*
    and the vector normal to the plane defined with *n1* and *n2*.

    Parameters
    ----------
    data_outVTK : VTK output object from a VTK reader
        VTK output object associated to the file we want to extract a slice of data.

    orig : tuple(3)
        A point element of the line where to extract the data.

    n1 : tuple(3)
        A vector normal to the line where to extract the data.

    n2 : tuple(3)
        A vector normal to the line where to extract the data.

    Returns
    -------
    line_outVTK : VTK output object from a VTK reader
        Returned data are still not in an array form, but a VTK output object
        is returned where the data are stored.
        
    """

    # function display 
    print '---- DAEPy::getLine ----'

    # stop execution if data not consistent with the method
    if data_outVTK.GetCell(0).GetNumberOfPoints() < 4:
        raise ValueError("Error: cells in data from VTK output object are not 3D cells, be sure the data used here are 3D.")

    # Double slicing
    print '--> 1st slicing...'
    dataSlice1 = getSlice(data_outVTK, orig, n1)
    print '--> 2nd slicing...'
    dataSlice2 = getSlice(dataSlice1, orig, n2)
    
    print ''
    return dataSlice2


def getLineTangentialVector(line_outVTK, x_perc, cord_choice):

    """
    Compute the tangential vector at point X0 to the line meshed in the 1D file *wall_outVTK*.

    We define a cord L of the line in the *line_outVTK* file thanks to the maximum distant 
    between two cells center in the *cord_choice* direction. Then, the coordinates of the cell center 
    and the cell tangential vector from the cell at *x_perc* percent of this cord are returned.
    
    WARNING : the data in the VTK output object *line_outVTK* must be 1D, meaning that the cells 
              in it are defined with only 2 points (not 4 or 8 points), so that the tangential 
              vector definition is not ambiguous.

    Parameters
    ----------
    line_outVTK : VTK output object from a VTK reader
        The cell in this file are 1D cell, meaning they should be defined with only 2 points.
        We test the CellType of the first cell to check if we have 1D cells.      

    x_perc : float
        Gives the percentage of the cord where we want to select a cell

    cord_choice : int
        Gives the direction (0->x, 1->y, 2->z), to use to define the cord of the line.

    Returns
    -------
    ori_out : tuple(3)
        Coordinates of the cell center that defined the returned tangential vector.

    nor_out : tuple(3)
        Coordinates of the tangential vector.

    """

    # function display 
    print '---- DAEPy::getLineTangentialVector ----'

    # stop execution if data not consistent with the method see VTK_file_format.pdf for a list of cell type
    if line_outVTK.GetCellType(0) > 4: 
        raise ValueError("Error: cells in data from VTK output object are not 1D cells, be sure the data are 1D or sliced from a 2D data set.")

    # number of cell
    n_cell = line_outVTK.GetNumberOfCells()
    
    # we extract the boundary of each 2D cell to compute the tangential vector
    a = np.array([line_outVTK.GetCell(i).GetPoints().GetPoint(0) for i in range(n_cell)])
    b = np.array([line_outVTK.GetCell(i).GetPoints().GetPoint(1) for i in range(n_cell)])
    tangent_vector = b-a

    # centers of cells
    Centers = (a+b)/2

    # compute the max and min value of the cell centers in the chosen direction
    dMin = min(Centers[:, cord_choice])
    dMax = max(Centers[:, cord_choice])

    Coord = dMin + x_perc*(dMax-dMin)

    iMin = (abs(Centers[:, cord_choice]-Coord)).argmin()

    tangent_vector = tangent_vector[iMin]/np.linalg.norm(tangent_vector[iMin])
    orig = Centers[iMin]

    print '--> Coordinates of the point selected: ', orig

    print ''
    return [orig, tangent_vector]


def getCommonFieldAndWallSlice(data_outVTK, wall_outVTK, orig, nor):

    """
    This function uses simple slicing function, yet it assures the user to work
    with consistent field and wall data in *data_outVTK* and *wall_outVTK*.
    WARNING : in a standard use, this function is conceived to sliced a 3D field
              and a 2D wall field.

    Parameters
    ----------
    data_outVTK : VTK output object from a VTK reader
        The field data, normally a (huge) 3D field

    wall_outVTK : VTK output object from a VTK reader
        The wall data associated to the field data

    orig : tuple(3)
        A point element of the plane used to slice the data.

    nor : tuple(3)
        A vector normal to the plane used to slice the data.
    
    Returns
    -------
    dataSl_outVTK : VTK output object from a VTK reader
        Returned data are still not in an array form, but a VTK output object 
        is returned where the data are stored.

    wallSl_outVTK : VTK output object from a VTK reader
        Returned data are still not in an array form, but a VTK output object 
        is returned where the data are stored.

    """

    # function display 
    print '---- DAEPy::getCommonFieldAndWallSlice ----'

    # slice the data with the same plane
    print '--> Field and Wall slicing'
    dataSl_outVTK = getSlice(data_outVTK, orig, nor)
    wallSl_outVTK = getSlice(wall_outVTK, orig, nor)

    # function output
    print ''
    return [dataSl_outVTK, wallSl_outVTK]


def getScalarCoord(Vcoords, dir):

    """
    This function gives a parametric description of the points in *Vcoords*.
    The parameter *dir* is used to gives the reference direction to describe the points.

    Parameters
    ----------
    Vcoords : vector of tuple(3)
        A list of points on the same line.

    Returns
    -------
    X_off : vector of float
        A vector of float that can gives the distance to the original points in *Vcoords*.

    """

    # function display 
    print '---- DAEPy::getScalarCoord ----'

    X = Vcoords[:,dir]
    V_off = Vcoords - Vcoords[X.argsort()[0],:]
    X_off = [np.linalg.norm(V_off[i,:]) for i in range(len(V_off[:,0]))]

    # function output
    print ''
    return X_off


def getSkinFrictionAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf=None, Rhoinf=None):

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

    Cf : float
        The skin friction coefficient.

    Utau : float
        The shear velocity.

    """

    # function display 
    print '---- DAEPy::getSkinFrictionAtPosition ----'

    # test if the field U_AVG is present
    if data_outVTK.GetPointData().HasArray('U_AVG')!=1:
        raise ValueError("Error : field U_AVG not present")

    # test if the field RHO_AVG is present
    if data_outVTK.GetPointData().HasArray('RHO_AVG')!=1:
        raise ValueError("Error : field RHO_AVG not present")

    # test if the field MU_LAM_AVG is present
    if data_outVTK.GetPointData().HasArray('MU_LAM_AVG')!=1:
        raise ValueError("Error : field MU_LAM_AVG not present")

    # test if the field TauWallAvg is present
    if wall_outVTK.GetPointData().HasArray('TauWallAvg')!=1:
        raise ValueError("Error : field TauWallAvg not present")

    # stop execution if data not consistent with the method see VTK_file_format.pdf for a list of cell type
    if (wall_outVTK.GetCellType(0) > 4): 
        raise ValueError("Error: cells in data from wall VTK output object are not 1D cells, be sure the data are 1D or sliced from a 2D data set.")

    # TODO : Find a way to bound data1D (or output of getArrayFromPointData) to a certain length
    # Could be usefull when post-processing data in very large domain (airfoil cases for example)

    # get the vector used to slice the field data
    [ori, vec_tan] = getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    data1D = getSlice(data_outVTK, ori, vec_tan)
    wall1D = getSlice(wall_outVTK, ori, vec_tan)

    # extract the tangential velocity
    [Vcoords, Mu, Rho, U] = getArrayFromPointData(data1D, ['MU_LAM_AVG','RHO_AVG','U_AVG'])
    Ut = np.sum(U*vec_tan, axis=1)

    # extract the shear stress
    [_, Tw] = getArrayFromPointData(wall1D, ['TauWallAvg'])

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()
    y_coo = y_coo[id_sort]
    Rho = Rho[id_sort]
    Ut = Ut[id_sort]
    Mu = Mu[id_sort]

    # Automatic definition of the external (free stream) velocity based on the extracted velocity if Uext = None
    if Uinf == None:
        # Uinf = Ut[-1]
        Uinf = np.mean(Ut[-20:-1])

    # Automatic definition of the external (free stream) density based on the extracted density if RhoExt = None
    if Rhoinf == None:
        # Rhoinf = Rho[-1]
        Rhoinf = np.mean(Rho[-20:-1])

    Cf = np.sqrt(Tw[0,0]**2+Tw[0,1]**2+Tw[0,2]**2)/(0.5*Rhoinf*Uinf**2)
    Utau = np.sqrt(np.sqrt(Tw[0,0]**2+Tw[0,1]**2+Tw[0,2]**2)/Rho[0])
    Tw = np.sqrt(Tw[0,0]**2+Tw[0,1]**2+Tw[0,2]**2)

    # [_, _, _, dudy, Tw] = mod.getStressAtPosition(data_outVTK, wall_outVTK, x_perc, cord_choice, Uinf, Rhoinf)
    # Utau = np.sqrt(dudy[4]/Rho[0])

    Yplus = ((y_coo[1]-y_coo[0])*Rho[0]*Utau)/Mu[0]


    # function output
    print ''
    return [ori, Cf, Utau, Yplus, Tw]


def getReynoldsAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf=None, Rhoinf=None, Muinf=None):
    """
    Gives the Reynolds number based on theboundary layer thickness :math:`\delta`, displacement :math:`\delta^*`
    and momentum :math:`\\theta` thickness in x0. 
    How is define x0 : let C be the curve in the 3D space described by the 1D cells 
    in *wall_outVTK*. The cells center are denoted with Xc_i. The curve C is include in 
    a cubic box of bounds Xmin, Xmax, Ymin, Ymax, Zmin, Zmax. Depending on the *cord_choice* 
    value choice, we look for the cell of C which center coordinates Xc_i[cord_choice] is 
    the closest from (Xmin + x_perc*(Xmax-Xmin)).
    This definition is based on an intelligent choice for the reference axis.

    Example 2D si *cord_choice*=0 and x_perc=0.6, the curve is between the 2 o :

              Ymax _________________________  __o
                                             /  |
                            _______,---X~~~~¨   |
              Ymin ___ o___/           |        |
                       |               |        |
                       |               |        |
                       |               |        |
                       |               |        |
                     Xmin           THE CELL   Xmax
                     
    Parameters
    ----------
    data_outVTK : VTK output object from a VTK reader
        This field data MUST fulfill 2 conditions:
            - contain the field "U_AVG",
            - contain the field "RHO_AVG",
            - represent a 2D field data, as we want to extract a 1D velocity
            profile from it.

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
        (optional) Free stream velocity, default value is the farthest point from the wall at x_perc.

    Rhoinf : float
        (optional) Free stream velocity, default value is the farthest point from the wall at x_perc.

    Muinf : float
        (optional) Free stream viscosity, default value is the farthest point from the wall at x_perc.

    Returns
    -------
    ori : tuple(3)
        The point coordinates of the point at the wall where we extract the data.

    ReD : vector of floats
        Reynolds number based on the BL thickness.

    ReDS : vector of floats
        Reynolds number based on the compressible BL displacement thickness.

    ReO : vector of floats
        Reynolds number based on the compressible BL momentum thickness.

    ReT : vector of floats
        Reynolds number based on the shear velocity and BL thickness.

    """

    # function display 
    print '---- DAEPy::getReynoldsAtPosition ----'

    # test if the field MU_LAM_AVG is present
    if data_outVTK.GetPointData().HasArray('MU_LAM_AVG')!=1:
        raise ValueError("Error : field MU_LAM_AVG not present")

    # test if the field RHO_AVG is present
    if data_outVTK.GetPointData().HasArray('RHO_AVG')!=1:
        raise ValueError("Error : field RHO_AVG not present")

    # test if the field U_AVG is present
    if data_outVTK.GetPointData().HasArray('U_AVG')!=1:
        raise ValueError("Error : field U_AVG not present")

    # stop execution if data not consistent with the method see VTK_file_format.pdf for a list of cell type
    if (wall_outVTK.GetCellType(0) > 4): 
        raise ValueError("Error: cells in data from wall VTK output object are not 1D cells, be sure the data are 1D or sliced from a 2D data set.")

    # get the vector used to slice the field data
    [ori, vec_tan] = getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    data1D = getSlice(data_outVTK, ori, vec_tan)

    # extract the tangential velocity
    [Vcoords, Mu, Rho, U] = getArrayFromPointData(data1D, ['MU_LAM_AVG','RHO_AVG','U_AVG'])
    Ut = np.sum(U*vec_tan, axis=1)

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()

    # sort cordinates along the line orthogonal to the wall
    y_coo = y_coo[id_sort]

    # sort the other quantities according to the coordinates y_coo. Velocity at the wall Ut[0] is
    # set to 0 to respect no slip condition
    Mu = Mu[id_sort]
    Rho = Rho[id_sort]
    Ut = Ut[id_sort] ; Ut[0] = 0

    # Automatic definition of the external (free stream) viscosity based on the extracted viscosity if Muinf = None
    if Muinf == None:
        Muinf = np.mean(Mu[-20:-1])
        # Muinf = Mu[-1]

    # Automatic definition of the external (free stream) density based on the extracted density if Rhoinf = None
    if Rhoinf == None:
        Rhoinf = np.mean(Rho[-20:-1])
        # Rhoinf = Rho[-1]

    # Automatic definition of the external (free stream) velocity based on the extracted density if Uinf = None
    if Uinf == None:
        Uinf = np.mean(Ut[-20:-1])
        # Uinf = Ut[-1]

    [_, delta, deltaS, theta] = getDeltaAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf, Rhoinf)
    [_, Cf, Utau, _, _] = getSkinFrictionAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf, Rhoinf)
    
    ReD = Rhoinf*Uinf*delta/Muinf
    ReDS = Rhoinf*Uinf*deltaS/Muinf
    ReO = Rhoinf*Uinf*theta/Muinf
    ReT = Rho[0]*Utau*delta/Mu[0]

    # function output
    print ''
    return [ori, ReD, ReDS, ReO, ReT]

def getIncDeltaAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf=None, Rhoinf=None):

    """
    Gives the boundary layer thickness :math:`\delta`, displacement :math:`\delta^*` and
    momentum :math:`\\theta` thickness in x0. 
    How is define x0 : let C be the curve in the 3D space described by the 1D cells 
    in *wall_outVTK*. The cells center are denoted with Xc_i. The curve C is include in 
    a cubic box of bounds Xmin, Xmax, Ymin, Ymax, Zmin, Zmax. Depending on the *cord_choice* 
    value choice, we look for the cell of C which center coordinates Xc_i[cord_choice] is 
    the closest from (Xmin + x_perc*(Xmax-Xmin)).
    This definition is based on an intelligent choice for the reference axis.

    Example 2D si *cord_choice*=0 and x_perc=0.6, the curve is between the 2 o :

              Ymax _________________________  __o
                                             /  |
                            _______,---X~~~~¨   |
              Ymin ___ o___/           |        |
                       |               |        |
                       |               |        |
                       |               |        |
                       |               |        |
                     Xmin           THE CELL   Xmax
                     
    Regarding the calculation of the displacement :math:`\delta^*` and momentum :math:`\theta` thickness,
    the compressible flow definition based on the mass flow rate is used:

    .. math:: \delta^* = \int_{0}^{\infty } \\left(1 - \\frac{\\rho (y) u(y) }{\\rho_0 u_0 } \\right) \\mathrm{d}y
    
    and

    .. math:: \\theta = \int_{0}^{\infty } \\frac{\\rho (y) u(y) }{\\rho_0 u_0 } \\left(1 - \\frac{u(y) }{u_0 } \\right) \\mathrm{d}y

    Where :math:`\\rho_0` and :math:`u_0` are the density and velocity in the '*free stream*' outside the boundary layer, and :math:`y` is the coordinate normal to the wall.

    Parameters
    ----------
    data_outVTK : VTK output object from a VTK reader
        This field data MUST fulfill 2 conditions:
            - contain the field "U_AVG",
            - contain the field "RHO_AVG",
            - represent a 2D field data, as we want to extract a 1D velocity
            profile from it.

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
        (optional) Free stream velocity, default value is the farthest point from the wall at x_perc.

    Rhoinf : float
        (optional) Free stream velocity, default value is the farthest point from the wall at x_perc.

    Returns
    -------
    ori : tuple(3)
        The point coordinates of the point at the wall where we extract the data.

    delta : float
        The boundary layer thickness.

    deltaS : float
        The compressible boundary layer displacement thickness.

    theta : float
        The compressible boundary layer momentum thickness.

    """

    # function display 
    print '---- DAEPy::getIncDeltaAtPosition ----'

    # test if the field U_AVG is present
    if data_outVTK.GetPointData().HasArray('U_AVG')!=1:
        raise ValueError("Error : field U_AVG not present")

    # stop execution if data not consistent with the method see VTK_file_format.pdf for a list of cell type
    if (wall_outVTK.GetCellType(0) > 4): 
        raise ValueError("Error: cells in data from wall VTK output object are not 1D cells, be sure the data are 1D or sliced from a 2D data set.")

    # TODO : Find a way to bound data1D (or output of getArrayFromPointData) to a certain length
    # Could be usefull when post-processing data in very large domain (airfoil cases for example)

    # get the vector used to slice the field data
    [ori, vec_tan] = getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    data1D = getSlice(data_outVTK, ori, vec_tan)

    # extract the tangential velocity
    [Vcoords, U] = getArrayFromPointData(data1D, ['U_AVG'])

    Ut = np.sum(U*vec_tan, axis=1)

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()

    # sort cordinates along the line orthogonal to the wall
    y_coo = y_coo[id_sort]

    # sort the other quantities according to the coordinates y_coo. Velocity at the wall Ut[0] is
    # set to 0 to respect no slip condition
    Ut = Ut[id_sort]; Ut[0] = 0

    # Automatic definition of the external (free stream) velocity based on the extracted velocity if Uinf = None
    if Uinf == None:
        Uinf = np.mean(Ut[-20:-1])
        # Uinf = Ut[-1]

    # find the index of the 2 points around 99% of Uext
    ordre = abs(Ut-0.99*Uinf).argsort()[0:2]
    ordre = np.sort(ordre)

    # compute the BL thickness using a simple linear model
    xa = y_coo[ordre[0]]
    ua = Ut[ordre[0]]
    xb = y_coo[ordre[1]] 
    ub = Ut[ordre[1]]
    tau = (ub-ua)/(xb-xa)
    delta = xa + (Uinf*0.99-ua)/tau
    # theta = np.trapz((Rho[0:ordre[1]]*Ut[0:ordre[1]])/(Rhoinf*Uinf)*(1-Ut[0:ordre[1]]/(Uinf)),y_coo[0:ordre[1]])
    # deltaS = np.trapz((1-(Rho[0:ordre[1]]*Ut[0:ordre[1]])/(Rhoinf*Uinf)),y_coo[0:ordre[1]])

    theta1 = np.trapz((Ut[0:ordre[0]+1])/(Uinf)*(1-Ut[0:ordre[0]+1]/(Uinf)),y_coo[0:ordre[0]+1])
    theta2 = np.trapz((Ut[0:ordre[1]+1])/(Uinf)*(1-Ut[0:ordre[1]+1]/(Uinf)),y_coo[0:ordre[1]+1])
    tau = (theta2-theta1)/(xb-xa)
    theta = theta1 + (delta-xa)*tau

    deltaS1 = np.trapz((1-(Ut[0:ordre[0]+1])/(Uinf)),y_coo[0:ordre[0]+1])
    deltaS2 = np.trapz((1-(Ut[0:ordre[1]+1])/(Uinf)),y_coo[0:ordre[1]+1])
    tau = (deltaS2-deltaS1)/(xb-xa)
    deltaS = deltaS1 + (delta-xa)*tau

    # print ordre
    # print y_coo


    # theta = np.trapz((Rho*Ut)/(Rhoinf*Uinf)*(1-Ut/(Uinf)),y_coo)
    # deltaS = np.trapz((1-(Rho*Ut)/(Rhoinf*Uinf)),y_coo)

    # fig = plt.figure(figsize=(10, 6))
    # ax = fig.add_subplot(111)

    # ax.plot(Ut,y_coo,'b-')
    # ax.plot(Ut[ordre[0]],y_coo[ordre[0]],'gs')
    # ax.plot(Ut[ordre[1]],y_coo[ordre[1]],'ks')
    # ax.plot(Uinf*0.99,delta,'rs')
    # ax.set_title(r'$x = %1.5f \/-\/ \theta_{i} = %1.6f \/-\/ \delta^{*} = %1.6f$'%(ori[0],theta,deltaS),y=1.02)

    # plt.show()



    # function output
    print ''
    return [ori, delta, deltaS, theta]


def getDeltaAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf=None, Rhoinf=None):

    """
    Gives the boundary layer thickness :math:`\delta`, displacement :math:`\delta^*` and
    momentum :math:`\\theta` thickness in x0. 
    How is define x0 : let C be the curve in the 3D space described by the 1D cells 
    in *wall_outVTK*. The cells center are denoted with Xc_i. The curve C is include in 
    a cubic box of bounds Xmin, Xmax, Ymin, Ymax, Zmin, Zmax. Depending on the *cord_choice* 
    value choice, we look for the cell of C which center coordinates Xc_i[cord_choice] is 
    the closest from (Xmin + x_perc*(Xmax-Xmin)).
    This definition is based on an intelligent choice for the reference axis.

    Example 2D si *cord_choice*=0 and x_perc=0.6, the curve is between the 2 o :

              Ymax _________________________  __o
                                             /  |
                            _______,---X~~~~¨   |
              Ymin ___ o___/           |        |
                       |               |        |
                       |               |        |
                       |               |        |
                       |               |        |
                     Xmin           THE CELL   Xmax
                     
    Regarding the calculation of the displacement :math:`\delta^*` and momentum :math:`\theta` thickness,
    the compressible flow definition based on the mass flow rate is used:

    .. math:: \delta^* = \int_{0}^{\infty } \\left(1 - \\frac{\\rho (y) u(y) }{\\rho_0 u_0 } \\right) \\mathrm{d}y
    
    and

    .. math:: \\theta = \int_{0}^{\infty } \\frac{\\rho (y) u(y) }{\\rho_0 u_0 } \\left(1 - \\frac{u(y) }{u_0 } \\right) \\mathrm{d}y

    Where :math:`\\rho_0` and :math:`u_0` are the density and velocity in the '*free stream*' outside the boundary layer, and :math:`y` is the coordinate normal to the wall.

    Parameters
    ----------
    data_outVTK : VTK output object from a VTK reader
        This field data MUST fulfill 2 conditions:
            - contain the field "U_AVG",
            - contain the field "RHO_AVG",
            - represent a 2D field data, as we want to extract a 1D velocity
            profile from it.

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
        (optional) Free stream velocity, default value is the farthest point from the wall at x_perc.

    Rhoinf : float
        (optional) Free stream velocity, default value is the farthest point from the wall at x_perc.

    Returns
    -------
    ori : tuple(3)
        The point coordinates of the point at the wall where we extract the data.

    delta : float
        The boundary layer thickness.

    deltaS : float
        The compressible boundary layer displacement thickness.

    theta : float
        The compressible boundary layer momentum thickness.

    """

    # function display 
    print '---- DAEPy::getDeltaAtPosition ----'

    # test if the field U_AVG is present
    if data_outVTK.GetPointData().HasArray('U_AVG')!=1:
        raise ValueError("Error : field U_AVG not present")

    # test if the field RHO_AVG is present
    if data_outVTK.GetPointData().HasArray('RHO_AVG')!=1:
        raise ValueError("Error : field RHO_AVG not present")

    # stop execution if data not consistent with the method see VTK_file_format.pdf for a list of cell type
    if (wall_outVTK.GetCellType(0) > 4): 
        raise ValueError("Error: cells in data from wall VTK output object are not 1D cells, be sure the data are 1D or sliced from a 2D data set.")

    # TODO : Find a way to bound data1D (or output of getArrayFromPointData) to a certain length
    # Could be usefull when post-processing data in very large domain (airfoil cases for example)

    # get the vector used to slice the field data
    [ori, vec_tan] = getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    data1D = getSlice(data_outVTK, ori, vec_tan)

    # extract the tangential velocity
    [Vcoords, U, Rho] = getArrayFromPointData(data1D, ['U_AVG','RHO_AVG'])

    Ut = np.sum(U*vec_tan, axis=1)

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()

    # sort cordinates along the line orthogonal to the wall
    y_coo = y_coo[id_sort]

    # sort the other quantities according to the coordinates y_coo. Velocity at the wall Ut[0] is
    # set to 0 to respect no slip condition
    Ut = Ut[id_sort]; Ut[0] = 0
    Rho = Rho[id_sort]

    # Automatic definition of the external (free stream) velocity based on the extracted velocity if Uinf = None
    if Uinf == None:
        Uinf = np.mean(Ut[-20:-1])
        # Uinf = Ut[-1]

    # Automatic definition of the external (free stream) density based on the extracted density if Rhoinf = None
    if Rhoinf == None:
        Rhoinf = np.mean(Rho[-20:-1])
        # Rhoinf = Rho[-1]

    # find the index of the 2 points around 99% of Uext
    ordre = abs(Ut-0.99*Uinf).argsort()[0:2]
    ordre = np.sort(ordre)

    # compute the BL thickness using a simple linear model
    xa = y_coo[ordre[0]]
    ua = Ut[ordre[0]]
    xb = y_coo[ordre[1]] 
    ub = Ut[ordre[1]]
    tau = (ub-ua)/(xb-xa)
    delta = xa + (Uinf*0.99-ua)/tau
    # theta = np.trapz((Rho[0:ordre[1]]*Ut[0:ordre[1]])/(Rhoinf*Uinf)*(1-Ut[0:ordre[1]]/(Uinf)),y_coo[0:ordre[1]])
    # deltaS = np.trapz((1-(Rho[0:ordre[1]]*Ut[0:ordre[1]])/(Rhoinf*Uinf)),y_coo[0:ordre[1]])


    theta1 = np.trapz((Rho[0:ordre[0]+1]*Ut[0:ordre[0]+1])/(Rhoinf*Uinf)*(1-Ut[0:ordre[0]+1]/(Uinf)),y_coo[0:ordre[0]+1])
    theta2 = np.trapz((Rho[0:ordre[1]+1]*Ut[0:ordre[1]+1])/(Rhoinf*Uinf)*(1-Ut[0:ordre[1]+1]/(Uinf)),y_coo[0:ordre[1]+1])
    tau = (theta2-theta1)/(xb-xa)
    theta = theta1 + (delta-xa)*tau

    deltaS1 = np.trapz((1-(Rho[0:ordre[0]+1]*Ut[0:ordre[0]+1])/(Rhoinf*Uinf)),y_coo[0:ordre[0]+1])
    deltaS2 = np.trapz((1-(Rho[0:ordre[1]+1]*Ut[0:ordre[1]+1])/(Rhoinf*Uinf)),y_coo[0:ordre[1]+1])
    tau = (deltaS2-deltaS1)/(xb-xa)
    deltaS = deltaS1 + (delta-xa)*tau


    # theta = np.trapz((Rho*Ut)/(Rhoinf*Uinf)*(1-Ut/(Uinf)),y_coo)
    # deltaS = np.trapz((1-(Rho*Ut)/(Rhoinf*Uinf)),y_coo)

    # fig = plt.figure(figsize=(10, 6))
    # ax = fig.add_subplot(111)

    # ax.plot(Ut,y_coo,'b-')
    # ax.plot(Ut[ordre[0]],y_coo[ordre[0]],'gs')
    # ax.plot(Ut[ordre[1]],y_coo[ordre[1]],'ks')
    # ax.plot(Uinf*0.99,delta,'rs')
    # ax.set_title(r'$x = %1.5f \/-\/ \theta_{i} = %1.6f \/-\/ \delta^{*} = %1.6f$'%(ori[0],theta,deltaS),y=1.02)

    # plt.show()



    # function output
    print ''
    return [ori, delta, deltaS, theta]


def getReynoldsBetweenPosition(data_outVTK, wall_outVTK, cord_choice, x_p0, x_p1, Npts, Uinf=None, Rhoinf=None, Muinf=None):
    """
    Compute the boundary layer thickness for *Npts* equally distributed between the 2 position
    defined thanks to *x_p0*, *x_p1* and *cord_choice*. See the documentation of the function
    getDeltaAtPosition() for more information.

    Parameters
    ----------
    data_outVTK : VTK output object from a VTK reader
        This field data MUST fulfill 2 conditions:
            - contain the field "U_AVG",
            - represent a 2D field data, as we want to extract a 1D velocity
            profile from it.

    wall_outVTK : VTK output object from a VTK reader
        This field is supposed to describe a curve with 1D cell type, like in getLineTangentialVector()

    cord_choice : integer, 0, 1 or 2
        Gives the axis that is going to be used to define the point where the velocity
        profile is extracted. Convention :
            - 0 = x axis
            - 1 = y axis
            - 2 = z axis

    x_p0, x_p1 : float between 0 and 1
        gives the bound of the line where to compute the BL thickness.

    Npts: integer
        Number of points where to compute delta wanted.

    Uinf : float 
        (optional) Free stream velocity

    Rhoinf : float 
        (optional) Free stream density

    Muinf : float 
        (optional) Free stream viscosity

    Returns
    -------
    pos : vector of tuple(3)
        the coordinates of the points where delta have been computed.

    ReD : vector of floats
        Reynolds number based on the BL thickness at the *Npts* different points.

    ReDS : vector of floats
        Reynolds number based on the compressible BL displacement thickness at the *Npts* different points.

    ReO : vector of floats
        Reynolds number based on the compressible BL momentum thickness at the *Npts* different points.

    ReT : vector of floats
        Reynolds number based on BL thickness and shear velocity at the *Npts* different points.

    """

    # function display 
    print '---- DAEPy::getReynoldsBetweenPosition ----'

    ReD = np.zeros(Npts)
    ReO = np.zeros(Npts)
    ReDS = np.zeros(Npts)
    ReT = np.zeros(Npts)
    pos = np.zeros((Npts, 3))

    for i in range(Npts):
        
        x_p_temp = x_p0 + (x_p1-x_p0)/(Npts-1)*i

        [pos[i,:], ReD[i], ReDS[i], ReO[i], ReT[i]] = getReynoldsAtPosition(data_outVTK, wall_outVTK, cord_choice, x_p_temp, Uinf, Rhoinf, Muinf)

    return [pos, ReD, ReDS, ReO, ReT]


def getSkinFrictionBetweenPosition(data_outVTK, wall_outVTK, cord_choice, x_p0, x_p1, Npts, Uinf=None, Rhoinf=None):

    """
    Compute the boundary layer thickness for *Npts* equally distributed between the 2 position
    defined thanks to *x_p0*, *x_p1* and *cord_choice*. See the documentation of the function
    getDeltaAtPosition() for more information.

    Parameters
    ----------
    data_outVTK : VTK output object from a VTK reader
        This field data MUST fulfill 2 conditions:
            - contain the field "U_AVG",
            - represent a 2D field data, as we want to extract a 1D velocity
            profile from it.

    wall_outVTK : VTK output object from a VTK reader
        This field is supposed to describe a curve with 1D cell type, like in getLineTangentialVector()

    cord_choice : integer, 0, 1 or 2
        Gives the axis that is going to be used to define the point where the velocity
        profile is extracted. Convention :
            - 0 = x axis
            - 1 = y axis
            - 2 = z axis

    x_p0, x_p1 : float between 0 and 1
        gives the bound of the line where to compute the BL thickness.

    Npts: integer
        Number of points where to compute delta wanted.

    Uinf : float 
        (optional) Free stream velocity

    Rhoinf : float 
        (optional) Free stream density

    Returns
    -------
    pos : vector of tuple(3)
        the coordinates of the points where delta have been computed.

    Cf : vector of floats
        Skin friction coefficient at the *Npts* different points.

    Utau : vector of floats
        Shear velocity at the *Npts* different points.

    """

    # function display 
    print '---- DAEPy::getSkinFrictionBetweenPosition ----'

    Cf = np.zeros(Npts)
    Utau = np.zeros(Npts)
    Yplus = np.zeros(Npts)
    Tw = np.zeros(Npts)
    pos = np.zeros((Npts, 3))

    for i in range(Npts):
        
        x_p_temp = x_p0 + (x_p1-x_p0)/(Npts-1)*i

        [pos[i,:], Cf[i], Utau[i], Yplus[i], Tw[i]] = getSkinFrictionAtPosition(data_outVTK, wall_outVTK, cord_choice, x_p_temp, Uinf, Rhoinf)

    return [pos, Cf, Utau, Yplus, Tw]


def getDeltaBetweenPosition(data_outVTK, wall_outVTK, cord_choice, x_p0, x_p1, Npts, Uinf=None, Rhoinf=None):

    """
    Compute the boundary layer thickness for *Npts* equally distributed between the 2 position
    defined thanks to *x_p0*, *x_p1* and *cord_choice*. See the documentation of the function
    getDeltaAtPosition() for more information.

    Parameters
    ----------
    data_outVTK : VTK output object from a VTK reader
        This field data MUST fulfill 2 conditions:
            - contain the field "U_AVG",
            - represent a 2D field data, as we want to extract a 1D velocity
            profile from it.

    wall_outVTK : VTK output object from a VTK reader
        This field is supposed to describe a curve with 1D cell type, like in getLineTangentialVector()

    cord_choice : integer, 0, 1 or 2
        Gives the axis that is going to be used to define the point where the velocity
        profile is extracted. Convention :
            - 0 = x axis
            - 1 = y axis
            - 2 = z axis

    x_p0, x_p1 : float between 0 and 1
        gives the bound of the line where to compute the BL thickness.

    Npts: integer
        Number of points where to compute delta wanted.

    Uinf : float 
        (optional) Free stream velocity

    Rhoinf : float 
        (optional) Free stream density

    Returns
    -------
    pos : vector of tuple(3)
        the coordinates of the points where delta have been computed.

    delta : vector of floats
        The BL thickness at the *Npts* different points.

    deltaS : vector of floats
        The compressible BL displacement thickness at the *Npts* different points.

    theta : vector of floats
        The compressible BL momentum thickness at the *Npts* different points.

    """

    # function display 
    print '---- DAEPy::getDeltaBetweenPosition ----'

    delta = np.zeros(Npts)
    deltaS = np.zeros(Npts)
    theta = np.zeros(Npts)
    pos = np.zeros((Npts, 3))

    for i in range(Npts):
        
        x_p_temp = x_p0 + (x_p1-x_p0)/(Npts-1)*i

        [pos[i,:], delta[i], deltaS[i], theta[i]] = getDeltaAtPosition(data_outVTK, wall_outVTK, cord_choice, x_p_temp, Uinf, Rhoinf)

    return [pos, delta, deltaS, theta]

def getIncDeltaBetweenPosition(data_outVTK, wall_outVTK, cord_choice, x_p0, x_p1, Npts, Uinf=None, Rhoinf=None):

    """
    Compute the boundary layer thickness for *Npts* equally distributed between the 2 position
    defined thanks to *x_p0*, *x_p1* and *cord_choice*. See the documentation of the function
    getDeltaAtPosition() for more information.

    Parameters
    ----------
    data_outVTK : VTK output object from a VTK reader
        This field data MUST fulfill 2 conditions:
            - contain the field "U_AVG",
            - represent a 2D field data, as we want to extract a 1D velocity
            profile from it.

    wall_outVTK : VTK output object from a VTK reader
        This field is supposed to describe a curve with 1D cell type, like in getLineTangentialVector()

    cord_choice : integer, 0, 1 or 2
        Gives the axis that is going to be used to define the point where the velocity
        profile is extracted. Convention :
            - 0 = x axis
            - 1 = y axis
            - 2 = z axis

    x_p0, x_p1 : float between 0 and 1
        gives the bound of the line where to compute the BL thickness.

    Npts: integer
        Number of points where to compute delta wanted.

    Uinf : float 
        (optional) Free stream velocity

    Rhoinf : float 
        (optional) Free stream density

    Returns
    -------
    pos : vector of tuple(3)
        the coordinates of the points where delta have been computed.

    delta : vector of floats
        The BL thickness at the *Npts* different points.

    deltaS : vector of floats
        The compressible BL displacement thickness at the *Npts* different points.

    theta : vector of floats
        The compressible BL momentum thickness at the *Npts* different points.

    """

    # function display 
    print '---- DAEPy::getDeltaBetweenPosition ----'

    delta = np.zeros(Npts)
    deltaS = np.zeros(Npts)
    theta = np.zeros(Npts)
    pos = np.zeros((Npts, 3))

    for i in range(Npts):
        
        x_p_temp = x_p0 + (x_p1-x_p0)/(Npts-1)*i

        [pos[i,:], delta[i], deltaS[i], theta[i]] = getIncDeltaAtPosition(data_outVTK, wall_outVTK, cord_choice, x_p_temp, Uinf, Rhoinf)

    return [pos, delta, deltaS, theta]


def getVelocityAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf=None, Rhoinf=None):

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

    y_coo : vector of floats
        wall normal coordinates at x_perc.

    yp : vector of floats
        dimensionless wall normal coordinates at x_perc.

    Ut : vector of floats
        tangential velocity along y_coo coordinates.

    Up : vector of floats
        dimensionless tangential velocity along y_coo coordinates.

    """

    # function display 
    print '---- DAEPy::getVelocityAtPosition ----'

    # test if the field RHO_AVG is present
    if data_outVTK.GetPointData().HasArray('RHO_AVG')!=1:
        raise ValueError("Error : field RHO_AVG not present")

    # test if the field U_AVG is present
    if data_outVTK.GetPointData().HasArray('U_AVG')!=1:
        raise ValueError("Error : field U_AVG not present")

    # test if the field MU_LAM_AVG is present
    if data_outVTK.GetPointData().HasArray('MU_LAM_AVG')!=1:
        raise ValueError("Error : field MU_LAM_AVG not present")

    # get the vector used to slice the field data
    [ori, vec_tan] = getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    data1D = getSlice(data_outVTK, ori, vec_tan)

    # extract the tangential velocity
    [Vcoords, Rho, U, Mu] = getArrayFromPointData(data1D, ['RHO_AVG','U_AVG','MU_LAM_AVG'])
    Ut = np.sum(U*vec_tan, axis=1)


    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()
    
    # sort cordinates along the line orthogonal to the wall
    y_coo = y_coo[id_sort]

    # sort the other quantities according to the coordinates y_coo. Velocity at the wall Ut[0] is
    # set to 0 to respect no slip condition
    Mu = Mu[id_sort]
    Rho = Rho[id_sort]
    Ut = Ut[id_sort] ; Ut[0] = 0

    [_, Cf, Utau, _, _] = getSkinFrictionAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf, Rhoinf)

    Up = Ut/Utau
    yp = y_coo*Utau*Rho[0]/Mu[0]


    # function output
    print ''
    return [ori, y_coo, yp, Ut, Up]


def getVanDriestIIAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf=None, Rhoinf=None, Muinf=None, Tinf=None):

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

    Muinf : float 
        (optional) Free stream viscosity

    Tinf : float 
        (optional) Free stream temperature

    Returns
    -------
    ori : tuple(3)
        The point coordinates of the point at the wall where we extract the data.

    Cfi : float
        incompressible skin friction coefficient (using Van Driest II transformation) at position.

    ReOi : float
        incompressible Reynolds theta (using Van Driest II transformation) at position.

    ReT : float
        Reynolds tau based on shear stress velocity at position.

    """

    # function display 
    print '---- DAEPy::getVanDriestIIAtPosition ----'

    # test if the field T_AVG is present
    if data_outVTK.GetPointData().HasArray('T_AVG')!=1:
        raise ValueError("Error : field T_AVG not present")

    # test if the field MU_LAM_AVG is present
    if data_outVTK.GetPointData().HasArray('MU_LAM_AVG')!=1:
        raise ValueError("Error : field MU_LAM_AVG not present")

    # get the vector used to slice the field data
    [ori, vec_tan] = getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    data1D = getSlice(data_outVTK, ori, vec_tan)

    # extract the tangential velocity
    [Vcoords, T, Mu] = getArrayFromPointData(data1D, ['T_AVG','MU_LAM_AVG'])

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()
    
    # sort cordinates along the line orthogonal to the wall
    y_coo = y_coo[id_sort]

    # sort the other quantities according to the coordinates y_coo. Velocity at the wall Ut[0] is
    # set to 0 to respect no slip condition
    T = T[id_sort]
    Mu = Mu[id_sort]

    # Automatic definition of the external (free stream) viscosity based on the extracted viscosity if Uinf = None
    if Muinf == None:
        Muinf = np.mean(Mu[-10:-1])
        # Muinf = Mu[-1]

    # Automatic definition of the external (free stream) temperature based on the extracted temperature if Tinf = None
    if Tinf == None:
        Tinf = np.mean(T[-10:-1])
        # Tinf = T[-1]


    [ori, ReD, ReDS, ReO, ReT] = getReynoldsAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf, Rhoinf, Muinf)

    [_, Cf, Utau, _, _] = getSkinFrictionAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf, Rhoinf)


    # Compute the Van Driest II incompressible Reynolds theta
    ReOi = ReO*(Muinf/Mu[0])

    # Compute the Van Driest II incompressible skin friction parameters
    alpha = (T[0]/Tinf - 1)/(np.sqrt(T[0]/Tinf*(T[0]/Tinf - 1)))
    Fc = (T[0]/Tinf - 1)/(np.arcsin(alpha)**2)

    # Compute the Van Driest II incompressible skin friction coefficient
    Cfi = Fc*Cf

    # function output
    print ''
    return [ori, Cfi, ReOi, ReT]


def getVanDriestIIBetweenPosition(data_outVTK, wall_outVTK, cord_choice, x_p0, x_p1, Npts, Uinf=None, Rhoinf=None, Muinf=None, Tinf=None):

    """
    Compute the Van Driest II transformation of the Skin friction coefficient and Reynolds theta
    for *Npts* equally distributed between the 2 position defined thanks to *x_p0*, *x_p1* and *cord_choice*.
    See the documentation of the function getVanDriestIIAtPosition() for more information.

    Parameters
    ----------
    data_outVTK : VTK output object from a VTK reader
        This field data MUST fulfill 2 conditions:
            - contain the field "U_AVG",
            - represent a 2D field data, as we want to extract a 1D velocity
            profile from it.

    wall_outVTK : VTK output object from a VTK reader
        This field is supposed to describe a curve with 1D cell type, like in getLineTangentialVector()

    cord_choice : integer, 0, 1 or 2
        Gives the axis that is going to be used to define the point where the velocity
        profile is extracted. Convention :
            - 0 = x axis
            - 1 = y axis
            - 2 = z axis

    x_p0, x_p1 : float between 0 and 1
        gives the bound of the line where to compute the BL thickness.

    Npts: integer
        Number of points where to compute delta wanted.

    Uinf : float 
        (optional) Free stream velocity

    Rhoinf : float 
        (optional) Free stream density

    Muinf : float 
        (optional) Free stream viscosity

    Tinf : float 
        (optional) Free stream temperature

    Returns
    -------
    pos : vector of tuple(3)
        the coordinates of the points where delta have been computed.

    Cfi : vector of floats
        incompressible skin friction coefficient (using Van Driest II transformation) at the *Npts* different points.

    ReOi : vector of floats
        incompressible Reynolds theta (using Van Driest II transformation) at the *Npts* different points.

    ReT : vector of floats
        Reynolds tau based on shear stress velocity.

    """

    # function display 
    print '---- DAEPy::getVanDriestIIBetweenPosition ----'

    Cfi = np.zeros(Npts)
    ReOi = np.zeros(Npts)
    ReT = np.zeros(Npts)
    pos = np.zeros((Npts, 3))

    for i in range(Npts):
        
        x_p_temp = x_p0 + (x_p1-x_p0)/(Npts-1)*i

        [pos[i,:], Cfi[i], ReOi[i], ReT[i]] = getVanDriestIIAtPosition(data_outVTK, wall_outVTK, cord_choice, x_p_temp, Uinf, Rhoinf, Muinf, Tinf)

    return [pos, Cfi, ReOi, ReT]


def getVanDriestVelocityAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf=None, Rhoinf=None, Muinf=None, Tinf=None):

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

    Muinf : float 
        (optional) Free stream viscosity

    Tinf : float 
        (optional) Free stream temperature

    Returns
    -------
    ori : tuple(3)
        The point coordinates of the point at the wall where we extract the data.

    yp : vector of floats
        dimensionless wall normal coordinates at x_perc.

    Uvd : vector of floats
        dimensionless Van Driest velocity along yp coordinates.

    """

    # function display 
    print '---- DAEPy::getVanDriestVelocityAtPosition ----'

    # test if the field T_AVG is present
    if data_outVTK.GetPointData().HasArray('T_AVG')!=1:
        raise ValueError("Error : field T_AVG not present")

    # get the vector used to slice the field data
    [ori, vec_tan] = getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    data1D = getSlice(data_outVTK, ori, vec_tan)

    # extract the tangential velocity
    [Vcoords, T] = getArrayFromPointData(data1D, ['T_AVG'])

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()
    
    # sort cordinates along the line orthogonal to the wall
    y_coo = y_coo[id_sort]

    # sort the other quantities according to the coordinates y_coo. Velocity at the wall Ut[0] is
    # set to 0 to respect no slip condition
    T = T[id_sort]


    [ori, y_coo, yp, Ut, Up] = getVelocityAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf, Rhoinf)

    [_, Cf, Utau, _, _] = getSkinFrictionAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf, Rhoinf)

    # Van Driest I Transformation parameters
    cp = 1005.25
    Prt = 0.9
    C1 = 1.0
    Gamma = 1.4 

    # Adiabatic case H = 0
    H = 0

    # Compute Van Driest I coefficients
    Mt = Utau/(np.sqrt((1.4-1.0)*cp*T[0]))
    R = Mt*np.sqrt((Gamma-1.)*Prt/2)
    D = np.sqrt(C1 + (R**2)*(H**2))

    # Compute the Van Driest I velocity
    Uvd = 1/R*(np.arcsin(R*(Up+H)/D)-np.arcsin(R*H/D))

    # function output
    print ''
    return [ori, yp, Uvd]


def getRMSAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf=None, Rhoinf=None):

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
    [ori, vec_tan] = getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    data1D = getSlice(data_outVTK, ori, vec_tan)

    # extract the tangential velocity
    [Vcoords, Urms, Urey, Rho, Mu] = getArrayFromPointData(data1D, ['U_RMS','U_REY','RHO_AVG','MU_LAM_AVG'])

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()
    
    # sort cordinates along the line orthogonal to the wall
    y_coo = y_coo[id_sort]
    Rho = Rho[id_sort]
    Mu = Mu[id_sort]

    [_, Cf, Utau, _, _] = getSkinFrictionAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf, Rhoinf)

    # sort the RMS and REY quantities according to the coordinates y_coo. Normalisation using Utau is also performed
    Urms = Urms[id_sort]**2/Utau**2
    Urey = np.abs(Urey[id_sort])/Utau**2

    yp = y_coo*Utau*Rho[0]/Mu[0]

    # Compute scaling factor for comparison with incompressible data
    scaling = np.sqrt(Rho/Rho[0])

    # function output
    print ''
    return [ori, yp, Urms.T, Urey.T, scaling]

def getRMSBetweenPosition(data_outVTK, wall_outVTK, cord_choice, x_p0, x_p1, Npts, Uinf=None, Rhoinf=None, Muinf=None, Tinf=None):

    """

    """

    # function display 
    print '---- DAEPy::getRMSBetweenPosition ----'

    Cfi = np.zeros(Npts)
    ReOi = np.zeros(Npts)
    ReT = np.zeros(Npts)
    pos = np.zeros((Npts, 3))
    Urms = np.zeros((Npts, 3))
    Urey = np.zeros((Npts, 3))

    for i in range(Npts):
        
        x_p_temp = x_p0 + (x_p1-x_p0)/(Npts-1)*i

        [pos[i,:], Cfi[i], ReOi[i], ReT[i]] = getVanDriestIIAtPosition(data_outVTK, wall_outVTK, cord_choice, x_p_temp, Uinf, Rhoinf, Muinf, Tinf)

    return [pos, Cfi, ReOi, ReT]


def getInnerScalingAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf=None, Rhoinf=None):

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
    [ori, vec_tan] = getLineTangentialVector(wall_outVTK, x_perc, cord_choice)

    # slice the 2D field data using the tangential vector to the wall
    data1D = getSlice(data_outVTK, ori, vec_tan)

    # extract the tangential velocity
    [Vcoords, Rho, Prms, Rhorms, Rhoavg, Tavg, Trms, T, Mu, U] = getArrayFromPointData(data1D, ['RHO','P_RMS','RHO_RMS','RHO_AVG','T_AVG','T_RMS','T','MU_LAM_AVG','U_AVG'])
    Ut = np.sum(U*vec_tan, axis=1)

    # define a new scalar coordinates along the line orthogonal to the wall
    y_coo = np.array(getScalarCoord(Vcoords, 1))
    id_sort = y_coo.argsort()
    
    # sort cordinates along the line orthogonal to the wall
    y_coo = y_coo[id_sort]
    Rho = Rho[id_sort]
    Rhorms = Rhorms[id_sort]
    Rhoavg = Rhoavg[id_sort]
    Prms = Prms[id_sort]
    Tavg = Tavg[id_sort]
    Trms = Trms[id_sort]
    T = T[id_sort]
    Mu = Mu[id_sort]
    Ut = Ut[id_sort]

    # Automatic definition of the external (free stream) velocity based on the extracted velocity if Uext = None
    if Uinf == None:
        # Uinf = Ut[-1]
        Uinf = np.mean(Ut[-35:-1])

    # Automatic definition of the external (free stream) density based on the extracted density if RhoExt = None
    if Rhoinf == None:
        Rhoinf = Rhoavg[-1]
        # Rhoinf = np.mean(Rho[-35:-1])

    [_, Cf, Utau, _, _] = getSkinFrictionAtPosition(data_outVTK, wall_outVTK, cord_choice, x_perc, Uinf, Rhoinf)

    scaling = np.sqrt(Rhoavg/Rhoavg[0])

    Tw = (Rho*T-Rhorms*Trms)/Rhoavg
    Mtau = Utau/np.sqrt((1.4*287.2142857*Tw[0]))

    prmsp = Prms/(Cf*(0.5*Rhoinf*Uinf**2))
    rhormsp = Rhorms/(1.4*Rhoavg[0]*Mtau**2)

    yp = y_coo*Utau*Rhoavg[0]/Mu[0]

    # function output
    print ''
    return [ori, yp, Mtau, prmsp, rhormsp, scaling]


def fctTemp():

    """
    Descritpion...
    
    .. math:: X(e^{j\omega } ) = x(n)e^{ - j\omega n}

    The value of :math:`\omega` and :math:`\\theta` is larger than 5.

    .. math::

     x(n) * y(n) \Leftrightarrow X(e^{j\omega } )Y(e^{j\omega } )\\
     another equation here
    
    

    Parameters
    ----------
    
    Returns
    -------
    

    """

    # function display 
    print '---- DAEPy::fctTemp ----'
    print '--> Some comments...'


    # function output
    print ''
    return []