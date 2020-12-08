
# -*- coding: utf-8 -*-
"""
dataManagment.py
"""

import vtk
import os, sys
import toolCase as util

import numpy as np

from vtk.util.numpy_support import vtk_to_numpy

sys.dont_write_bytecode = True

"#################### FIELDS FROM .vtk  ####################"

def getFields(out_f,fields_names):
    
    fields_dictionary = {}
    
    fields_number = out_f.GetPointData().GetNumberOfArrays()
    fields_avails = [out_f.GetPointData().GetArrayName(i) for i in range(fields_number)]

#    [GridCoords, field_array] = getArrayFromPointData(out_f, fields_avails[0]) # just to pick the Grid, anyway

    "###############  GET THE FIELDS  ###############"

    if not fields_names: fields_names = fields_avails ; print("Fields list is empty: looking for available fields")

    for field_name in fields_names: 

        if field_name not in fields_avails: continue
        
        [GridCoords, field_array] = getArrayFromPointData(out_f, field_name)
        
        if field_array.ndim > 1:
        
            field_array_comp = field_array[:,0]
            fields_dictionary.update({field_name + '_t': field_array_comp})
            field_array_comp = field_array[:,1]
            fields_dictionary.update({field_name + '_r': field_array_comp})
            
            if np.shape(field_array)[1] == 3:
            
                field_array_comp = field_array[:,2]
                fields_dictionary.update({field_name + '_z': field_array_comp})
            
            field_array_magn = np.linalg.norm(field_array,axis=1)
            fields_dictionary.update({field_name + '_m': field_array_magn})
        
        else:
            fields_dictionary.update({field_name : field_array})

    return GridCoords, fields_dictionary

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
    print ('\n---- DAEPy::getOutputVTKwithPointDataFromFile ----')
    
    # test if the file exists
    print ('\n--> Reading', fileName)
    if not os.path.isfile(fileName):
        raise ValueError("Error : file does not exists")

    extension = os.path.splitext(fileName)[-1]
    if extension == '.vtu':
        reader = vtk.vtkUnstructuredGridReader()
    elif extension == '.vtk':
        reader = vtk.vtkStructuredGridReader()
    elif extension == '.pvtu':
        reader = vtk.vtkXMLPUnstructuredGridReader()
    elif extension == '.vtp':
        reader = vtk.vtkXMLPolyDataReader()
    elif extension == '.vtm':
        # TODO : To check
        reader = vtk.vtkXMLMultiBlockDataReader()
        reader = vtk.MergeBlocks(reader)
    else:
        raise ValueError("Error: unknown extension of file " + fileName)

#    reader.Initialize()
    reader.SetFileName(fileName)

    reader.ReadAllScalarsOn()
#    reader.ReadAllVectorsOn()

    reader.Update()
    data_outVTK = reader.GetOutput()

#    # All the data are transfered to the nodes
#    c2p = vtk.vtkCellDataToPointData()
#    c2p.SetInputData(data_outVTK)
#    c2p.Update()
#    data_outVTK = c2p.GetOutput()

    # list the fields available

    n_fields = data_outVTK.GetPointData().GetNumberOfArrays()

    print ('\n--> Available:', n_fields, 'fields\n')

    for i in range(n_fields):
        print ('        -', data_outVTK.GetPointData().GetArrayName(i))


    print ('')
    return data_outVTK


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
    print ('\n---- DAEPy::getArrayFromPointData ----')

    coord = np.array(
        [data_outVTK.GetPoint(i) for i in range(data_outVTK.GetNumberOfPoints())],
        )
#        dtype=np.float32)

    if isinstance(field_name, list):
        for f in field_name:
            print ('\n--> extract fields', f)
            data_arr = [vtk_to_numpy(data_outVTK.GetPointData().GetArray(f)) for f in field_name]
    else:
        print ('\n--> extract field', field_name)
        data_arr = vtk_to_numpy(data_outVTK.GetPointData().GetArray(field_name))

#    print '\n--> extract fields', [f for f in field_name]
#    data_arr = [vtk_to_numpy(data_outVTK.GetPointData().GetArray(f)) for f in field_name]


    print ('')
#    return [coord] + data_arr
    return [coord , data_arr]


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
    print ('---- DAEPy::getSlice ----')

    # Plane of the Slice
    plane = vtk.vtkPlane()
    plane.SetOrigin(orig)
    plane.SetNormal(nor)
    print ('--> normal used to slice: ', nor)

    # Cutter plane
    planeCut = vtk.vtkCutter()
#    planeCut.SetInputData(data_outVTK) # VTK6 SetInput() repl. with SetInputData() vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
    planeCut.SetInput(data_outVTK)    # sometimes to switch one another
    planeCut.SetCutFunction(plane)

    # Make the slice
    planeCut.Update()
    slice_outVTK = planeCut.GetOutput()

    print ('')
    return slice_outVTK


"#################### FIELDS FROM Slice.vtk  ####################"

def getFieldsFromSlice(field_slice,fields_names):

    out_f = field_slice ; fields_dictionary = {}
    
    fields_number = out_f.GetPointData().GetNumberOfArrays()
    fields_avails = [out_f.GetPointData().GetArrayName(i) for i in range(fields_number)]

    [GridCoords, field_array] = getArrayFromPointData(out_f, fields_avails[0]) # just to pick the Grid, anyway

    "###############  GET THE FIELDS  ###############"

    if not fields_names: fields_names = fields_avails ; print("Fields list is empty: looking for available fields")

    for field_name in fields_names: 

        if field_name not in fields_avails: continue

        [GridCoords, field_array] = getArrayFromPointData(out_f, field_name)

        if field_array.ndim > 1:

            field_array_comp = field_array[:,0]
            fields_dictionary.update({field_name + '_t': field_array_comp})
            field_array_comp = field_array[:,1]
            fields_dictionary.update({field_name + '_r': field_array_comp})

            if np.shape(field_array)[1] == 3:
		
                field_array_comp = field_array[:,2]
                fields_dictionary.update({field_name + '_z': field_array_comp})

            field_array_magn = np.linalg.norm(field_array,axis=1)
            fields_dictionary.update({field_name + '_m': field_array_magn})

        else:
            fields_dictionary.update({field_name : field_array})

    return GridCoords, fields_dictionary

