
# -*- coding: utf-8 -*-
"""
plotFields.py
"""

import os, sys, re, copy #, shutil
import dataManagement as rosie
import toolCase as util
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm     as cm
import matplotlib.ticker as ticker


from scipy.interpolate import griddata #, interpn, interp2d
from matplotlib import colors

sys.dont_write_bytecode = True

"#################### INSTANTANEOUS CONTOUR PLOTS #########################" 

def snapshots(run_dir,run_ver,res_dir,fn_head,name_list):

    nan_list = ['Crit_rad']

    res_dir = res_dir + '2D_snapshot/'

    if util.chk_dir(res_dir) == False: os.makedirs(res_dir) 

    zones_dictionary = {} ; zones_fields = {}
    
    orig = [0,0,0] ; norm = [1,0,0]

    for root, dirs, files in os.walk("%s"%run_dir):

        field_numb = - 69
    
        for filename in files:

            if re.match('%s%s*.*.vtk'%(run_ver,fn_head), filename):

                print ('\n--> looking for %s'%filename) 
    
                filename = filename.lstrip('%s%s'%(run_ver,fn_head))
                filename = filename.rstrip('.vtk')

                if int(filename) > field_numb:

                    field_numb = int(filename)
                    last_field = filename

        break   #prevent decending into subfolders
    
    print ('\n--> last field found: %s'%last_field) 

    for root, dirs, files in os.walk("%s"%run_dir):
    
        for filename in files:

            if re.match('%s*.*%s.vtk'%(run_ver,last_field), filename):
    
                zonename = filename.lstrip('%s'%run_ver)
                zonename = filename.rstrip('%s.vtk'%last_field)

                print ('\n--> zone found: %s'%zonename) 

                zones_dictionary.update({zonename : filename})

        break   #prevent decending into subfolders

    for zone in zones_dictionary:   

        file_d = '%s'%run_dir + zones_dictionary[zone] ; out_f = rosie.getOutputVTKwithPointDataFromFile(file_d)
        
        inst_slice = rosie.getSlice(out_f, orig, norm) ; print ('\n--> Slicing %s'%file_d) 
        
        grid, inst_fields = rosie.getFieldsFromSlice(inst_slice,name_list)
        
        zones_fields.update({zone : [grid, inst_fields]})

    fields_list = []

    for zone in zones_fields:   

        [grid, inst_fields] = zones_fields[zone]

        fields_list = list(set().union(fields_list, list(inst_fields.keys())))
    
    nptsx, nptsy = 1000, 2000 ; zmin, zmax = 0.0, 0.0

    for field in fields_list:   

        print ('\n--> plotting field %s'%field)

        fig = plt.figure(figsize=(5, 8)) ; ax = fig.add_subplot(111)

        ax.set_aspect(1, adjustable = 'box') ; cmap = cm.get_cmap('Blues_r')

        for zone in ['%sJE'%run_ver,'%sPI'%run_ver]:   # zones_fields: # not to match the zones 

            if   zone == '%sJE'%run_ver: nptsx, nptsy = 1000, 2000 #125 #2000
            elif zone == '%sPI'%run_ver: nptsx, nptsy =   50,  100 #500 #1000
            
            [grid, inst_fields] = zones_fields[zone]
 
            x, y = grid[:,1], grid[:,2]
    
            xmin, xmax = min(x), max(x) ; ymin, ymax = min(y), max(y)
            
            xi = np.linspace(xmin, xmax, nptsx) ; yi = np.linspace(ymin, ymax, nptsy)
            
            xi, yi = np.meshgrid(xi, yi)#, indexing='ij')#, sparse=True)
            
            ax.set_xlim(-25.0, 25.0) ; ax.set_ylim( -2*np.pi, 100.0)
#            ax.set_xlim(xmin, xmax)  ; ax.set_ylim(ymin, ymax)

            func = lambda x, pos: "" if np.isclose(x,0) else "%.f" % x
            
            plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(func))
            plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(func))

            if field in inst_fields :
            
                z = inst_fields[field] ; zmin, zmax = min(z), max(z)
                
                lmin = min(z) ; lmax = max(z)
            
                if   field in nan_list:
            
                    if   field == 'Crit_rad' : z[z <= 1.e-8 ] = np.nan
                    else                     : z[z <= 1.e-16] = np.nan
                    
                    lmin = np.nanmin(z) ; lmax = np.nanmax(z)
            
                zi = griddata((x, y), z, (xi, yi), method='linear')
            
                if   field in nan_list: zi = np.nan_to_num(zi)#, nan=-9999) # need upgrade numpy 
                
                if zone == '%sJE'%run_ver : zpipe = util.ProbeAtLocation(zi, xi, yi, 0.0, 0.1)

            else:

                zi = griddata((x, y), inst_fields['Uf_trz_m'], (xi, yi), method='linear') ; zi.fill(zpipe)
            
            levels_n  = 8 ; cm_format = None ; my_cmap = copy.copy(cmap) # to init default values 
                                
            if field in nan_list: my_cmap.set_under(color='lightyellow')
            
            if   field == 'Temperature' : lmin =  1.05    ; lmax =   1.35   ; levels_n  = 6
            elif field == 'Rho'         : lmin =  0.65    ; lmax =   0.95   ; levels_n  = 6
            elif field == 'Sat_ratio'   : lmin =  0.0     ; lmax = 900.0    ; levels_n  = 9
            elif field == 'Uf_trz_t'    : lmin = -1.0     ; lmax =   1.0    ; levels_n  = 10
            elif field == 'Uf_trz_r'    : lmin = -0.5     ; lmax =   0.5    ; levels_n  = 10
            elif field == 'Uf_trz_z'    : lmin = -0.3     ; lmax =   2.4    ; levels_n  = 18
            elif field == 'Uf_trz_m'    : lmin =  0.0     ; lmax =   2.4    ; levels_n  = 16
            elif field == 'Y_vapour'    : lmin =  1.0e-9  ; lmax =   3.2e-3 ; levels_n  = 16 ; oo_magn = -3
            elif field == 'Y_saturat'   : lmin =  1.0e-9  ; lmax =   3.2e-6 ; levels_n  = 16 ; oo_magn = -6
            
            elif re.match('fp_.+',field) and zone == '%sJE'%run_ver :

                lmin = util.OOMRoundUp(lmin) ; lmax = util.OOMRoundUp(lmax)

                if   field == 'fp_vap' or re.match('fp_trz_.+',field) : lmin = lmin/1.e+4 ; lmax = lmax/1.e+4
                else                                                  : lmin = lmin/1.e+1 ; lmax = lmax/1.e+1
                
                labs = min(abs(lmin),abs(lmax)) if min(abs(lmin),abs(lmax)) > 1.e-16 else max(abs(lmin),abs(lmax))
                
                lmin = - labs ; lmax = labs ; levels_n  = 5
            
            if   zone == '%sJE'%run_ver and np.isclose(max(inst_fields['Y_vapour']),1.8e-3, rtol=1e-4, atol=1e-4):
            
                if   field == 'J_rate'      : lmin =  0.0     ; lmax =  2.4    ; levels_n  = 9
                elif field == 'Rate_nucl'   : lmin =  1.0e-16 ; lmax =  1.0e-4
                elif field == 'Crit_rad'    : lmin =  1.0e-7  ; lmax =  9.0e-7                 ; oo_magn = -7
            
            elif zone == '%sJE'%run_ver and np.isclose(max(inst_fields['Y_vapour']),2.1e-3, rtol=1e-4, atol=1e-4):
            
                if   field == 'J_rate'      : lmin =  0.0     ; lmax =  3.6    ; levels_n  = 9
                elif field == 'Rate_nucl'   : lmin =  1.0e-16 ; lmax =  1.0e-4
                elif field == 'Crit_rad'    : lmin =  1.0e-7  ; lmax =  9.0e-7                 ; oo_magn = -7
            
            elif zone == '%sJE'%run_ver and np.isclose(max(inst_fields['Y_vapour']),2.5e-3, rtol=1e-4, atol=1e-4):
            
                if   field == 'J_rate'      : lmin =  0.0     ; lmax =  7.2    ; levels_n  = 9
                elif field == 'Rate_nucl'   : lmin =  1.0e-16 ; lmax =  1.0e-4
                elif field == 'Crit_rad'    : lmin =  1.0e-7  ; lmax =  9.0e-7                 ; oo_magn = -7
            
            lcountour = np.linspace(lmin, lmax, levels_n + 1)
            
            if   field == 'Sat_ratio'   : lcountour[0] = 1.0 ; lmin = 1.0
            elif field == 'J_rate'      : lcountour[0] = 1.0e-3
            
            elif field in ['Y_vapour','Y_saturat','Crit_rad']:
            
                    cm_format = util.OOMFormatter(oo_magn, mathText=False)
            
            if field not in name_list: plt.title(field)
            
            # CONTOUR: draws the boundaries of the isosurfaces & fill the contour plots
            
            cp = plt.contourf(xi, yi, zi, levels=lcountour , vmin=lmin, vmax=lmax, extend='both',cmap=my_cmap)

        clb = plt.colorbar(cp,format=cm_format) ; fig.tight_layout() #; clb.ax.set_title('NaN', pad=5.0) # upgrade plt
    
        plt.savefig(res_dir + '/%s.png'%field, format='png', dpi=900) ; plt.close('all')

    return []

