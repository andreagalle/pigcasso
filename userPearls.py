
# -*- coding: utf-8 -*-
"""
userPearls.py
"""

import os, sys, re, copy, csv #, math, shutil
import dataManagement as rosie
import toolCase as util
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm     as cm
import matplotlib.ticker as ticker

import matplotlib
matplotlib.rcParams['backend'] = 'Agg' # configure backend here

from scipy.interpolate import griddata #, interpn, interp2d
from matplotlib import colors
from itertools import combinations 

sys.dont_write_bytecode = True

"###############################################################"

"""

Put your code here below and it will be run as a custom script 
in addition to what is run by default every time pigcasso runs.

"""

def throw(run_directory,run_version,res_directory):

    check = False

    statistics(run_directory,run_version,res_directory) ; check = True 

    return check

"#############################################################"

def statistics(run_directory,run_version,res_directory):

    fields_list = []

    contour_plots     (run_directory, "flu_mean",             run_version, res_directory, fields_list)
    contour_plots     (run_directory, "prt_mean",             run_version, res_directory, fields_list)
    profile_plots     (run_directory, "fcl_mean",             run_version, res_directory, fields_list)
    profile_plots     (run_directory, "fra_mean",             run_version, res_directory, fields_list)
    profile_plots     (run_directory, "prad_pdf",             run_version, res_directory, fields_list)
    profile_plots     (run_directory, "pvap_pdf",             run_version, res_directory, fields_list)
    profile_plots     (run_directory, "pexp_pdf",             run_version, res_directory, fields_list)
    profile_plots     (run_directory, "pvol_pdf",             run_version, res_directory, fields_list)
    scatter_plots     (run_directory,                         run_version, res_directory, fields_list)
   
    norm_profile_plots(run_directory, "fra_mean", "fcl_mean", run_version, res_directory)

#     run_directory = "cfr-pdf-rad/"
#
#    cfr_profile_plots (run_directory, "prad_pdf",             run_version, res_directory)
#    cfr_profile_plots (run_directory, "pexp_pdf",             run_version, res_directory)
#    cfr_profile_plots (run_directory, "pvol_pdf",             run_version, res_directory)

    run_directory = "../../doc/cases/"

    cfr_DNSvsExp     (run_directory, "fried2000F10.11",      run_version, res_directory)

#    mdot_HKvsMA      (run_directory,                         run_version, res_directory)
    Jrate_sigma      (run_directory,                         run_version, res_directory)

"#################### PHYSICAL CONSTANTS #############################################################"

Sh = 2. ; Re = 1500. ; Sc = 1. ; Ma = 0.046 ; alpha = 1. ; theta_ref = 299. ; rho_liq = 910.586 ; rho_ref = 1.1531

Wvap = 278.34e-3 ; Wgas = 28.29e-3 ; Lat_heat = 91.70e+3 ; theta_boi = 613. ; Runi = 8.314 ; Rvap = 0.102 

Navg = 6.022e+23 ; Kbol = 1.381e-23 ; Mmol = Wvap/Navg ; theta_zero = 273.15 ; L_ref = 1.75e-3 ; U_ref = 16.212 

"#################### MEAN CONTOUR PLOTS ###############################################################################" 

def contour_plots(run_dir,run_out,run_ver,res_dir,name_list):

    nan_list = ['Crit_Radius','Crit_Radius_rms','part_number','part_radius','mass_loading','dexp_number']

    res_dir = res_dir + '2D_contour/'

    if util.chk_dir(res_dir) == False: os.makedirs(res_dir) 

    file_d = run_dir + '%s%s.vtk'%(run_ver,run_out) ; out_f = rosie.getOutputVTKwithPointDataFromFile(file_d)

    grid, mean_fields = rosie.getFields(out_f,name_list)

    x, y = grid[:,1], grid[:,2] ; npts = 500
    
    xmin, xmax = min(x), max(x) ; ymin, ymax = min(y), max(y)

    xi = np.linspace(xmin, xmax, npts) ; yi = np.linspace(ymin, ymax, npts)
    
    xi, yi = np.meshgrid(xi, yi)#, sparse=True)

    for field in mean_fields:   

        print '\n--> plotting field %s'%field

        fig = plt.figure(figsize=(9, 3)) ; ax = fig.add_subplot(111)

        ax.set_aspect(1, adjustable = 'box') ; cmap = cm.get_cmap('Blues_r')

        ax.set_xlim(0.0, 99.0) ; ax.set_ylim(0.0, 24.0)

        func = lambda x, pos: "" if np.isclose(x,0) else "%.f" % x

        plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(func))
        plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(func))

        z = mean_fields[field] ; lmin = min(z) ; lmax = max(z)

        if field in nan_list: 

           if   field == 'part_number' : z[z <= 1.e-1 ] = np.nan
           elif field == 'dexp_number' : z[z <= 1.e-1 ] = np.nan
           elif field == 'part_radius' : z[z <= 1.e-7 ] = np.nan
           else                        : z[z <= 1.e-16] = np.nan

           lmin = np.nanmin(z) ; lmax = np.nanmax(z)

	zi = griddata((x, y), z, (xi, yi), method='linear')

        levels_n  = 8 ; cm_format = None ; my_cmap = copy.copy(cmap) # to init default values 

        if field in nan_list: my_cmap.set_under(color='lightyellow') ; zi = np.nan_to_num(zi)#, nan=-9999) # upgrade numpy 

        if   field == 'Temperature'    : lmin =   0.9     ; lmax =   1.3
        elif field == 'Temperature_rms': lmin =   0.0     ; lmax =   0.1    ; levels_n  = 5
        elif field == 'Y'              : lmin =   0.0     ; lmax =   3.2e-3                 ; oo_magn = -3
        elif field == 'Sat_Ratio '     : lmin =   0.0     ; lmax = 900.0    ; levels_n  = 6
        elif field == 'Rho'            : lmin =   0.7     ; lmax =   1.0    ; levels_n  = 6
        elif field == 'U_r'            : lmin =  -0.03    ; lmax =   0.03   ; levels_n  = 6 ; oo_magn = -2
        elif field == 'fp_div'         : lmin = - 1.0e-7  ; lmax =   0.0    ; levels_n  = 5 ; oo_magn = -7
        elif field == 'fp_vap'         : lmin =   0.0     ; lmax =   5.0e-8 ; levels_n  = 5 ; oo_magn = -8
        elif field == 'Crit_Radius'    : lmin =   1.0e-7  ; lmax =   9.0e-7                 ; oo_magn = -7
        elif field == 'Crit_Radius_rms': lmin =   1.0e-7  ; lmax =   3.0e-7 ; levels_n  = 6 ; oo_magn = -7
        elif field == 'Nucl_rate'      : lmin =   1.0e-16 ; lmax =   1.0e-5 ; levels_n  = 6 ; oo_magn = -5

        if run_out == 'flu_mean' and np.isclose(max(mean_fields['Y']),1.8e-3, rtol=1e-4, atol=1e-4):

            if   field == 'J_rate'         : lmin =   0.0     ; lmax =   2.4    ; levels_n  = 6 ; oo_magn =  0
            elif field == 'fake_J_rate'    : lmin =   0.0     ; lmax =   2.4e-1 ; levels_n  = 6 ; oo_magn = -1
            elif field == 'Jrate_rms'      : lmin =   0.0     ; lmax =   0.8    ; levels_n  = 6 ; oo_magn =  0

        if run_out == 'flu_mean' and np.isclose(max(mean_fields['Y']),2.1e-3, rtol=1e-4, atol=1e-4):

            if   field == 'J_rate'         : lmin =   0.0     ; lmax =   3.6    ; levels_n  = 6 ; oo_magn =  0
            elif field == 'fake_J_rate'    : lmin =   0.0     ; lmax =   3.6    ; levels_n  = 6 ; oo_magn =  0
            elif field == 'Jrate_rms'      : lmin =   0.0     ; lmax =   1.2    ; levels_n  = 6 ; oo_magn =  0

        if   field == 'part_number' : lmin = 0.0 ; lmax = util.OOMRoundDown(lmax)    ; levels_n  = 5
        elif field == 'dexp_number' : lmin = 0.0 ; lmax = util.OOMRoundDown(lmax)    ; levels_n  = 5
        elif field == 'part_radius' : lmin = 0.0 ; lmax = util.OOMRoundUp  (lmax)    ; levels_n  = 5 ; oo_magn = util.OOMUp(lmax)
        elif field == 'mass_loading': lmin = 0.0 ; lmax = util.OOMRoundUp  (lmax)/10 ; levels_n  = 5 ; oo_magn = util.OOMUp(lmax)

        lcountour = np.linspace(lmin, lmax, levels_n + 1)

        if   field == 'Sat_Ratio'       : lcountour[0] = 1.0

        elif field == 'Temperature_rms' : lcountour    = np.insert(lcountour, 1, 0.01) #, axis = 0)

        elif field in ['Y','U_r','part_radius','mass_loading','fp_div'] : 

            cm_format = util.OOMFormatter    (oo_magn, mathText=False)

        elif field == 'fp_vap' :

            cm_format = util.OOMFormatter_eng(oo_magn, mathText=False)

        elif field in ['J_rate','Jrate_rms','fake_J_rate'] : 

            lcountour[0] = 1.e-3 ; cm_format = util.OOMFormatter    (oo_magn, mathText=False)

        elif field in ['Crit_Radius','Crit_Radius_rms'] : 

            cm_format = util.OOMFormatter(oo_magn, mathText=False)
    
        if field not in name_list: plt.title(field)
    
	# CONTOUR: draws the boundaries of the isosurfaces & fill the contour plots

        if run_out != 'prt_mean': cp = plt.contour (yi,xi,zi,levels=lcountour,linewidths=0.5,colors='black')

        cp = plt.contourf(yi,xi,zi,levels=lcountour,vmin=lmin,vmax=lmax,extend='both',cmap=my_cmap) #,cmap=cm.viridis) 
    
        plt.colorbar(cp,format=cm_format,orientation='horizontal') ; fig.tight_layout() # upgrade plt: location='bottom'
    
        plt.savefig(res_dir + '/%s.png'%field, format='png', dpi=300) ; plt.close('all')

    return []

"#################### MEAN PROFILE PLOT ####################" 

def profile_plots(run_dir,run_out,run_ver,res_dir,name_list):

    res_dir  = res_dir + '/1D_plots/'
    
    if util.chk_dir(res_dir) == False: os.makedirs(res_dir) # make sure the results directory exists 

    file_d = run_dir + '%s%s.vtk'%(run_ver,run_out) ; out_f = rosie.getOutputVTKwithPointDataFromFile(file_d)

    grid, mean_fields = rosie.getFields(out_f,name_list)

    for field in mean_fields:   

        print '\n--> plotting field %s'%field

        fig = plt.figure(figsize=(10, 5)) ; ax = fig.add_subplot(111) # ; ax.set_aspect(1, adjustable = 'box')
        
        x = grid[:,1] if run_out == "fra_mean" else grid[:,2] ; y = mean_fields[field]

        if run_out == "pvol_pdf": x = 2.*L_ref*x # expressed as (dimensional) diameter 

        xmin, xmax = min(x), max(x) ; ax.set_xlim(xmin, xmax)
        ymin, ymax = min(y), max(y) ; ax.set_ylim(ymin, ymax)

        field_name = field ; plt.title('%s'%field_name )

        if   run_out == "fra_mean": ax.set_xlabel('radial distance') ; ax.set_xlim  (xmin, 10.0)
        elif run_out == "fcl_mean": ax.set_xlabel('axial  distance')

        elif run_out == "prad_pdf": 
            field_name = "rad_" + field ; ax.set_xscale('log')
            ax.set_xlabel('non-dimensional particle radius')
            ax.set_ylabel('normalized pdf')
            plt.title    ('particle size distribution')

        elif run_out == "pexp_pdf": 
            field_name = "exp_" + field ; ax.set_xscale('log')
            ax.set_xlabel('non-dimensional particle radius')
            ax.set_ylabel('(filtered) normalized pdf')
            plt.title    ('particle size distribution - experimentally detected')

        elif run_out == "pvol_pdf": 
            field_name = "vol_" + field ; ax.set_xscale('log')
            ax.set_xlabel('dimensional particle diameter [m]')
            ax.set_ylabel('normalized volume fraction')
            plt.title    ('particle volume fraction distribution')

        elif run_out == "pvap_pdf": 
            field_name = "vap_" + field ; ax.set_xscale('log')
            ax.set_xlabel('non-dimensional particle radius')
            ax.set_ylabel('normalized pdf')
            plt.title    ('particle mass-flux distribution')

        ax.plot(x,y,'b-',label=r'%s'%field_name) 

        plt.legend(loc='best') ; plt.grid(which='both', axis='both',color='darkgrey')

        plt.savefig(res_dir + '/%s.png'%field_name, format='png', dpi=300) ; plt.close('all')

    return []

"#################### MEAN SCATTER PLOT ####################" 

def scatter_plots(run_dir,run_ver,res_dir,name_list):
                
    res_dir  = res_dir + '/2D_scatter/'
                    
    if util.chk_dir(res_dir) == False: os.makedirs(res_dir) 

    z_name_list  = ['J_rate','Rho'] ; xy_name_list = ['Temperature','Sat_ratio','Y_vapour','Rho']

    for z_name in z_name_list:
        
        for xy_name in list(combinations(xy_name_list, 2)): 

            x_name = xy_name[0] ; y_name = xy_name[1]

            if x_name == z_name or y_name == z_name: continue

            fig = plt.figure(figsize=(8, 5)) ; ax = fig.add_subplot(111) # ; ax.set_aspect(1, adjustable = 'box')
                
            name_list = [x_name, y_name, z_name]

            for root, dirs, files in os.walk("%s"%run_dir):
            
                for filename in files:

                    if re.match('%sJE*.*.vtk'%run_ver, filename):
            
                        file_d = run_dir + '%s'%(filename) ; out_f = rosie.getOutputVTKwithPointDataFromFile(file_d)
                   
                        grid, inst_fields = rosie.getFields(out_f,name_list)
            
                        x = inst_fields[x_name] ; y = inst_fields[y_name] ; z = inst_fields[z_name]

#                        if x_name == 'Y_vapour' or x_name == 'Temperature': x  = (x - min(x))/(max(x) - min(x)) # to restore asap 
#                        if y_name == 'Y_vapour' or y_name == 'Temperature': y  = (y - min(y))/(max(y) - min(y)) # to restore asap
                                        
                        xmin, xmax = min(x), max(x) # xmin, xmax = 1.e-3, 1.  

                        ax.set_xlim(xmin, xmax) ; ax.set_xlabel(r'%s'%x_name)

                        ymin, ymax = min(y), max(y) # ymin, ymax = 1.e-3, 2.e+3 

                        ax.set_ylim(ymin, ymax) ; ax.set_ylabel(r'%s'%y_name) # ; ax.set_yscale('log')

                        if y_name == 'Sat_ratio': ax.set_ylabel(r'Saturation ratio') ; ax.set_xscale('log')

                        vmin = min(z) ; vmax = max(z) ; norm = colors.Normalize(vmin=vmin, vmax=vmax)

                        if   z_name in ['Rate_nucl','J_rate','Crit_rad']:
                       
                            if   z_name == 'Rate_nucl': vmin = 1.0e-16 ; vmax = 1.0e-4
                            elif z_name == 'J_rate':    vmin = 1.0e-6  ; vmax = 1.0e+2
                            elif z_name == 'Crit_rad':  vmin = 1.0e-8

                            norm = colors.LogNorm  (vmin=vmin, vmax=vmax)
                                        
                        sc = plt.scatter(x, y, s=0.05, c=z, norm=norm, alpha=1.0, linewidths=0.0, edgecolors='face', rasterized=True, cmap = plt.cm.get_cmap('jet', 8))

                cbar = plt.colorbar(sc) ; cbar.set_label('%s'%z_name, rotation=270, labelpad=20)

                if z_name == 'J_rate': cbar.set_label('Nucleation rate', rotation=270, labelpad=15)
                
                plt.savefig(res_dir+'/%s_%s_%s.png'%(x_name[0:3:1],y_name[0:3:1],z_name[0:3:1]), format='png', dpi=300) ; plt.close('all')

                break   #prevent decending into subfolders
   
    return []

"#################### NORM PROFILE PLOT ####################" 

def norm_profile_plots(run_dir,run_out,norm_out,run_ver,res_dir):

    res_dir  = res_dir + '/1D_plots/'
    
    if util.chk_dir(res_dir) == False: os.makedirs(res_dir) # make sure the results directory exists 

    name_list  = ['Uz_centline','R_halfw_Uz','Y_centline','R_halfw_Y']

    file_d = run_dir + '%s%s.vtk'%(run_ver,norm_out) ; out_f = rosie.getOutputVTKwithPointDataFromFile(file_d)

    axi_grid, mean_axi_fields = rosie.getFields(out_f,name_list)

    name_list  = []

    file_d = run_dir + '%s%s.vtk'%(run_ver,run_out) ; out_f = rosie.getOutputVTKwithPointDataFromFile(file_d)

    rad_grid, mean_rad_fields = rosie.getFields(out_f,name_list)

    for field in mean_rad_fields:   

        fig = plt.figure(figsize=(9, 3)) ; ax = fig.add_subplot(111) # ; ax.set_aspect(1, adjustable = 'box')
        
        if run_out == "fra_mean": x = rad_grid[:,1]

        y = mean_rad_fields[field]
        
        rad_field = field.split("_",2) ; rad_field_name = rad_field[0] ; rad_field_dist = float(rad_field[1])

        axi_norm = mean_axi_fields[rad_field_name + "_centline"] ; rad_norm = mean_axi_fields["R_halfw_" + rad_field_name]

        for i in range(len(axi_norm)):

            if np.isclose(axi_grid[i,2],rad_field_dist, rtol=1e-2, atol=1e-2): # use util.find_nearest

                y = y/axi_norm[i] ; x = x/rad_norm[i]

                break

        xmin, xmax = min(x), max(x) ; ax.set_xlim(xmin, xmax)
        ymin, ymax = min(y), max(y) ; ax.set_ylim(ymin, ymax)

        field_name = field

        if run_out == "fra_mean": ax.set_xlim(xmin, 3.0)

#        ax.set_ylim(xmin, 1.0)

        plt.title('norm_%s'%field_name ) ; ax.plot(x,y,'b-',label=r'norm_%s'%field_name) 

        ax.set_xlabel(r'%s'%field) ; ax.set_ylabel(field_name)

        plt.legend(loc='best') ; plt.grid(which='both', axis='both',color='darkgrey')

        plt.savefig(res_dir + '/norm_%s.png'%field_name, format='png', dpi=300) ; plt.close('all')

    return []


"#################### MEAN PROFILE PLOT ####################" 

def cfr_profile_plots(run_dir,run_out,run_ver,res_dir):

    res_dir = res_dir + '/1D_plots/'
    
    if util.chk_dir(res_dir) == False: os.makedirs(res_dir) # make sure the results directory exists 

    fig = plt.figure(figsize=(8, 5)) ; ax = fig.add_subplot(111) # ; ax.set_aspect(1, adjustable = 'box')

    xmin =   1.e+16 ; xmax = - 1.e+16
    ymin =   1.e+16 ; ymax = - 1.e+16

    name_list  = []

    for root, dirs, files in os.walk("%s"%run_dir):

        for filename in files:

            if re.match('%s%s*.*.vtk'%(run_ver,run_out), filename):

                print '\n--> looking for %s'%filename 

                fname_end = filename.split("_") ; prof_case = fname_end[-1] ; prof_case = prof_case.rstrip('.vtk')

                print '\n--> case %s'%prof_case 

                file_d = run_dir + '%s'%(filename) ; out_f = rosie.getOutputVTKwithPointDataFromFile(file_d)

                grid, mean_fields = rosie.getFields(out_f,name_list)

                for field in mean_fields:   

                    if field == "pdf_38.0":

                        x = grid[:,2] ; xmin, xmax = 2.e-7, 2.e-6 # ; xmin, xmax = min(x), max(x)

                        ax.set_xlim(xmin, xmax) ; ax.set_xscale('log')# ; ax.set_xlabel(r'%s'%field)

                        y = mean_fields[field] ; y = y*(4/3*3.14*np.power(x,4)) ; x = 2*x*0.00175

#                        for i in range(len(y)): 
#
#                        	if i < (len(y)-1):
#				
#					delta_x = x[i+1] - x[i]
#
#                        	else:
#
#					delta_x = x[i] - x[i-1]
#
#				y[i] = y[i]/(math.log10(x[i] + delta_x) - math.log10(x[i] - delta_x)) 
#
##				print y[i]

                        ymin = min(np.append(y,ymin)) ; ymax = max(np.append(y,ymax))

                        ax.set_ylim(ymin, ymax) # ; ax.set_yscale('log')

                        field_name = field

                        if   filename == "fra_mean": ax.set_xlim(xmin, 10.0)
                        elif filename == "prad_pdf" + prof_case + ".vtk": field_name = "rad_" + field ; ax.set_xscale('log') # ; ax.set_yscale('log')
                        elif filename == "pexp_pdf" + prof_case + ".vtk": field_name = "exp_" + field ; ax.set_xscale('log') # ; ax.set_yscale('log')
                        elif filename == "pvol_pdf" + prof_case + ".vtk": field_name = "vol_" + field ; ax.set_xscale('log') # ; ax.set_yscale('log')
                        elif filename == "pvap_pdf": field_name = "vap_" + field ; ax.set_xscale('log') #  ax.plot(x,y,'b-',label=r'%s'%prof_case) 

                        if prof_case == "x18": ax.plot(x,y,'#0d4186',label=r'1.8',linewidth=2.0) # ax.plot(x,y,'#3470F0',label=r'%s'%prof_case) 
                        if prof_case == "x21": ax.plot(x,y,'#0d4186',label=r'2.1',linewidth=2.0) # ax.plot(x,y,'#3470F0',label=r'%s'%prof_case) 
                        if prof_case == "x25": ax.plot(x,y,'#0d4186',label=r'2.5',linewidth=2.0) # ax.plot(x,y,'#3470F0',label=r'%s'%prof_case) 
                        if prof_case == "x32": ax.plot(x,y,'#0d4186',label=r'3.2',linewidth=2.0) # ax.plot(x,y,'#3470F0',label=r'%s'%prof_case) 

        ax.yaxis.set_major_formatter(util.OOMFormatter(-10, mathText=False))

	xvals = [8.0e-7,7.5e-7,8.0e-7,9.0e-7] ; util.labelLines(plt.gca().get_lines(),align=False,xvals=xvals,color='k') #,fontsize=14)

        plt.grid(which='both', axis='both',color='darkgrey') ; fig.tight_layout()

        plt.savefig(res_dir + '/cfr_%s.png'%field_name, format='png', dpi=300) ; plt.close('all')

        break   #prevent decending into subfolders

    return []

"#################### MEAN PART NUMBER PLOT ####################" 

def cfr_DNSvsExp(run_dir,plot_name,run_ver,res_dir):

    name_list = [] ; file_list = [] ; run_vers = [run_ver,"g"] #; run_vers = list(run_ver).append("g")

    res_dir = res_dir + '2D_scatter/'

    if util.chk_dir(res_dir) == False: os.makedirs(res_dir) 
    
    x_exp=[] ; y_exp=[]

    x_dns_2w_shrt = [] ; y_dns_2w_shrt = [] ; y_dns_2w_shrt_norm = [] ; y_dns_2w_shrt_dime = [] 
    x_dns_2w_long = [] ; y_dns_2w_long = [] ; y_dns_2w_long_norm = [] ; y_dns_2w_long_dime = []
    x_dns_1w_long = [] ; y_dns_1w_long = [] ; y_dns_1w_long_norm = [] ; y_dns_1w_long_dime = []

    x_dns_2w_long_okuyama = [] ; y_dns_2w_long_okuyama = [] ; y_dns_2w_long_norm_okuyama = [] ; y_dns_2w_long_dime_okuyama = []
    x_dns_2w_long_hameri  = [] ; y_dns_2w_long_hameri  = [] ; y_dns_2w_long_norm_hameri  = [] ; y_dns_2w_long_dime_hameri  = []


    Wvap = 278.34e-3 ; Wgas = 28.29e-3 ; Lref_exp = 0.00375/2. ; Lref_dns = 0.00175 # pipe dimensional radius (DNS)
    
    norm_exp  = (Lref_exp*1.e2)**3 ; norm_dns  = (Lref_dns*1.e2)**3

    if plot_name == "fried2000F10.11": exp_dataset = 'DBP_Re4700_d0375cm.csv' ; dns_dataset = 'prt_mean'

    for root, dirs, files in os.walk("%s"%run_dir):
    
        for filename in files:

            if re.match('%s'%exp_dataset, filename):
    
                exp_dataset = '%s'%os.path.join(root, filename); print '\n--> Reading', exp_dataset
#                break   #prevent decending into subfolders
    
                with open(exp_dataset, 'r') as csvfile:
                
                    plots = csv.reader(csvfile, delimiter=';')

                    for row in plots:
                        x_exp.append(float(row[0]))
                        y_exp.append(float(row[1]))
                
                y_exp_dime = [i * ((2.0)**3) for i in y_exp]
                
                y_exp_norm = [i * norm_exp for i in y_exp_dime]

            if re.match('%s%s.vtk'%(run_vers,dns_dataset), filename): file_list.append(os.path.join(root, filename))

    print "the file list is : ", file_list
    print "the run ver is : ", run_vers

    for dns_dataset in file_list: 


        file_d = dns_dataset ; out_f = rosie.getOutputVTKwithPointDataFromFile(file_d)

        grid, mean_fields = rosie.getFields(out_f,name_list)

        x, y = grid[:,1], grid[:,2] ; npts = 500 ; xmin, xmax = min(x), max(x) ; ymin, ymax = min(y), max(y)

        xi = np.linspace(xmin, xmax, npts) ; yi = np.linspace(ymin, ymax, npts) ; xi, yi = np.meshgrid(xi, yi)

        z = mean_fields['part_number'] ; zi = griddata((x, y), z, (xi, yi), method='linear')

        y_dns = util.ProbeAtLocation(zi, xi, yi, 0.5, 40.0) ; print "probed value is :", y_dns

        if   re.match('.+okuyama', file_d): y_dns_2w_long_okuyama.append(y_dns) ; y_dns_2w_long_norm_okuyama = y_dns_2w_long_okuyama ; \
                                            y_dns_2w_long_dime_okuyama = [(i / norm_dns )* ((0.8)**3) for i in y_dns_2w_long_norm_okuyama]

        if   re.match('.+hameri', file_d): y_dns_2w_long_hameri.append(y_dns) ; y_dns_2w_long_norm_hameri = y_dns_2w_long_hameri ; \
                                            y_dns_2w_long_dime_hameri = [(i / norm_dns )* ((0.8)**3) for i in y_dns_2w_long_norm_hameri]

        if   re.match('.+shrt.+', file_d):

            if   re.match('.+0w.+', file_d): y_dns_0w_shrt.append(y_dns) ; y_dns_0w_shrt_norm = y_dns_0w_shrt ; \
                                             y_dns_0w_shrt_dime = [(i / norm_dns )* ((0.8)**3) for i in y_dns_0w_shrt_norm]

            elif re.match('.+1w.+', file_d): y_dns_1w_shrt.append(y_dns) ; y_dns_1w_shrt_norm = y_dns_1w_shrt ; \
                                             y_dns_1w_shrt_dime = [(i / norm_dns )* ((0.8)**3) for i in y_dns_1w_shrt_norm]

            elif re.match('.+2w.+', file_d): y_dns_2w_shrt.append(y_dns) ; y_dns_2w_shrt_norm = y_dns_2w_shrt ; \
                                             y_dns_2w_shrt_dime = [(i / norm_dns )* ((0.8)**3) for i in y_dns_2w_shrt_norm]

        elif re.match('.+long.+', file_d):

            if   re.match('.+0w.+', file_d): y_dns_0w_long.append(y_dns) ; y_dns_0w_long_norm = y_dns_0w_long ; \
                                             y_dns_0w_long_dime = [(i / norm_dns )* ((0.8)**3) for i in y_dns_0w_long_norm]

            elif re.match('.+1w.+', file_d): y_dns_1w_long.append(y_dns) ; y_dns_1w_long_norm = y_dns_1w_long ; \
                                             y_dns_1w_long_dime = [(i / norm_dns )* ((0.8)**3) for i in y_dns_1w_long_norm]

            elif re.match('.+2w.+', file_d): y_dns_2w_long.append(y_dns) ; y_dns_2w_long_norm = y_dns_2w_long ; \
                                             y_dns_2w_long_dime = [(i / norm_dns )* ((0.8)**3) for i in y_dns_2w_long_norm]


        file_d = dns_dataset.replace('prt_mean', 'flu_mean') ; out_f = rosie.getOutputVTKwithPointDataFromFile(file_d)

        grid, mean_fields = rosie.getFields(out_f,name_list)

        x, y = grid[:,1], grid[:,2] ; npts = 500 ; xmin, xmax = min(x), max(x) ; ymin, ymax = min(y), max(y)

        xi = np.linspace(xmin, xmax, npts) ; yi = np.linspace(ymin, ymax, npts) ; xi, yi = np.meshgrid(xi, yi)

        z = mean_fields['Y'] ; zi = griddata((x, y), z, (xi, yi), method='linear')

        x_dns = util.ProbeAtLocation(zi, xi, yi, 0.0, 0.1) ; x_dns = x_dns/(x_dns + (1-x_dns)*Wvap/Wgas)

        if   re.match('.+okuyama', file_d): x_dns_2w_long_okuyama.append(x_dns)
        if   re.match('.+hameri', file_d): x_dns_2w_long_hameri.append(x_dns)

        if   re.match('.+shrt.+', file_d):

            if   re.match('.+0w.+', file_d): x_dns_0w_shrt.append(x_dns)
            elif re.match('.+1w.+', file_d): x_dns_1w_shrt.append(x_dns)
            elif re.match('.+2w.+', file_d): x_dns_2w_shrt.append(x_dns)

        elif re.match('.+long.+', file_d):

            if   re.match('.+0w.+', file_d): x_dns_0w_long.append(x_dns)
            elif re.match('.+1w.+', file_d): x_dns_1w_long.append(x_dns)
            elif re.match('.+2w.+', file_d): x_dns_2w_long.append(x_dns)


    fig = plt.figure(figsize=(12, 6))

    ###############################################################
    
    ax_dime  = fig.add_subplot(121) ; ax_dime.set_yscale('log')
    
    xmin_dime, xmax_dime = 0., 5.e-4 ; ymin_dime, ymax_dime = 1., 1.e+6
    
    ax_dime.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))

    ax_dime.set_xlim(xmin_dime, xmax_dime) ; ax_dime.set_ylim(ymin_dime, ymax_dime)
    
    plt.scatter(x_exp      , y_exp_dime      , c='b', label='experiments ')
#    plt.scatter(x_dns_2w_shrt, y_dns_2w_shrt_dime, c='r', label='2w DNS short')
    plt.scatter(x_dns_2w_long, y_dns_2w_long_dime, c='m', label='2w DNS long ')
#    plt.scatter(x_dns_1w_long, y_dns_1w_long_dime, c='g', label='1w DNS long ')

    plt.scatter(x_dns_2w_long_okuyama, y_dns_2w_long_dime_okuyama, c='orange', label='2w DNS long Okuyama')
    plt.scatter(x_dns_2w_long_hameri,  y_dns_2w_long_dime_hameri,  c='red',    label='2w DNS long Hameri')
    
    plt.grid(True) ; plt.legend(loc="lower right") ; plt.title('dimensional results')
    
    plt.xlabel('DBP inlet molar-franction') ; plt.ylabel('Particle Number Density (#/cc)')
    
    #plt.savefig(res_dir + '/dime_%s.eps'%plot_name, format='eps', dpi=300)
    
    ###############################################################
    
    #fig= plt.figure(figsize=(6, 6))
    
    ax_norm  = fig.add_subplot(122) ; ax_norm.set_yscale('log')
    
    xmin_norm, xmax_norm = 0.   , 5.e-4 ; ymin_norm, ymax_norm = 1.e-2, 1.e+4
    
    ax_norm.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))

    ax_norm.set_xlim(xmin_norm, xmax_norm) ; ax_norm.set_ylim(ymin_norm, ymax_norm)
    
    plt.scatter(x_exp      , y_exp_norm      , c='b', label='experiments ')
#    plt.scatter(x_dns_2w_shrt, y_dns_2w_shrt_norm, c='r', label='2w DNS short')
    plt.scatter(x_dns_2w_long, y_dns_2w_long_norm, c='m', label='2w DNS long ')
#    plt.scatter(x_dns_1w_long, y_dns_1w_long_norm, c='g', label='1w DNS long ')

    plt.scatter(x_dns_2w_long_okuyama, y_dns_2w_long_norm_okuyama, c='orange', label='2w DNS long Okuyama')
    plt.scatter(x_dns_2w_long_hameri,  y_dns_2w_long_norm_hameri,  c='red',    label='2w DNS long Hameri')
    
    plt.grid(True) ; plt.legend(loc="lower right") ; plt.title('non-dimensional results')
    
    plt.xlabel('DBP inlet molar-franction') ; plt.ylabel('Particle Number Density (#/dV)')
    
    fig.suptitle('DNS vs Experimental (Friedlander)', fontsize=14)
    
    plt.savefig(res_dir + '/cfrDNSvsExp_%s.png'%plot_name, format='png', dpi=300) ; plt.close('all')
    
    return []

"#################### MASS TRANSFER MODELS COMPARISON PLOT ####################" 

def mdot_HKvsMA(run_dir,run_ver,res_dir):

    res_dir = res_dir + '1D_plots/'

    if util.chk_dir(res_dir) == False: os.makedirs(res_dir) 


    xmin = -16 ; xmax = -3 ; npts = 14 ; Y_vap = np.logspace(xmin, xmax, npts)

    xmin = -7 ; xmax = -2 ; npts = 25 ; R_prt = np.logspace(xmin, xmax, npts)


    fig = plt.figure(figsize=(12, 6))

        
    ax_mdot = fig.add_subplot(121) ; ax_mdot.set_xscale('log') ; ax_mdot.set_yscale('log')
    
    xmin_mdot, xmax_mdot = 1.e-7, 1.e-2 ; ymin_mdot, ymax_mdot = 1.e-6, 1.e-2
    
    ax_mdot.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))

#    ax_mdot.set_xlim(xmin_mdot, xmax_mdot) ; ax_mdot.set_ylim(ymin_mdot, ymax_mdot)


    for Y in Y_vap:

        theta = 1. + 200.*Y

        X_sat = np.exp((Lat_heat/Runi)*(1./theta_boi - 1./(theta*theta_ref))) 

        Y_sat = X_sat/(X_sat + (1-X_sat)*Wgas/Wvap)

        X = Y/(Y + (1-Y)*Wvap/Wgas)

        mdot_MA = [] ; mdot_HK = []

        mdot_MA = [abs(2.*np.pi*R*Sh/(Re*Sc)*(Y - Y_sat)) for R in R_prt]

        mdot_HK = [abs((alpha/Ma)*(R**2)*np.sqrt(8.*np.pi/(Rvap*theta))*(X - X_sat)) for R in R_prt]

        ###############################################################

        R_cross = abs(2.*np.pi*Sh/(Re*Sc)*(Y - Y_sat)) / abs((alpha/Ma)*np.sqrt(8.*np.pi/(Rvap*theta))*(X - X_sat))
        m_cross = abs(2.*np.pi*R_cross*Sh/(Re*Sc)*(Y - Y_sat))

        ax_mdot.plot(R_cross, m_cross,'r*',label='Y = %f'%Y)

        ax_mdot.plot(R_prt, mdot_MA,'b-',label='Y = %f'%Y) ; ax_mdot.plot(R_prt, mdot_HK,'g-',label='Y = %f'%Y)

        
    plt.grid(True) ; plt.title('Particle mass RHS') ; plt.xlabel('DBP local mass-franction') 
        
    ax_rdot = fig.add_subplot(122) ; ax_rdot.set_xscale('log') ; ax_rdot.set_yscale('log')
   
    xmin_rdot, xmax_rdot = 1.e-7, 1.e-2 ; ymin_rdot, ymax_rdot = 1.e-10, 1.e-6
    
    ax_rdot.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))

#    ax_rdot.set_xlim(xmin_rdot, xmax_rdot) ; ax_rdot.set_ylim(ymin_rdot, ymax_rdot)
        

    for Y in Y_vap:

        theta = 1. + 200.*Y

        X_sat = np.exp((Lat_heat/Runi)*(1./theta_boi - 1./(theta*theta_ref))) 

        Y_sat = X_sat/(X_sat + (1-X_sat)*Wgas/Wvap)

        X = Y/(Y + (1-Y)*Wvap/Wgas)

        rdot_MA = [] ; rdot_HK = []

        rdot_MA = [abs(Sh/(Re*Sc)/(2.*rho_liq*R)*(Y - Y_sat)) for R in R_prt]

        rdot_HK = [abs((alpha/(Ma*rho_liq))*np.sqrt(1./(2.*np.pi)/(Rvap*theta))*(X - X_sat)) for R in R_prt]

        ###############################################################

        R_cross = abs(Sh/(Re*Sc)/(2.)*(Y - Y_sat)) / abs((alpha/(Ma))*np.sqrt(1./(2.*np.pi)/(Rvap*theta))*(X - X_sat)) 
        m_cross = abs(Sh/(Re*Sc)/(2.*rho_liq*R_cross)*(Y - Y_sat))

        ax_rdot.plot(R_cross, m_cross,'r*',label='Y = %f'%Y)

        ax_rdot.plot(R_prt, rdot_MA,'b-',label='Y = %f'%Y) ; ax_rdot.plot(R_prt, rdot_HK,'g-',label='Y = %f'%Y)
        
    plt.grid(True); plt.title('Particle radius RHS') ; plt.xlabel('DBP local mass-franction')
    
        
    fig.suptitle('Mass analogy vs Hertz-Knudsen mass transfer : approx. Temperature ', fontsize=14)
        
    plt.savefig(res_dir + '/mdot_HKvsMA.png', format='png', dpi=300) ; plt.close('all')

                

    res_dir  = res_dir + '../2D_scatter/'
                    
    if util.chk_dir(res_dir) == False: os.makedirs(res_dir) 


    fig = plt.figure(figsize=(18, 12))

    sc_mdot  = fig.add_subplot(221) ; sc_rdot = fig.add_subplot(222)
    sc_temp  = fig.add_subplot(223) ; sc_cros = fig.add_subplot(224)

    sc_mdot.set_xscale('log'); sc_mdot.set_yscale('log') ; sc_mdot.set_xlabel('particle radius') ; sc_mdot.set_ylabel('mass transfer')
    sc_rdot.set_xscale('log'); sc_rdot.set_yscale('log') ; sc_rdot.set_xlabel('particle radius') ; sc_rdot.set_ylabel('radius RHS')

#    print 'lines.markersize = ', matplotlib.rcParams['lines.markersize']

    name_list = ['Y_vapour','Temperature']

    for root, dirs, files in os.walk("%s"%run_dir):
    
        for filename in files:

            if re.match('%sJE*.*.vtk'%run_ver, filename):
    
                file_d = run_dir + '%s'%(filename) ; out_f = rosie.getOutputVTKwithPointDataFromFile(file_d)
           
                grid, inst_fields = rosie.getFields(out_f,name_list)
    
                Y_vap = inst_fields['Y_vapour'] ; Tempera = inst_fields['Temperature']


                X_sat = [np.exp((Lat_heat/Runi)*(1./theta_boi - 1./(theta*theta_ref))) for theta in Tempera] 

                Y_sat = [X/(X + (1-X)*Wgas/Wvap) for X in X_sat]
                X_vap = [Y/(Y + (1-Y)*Wvap/Wgas) for Y in Y_vap]

                sqrtInvTempera = [np.sqrt(1./theta) for theta in Tempera]

                number = 0 ; mdot_MA = [] ; mdot_HK = [] ; rdot_MA = [] ; rdot_HK = [] ; Rp = np.copy(Y_vap)

                for R in R_prt:

                    mdot_MA = abs(Sh/(Re*Sc)* 2.*np.pi  *R *np.subtract(Y_vap,Y_sat))
                    rdot_MA = abs(Sh/(Re*Sc)/(2.*rho_liq*R)*np.subtract(Y_vap,Y_sat))
    
                    mdot_HK = abs((alpha/ Ma)*(R**2)  *np.sqrt(8. * np.pi/Rvap) *np.multiply(sqrtInvTempera,np.subtract(X_vap,X_sat)))
                    rdot_HK = abs((alpha/(Ma*rho_liq))*np.sqrt(0.5/(np.pi*Rvap))*np.multiply(sqrtInvTempera,np.subtract(X_vap,X_sat)))

		    print 'min/max MA m_dot : ', np.amin(mdot_MA), np.amax(mdot_MA)
		    print 'min/max MA r_dot : ', np.amin(rdot_MA), np.amax(rdot_MA)
		    print 'min/max HK m_dot : ', np.amin(mdot_HK), np.amax(mdot_HK)
		    print 'min/max HK r_dot : ', np.amin(rdot_HK), np.amax(rdot_HK)

		    Rp.fill(R*1.05)

		    sc_mdot.scatter(Rp, mdot_MA, s=1.0, c='b', alpha=1.e-2, edgecolors='none', rasterized=True)
		    sc_rdot.scatter(Rp, rdot_MA, s=1.0, c='b', alpha=1.e-2, edgecolors='none', rasterized=True)
                                                                                                            
                    Rp.fill(R*0.95)

		    sc_mdot.scatter(Rp, mdot_HK, s=1.0, c='r', alpha=1.e-2, edgecolors='none', rasterized=True)
		    sc_rdot.scatter(Rp, rdot_HK, s=1.0, c='r', alpha=1.e-2, edgecolors='none', rasterized=True)

                    print 'radius : ', R , ' number : ', number ; number += 1 


                R_cross = [] ; R_cross = np.divide( abs(Sh/(Re*Sc)/(2.)*np.subtract(Y_vap,Y_sat)) , abs((alpha/Ma)*np.sqrt(1./(2.*np.pi*Rvap))*np.multiply(sqrtInvTempera,np.subtract(X_vap,X_sat))) ) 


                vmin = min(R_cross) ; vmax = max(R_cross) ; norm = colors.Normalize(vmin=vmin, vmax=vmax) #; norm = colors.LogNorm  (vmin=vmin, vmax=vmax)

                sc = sc_temp.scatter(Y_vap, Tempera, s=1.0, c=R_cross, norm=norm, alpha=0.1, edgecolors='none', rasterized=True)

                sc_temp.set_xlim(min(Y_vap  ), max(Y_vap  )) ; sc_temp.set_xlabel('Y')
		sc_temp.set_ylim(min(Tempera), max(Tempera)) ; sc_temp.set_ylabel('Temperature')

        	oo_magn = -5  ; cm_format = util.OOMFormatter(oo_magn, mathText=False)

                cbar = plt.colorbar(sc, ax=sc_temp, format=cm_format) ; cbar.set_label('cross-over radius', rotation=270, labelpad=20)


                vmin = min(Tempera) ; vmax = max(Tempera) ; norm = colors.Normalize(vmin=vmin, vmax=vmax)

                sc = sc_cros.scatter(Y_vap, R_cross, s=1.0, c=Tempera, norm=norm, alpha=0.1, edgecolors='none', rasterized=True)

                sc_cros.set_xlim(min(Y_vap  ), max(Y_vap  )) ; sc_cros.set_xlabel('Y')
		sc_cros.set_ylim(min(R_cross), max(R_cross)) ; sc_cros.set_ylabel('cross-over radius')

        	oo_magn = -5 ; sc_cros.yaxis.set_major_formatter(util.OOMFormatter(oo_magn, mathText=False))

                cbar = plt.colorbar(sc, ax=sc_cros) ; cbar.set_label('Temperature', rotation=270, labelpad=20)

        plt.savefig(res_dir + '/scatter_HKvsMA.png', format='png', dpi=300) ; plt.close('all')

        break   #prevent decending into subfolders
    
    return []

"#################### NUCLEATION RATE COMPARISON VARYING SIGMA PLOT ####################" 

def Jrate(Rho, Y_vap, Sigma, Tempera, Sat_R):

    SupSat = 1.1 ; Rho = Rho[Sat_R > SupSat] ; Y_vap = Y_vap[Sat_R > SupSat] ; Sigma = Sigma[Sat_R > SupSat] ; Tempera = Tempera[Sat_R > SupSat] ; Sup_S = Sat_R[Sat_R > SupSat]

    a_0 =     np.log(L_ref**4/U_ref) 
    a_1 = 2. *np.log(Rho*Y_vap)
    a_2 = 0.5*np.log(np.divide(Sigma , (np.pi*(Mmol**3)*((rho_liq*rho_ref)**2))/2))

    b_0 = 16./3.*np.pi*(Wvap/(np.sqrt(Kbol)*Runi*rho_liq*rho_ref))**2
    b_1 = np.power(np.divide(Sigma   , Tempera*theta_ref) , 3)
    b_2 = np.power(np.log(Sup_S), -2)

    print 'a_0 min/max', np.amin(a_0), np.amax(a_0), 'a_1 min/max', np.amin(a_1), np.amax(a_1), 'a_2 min/max', np.amin(a_2), np.amax(a_2)
    print 'b_0 min/max', np.amin(b_0), np.amax(b_0), 'b_1 min/max', np.amin(b_1), np.amax(b_1), 'b_2 min/max', np.amin(b_2), np.amax(b_2)

    J= np.zeros_like(Sat_R) ; J[Sat_R > SupSat] = np.exp(a_0 + a_1 + a_2 + ( - b_0*b_1*b_2))

    J[Sat_R <= SupSat] = 0

    J[J <= 1.e-16] = 0.

    print ''
    print 'Sigma min/max', np.amin(Sigma), np.amax(Sigma)
    print 'Jrate min/max', np.amin(J), np.amax(J)
    print ''

    return J 

def Jrate_sigma(run_dir,run_ver,res_dir):

    res_dir = res_dir + '2D_scatter/'

    if util.chk_dir(res_dir) == False: os.makedirs(res_dir) 


    fig = plt.figure(figsize=(18, 9))
        
    sc_Jrat = fig.add_subplot(111) ; sc_Jrat.set_yscale('log')
    
    xmin_Jrat, xmax_Jrat = 2.7e-2, 3.4e-2 ; ymin_Jrat, ymax_Jrat = 1.e-6, 50.
    
#    sc_Jrat.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))

    sc_Jrat.set_xlim(xmin_Jrat, xmax_Jrat) ; sc_Jrat.set_ylim(ymin_Jrat, ymax_Jrat)


    name_list = ['Y_vapour','Temperature','Rho']

    for root, dirs, files in os.walk("%s"%run_dir):
    
        for filename in files:

            if re.match('%sJE*.*.vtk'%run_ver, filename):
    
                file_d = run_dir + '%s'%(filename) ; out_f = rosie.getOutputVTKwithPointDataFromFile(file_d)
           
                grid, inst_fields = rosie.getFields(out_f,name_list)
    
                Y_vap = inst_fields['Y_vapour'] ; Tempera = inst_fields['Temperature'] ; Rho = inst_fields['Rho']


                X_sat = [np.exp((Lat_heat/Runi)*(1./theta_boi - 1./(theta*theta_ref))) for theta in Tempera] 

                Y_sat = [X/(X + (1-X)*Wgas/Wvap) for X in X_sat]
                X_vap = [Y/(Y + (1-Y)*Wvap/Wgas) for Y in Y_vap]

		Sat_R = np.divide(X_vap , X_sat)

    		print ''
		print 'calculated X_vap min/max', np.amin(X_vap), np.amax(X_vap)
		print 'calculated X_sat min/max', np.amin(X_sat), np.amax(X_sat)
		print 'calculated Sat_R min/max', np.amin(Sat_R), np.amax(Sat_R)
    		print ''

		Sigma_Okuyama = np.array([1.e-3*(35.30 - 0.0863*(theta*theta_ref - theta_zero)) for theta in Tempera])
		Sigma_Bedanov = np.array([1.e-3*(35.72 - 0.0894*(theta*theta_ref - theta_zero)) for theta in Tempera])

		theta_crit = 781. ; Sigma_AIChE   = np.array([0.059663*(1. - theta*theta_ref/theta_crit)**1.2457 for theta in Tempera])
		theta_crit = 784. ; Sigma_Ambrose = np.array([0.059663*(1. - theta*theta_ref/theta_crit)**1.2457 for theta in Tempera])
		theta_crit = 823. ; Sigma_Hameri  = np.array([0.059663*(1. - theta*theta_ref/theta_crit)**1.2457 for theta in Tempera])

		Jrate_AIChE    = Jrate(Rho, Y_vap, Sigma_AIChE  , Tempera, Sat_R)
		Jrate_Ambrose  = Jrate(Rho, Y_vap, Sigma_Ambrose, Tempera, Sat_R)
		Jrate_Okuyama  = Jrate(Rho, Y_vap, Sigma_Okuyama, Tempera, Sat_R)
		Jrate_Bedanov  = Jrate(Rho, Y_vap, Sigma_Bedanov, Tempera, Sat_R)
		Jrate_Hameri   = Jrate(Rho, Y_vap, Sigma_Hameri , Tempera, Sat_R)

    		print ''
    		print 'Jrate_AIChE   min/max', np.amin(Jrate_AIChE  ), np.amax(Jrate_AIChE  )
    		print 'Jrate_Ambrose min/max', np.amin(Jrate_Ambrose), np.amax(Jrate_Ambrose)
    		print 'Jrate_Okuyama min/max', np.amin(Jrate_Okuyama), np.amax(Jrate_Okuyama)
    		print 'Jrate_Bedanov min/max', np.amin(Jrate_Bedanov), np.amax(Jrate_Bedanov)
    		print 'Jrate_Hameri  min/max', np.amin(Jrate_Hameri ), np.amax(Jrate_Hameri )
    		print ''

                sc = sc_Jrat.scatter(Sigma_AIChE  , Jrate_AIChE  , s=0.5, c='green' , alpha=1., edgecolors='none', rasterized=True, label='AIChE  ')
                sc = sc_Jrat.scatter(Sigma_Ambrose, Jrate_Ambrose, s=0.5, c='orange', alpha=1., edgecolors='none', rasterized=True, label='Ambrose')
                sc = sc_Jrat.scatter(Sigma_Okuyama, Jrate_Okuyama, s=0.5, c='blue'  , alpha=1., edgecolors='none', rasterized=True, label='Okuyama')
                sc = sc_Jrat.scatter(Sigma_Bedanov, Jrate_Bedanov, s=0.5, c='red'   , alpha=1., edgecolors='none', rasterized=True, label='Bedanov')
                sc = sc_Jrat.scatter(Sigma_Hameri , Jrate_Hameri , s=0.5, c='black' , alpha=1., edgecolors='none', rasterized=True, label='Hameri ')

	plt.grid(True) ; plt.legend(loc="upper left", markerscale=5., scatterpoints=1, fontsize=10) ; plt.title('non-dimensional Nucleation rates varying the Surface tension')

	sc_Jrat.set_xlabel('Surface tension [N/m]') ; sc_Jrat.set_ylabel('Nucleation rate [#]')
	
        plt.savefig(res_dir + '/scatter_Jrate_Sigma.png', format='png', dpi=300) ; plt.close('all')

        break   #prevent decending into subfolders
    
    return []

