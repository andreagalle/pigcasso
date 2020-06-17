# -*- coding: utf-8 -*-
"""
PPchan.py
"""

import modUtils as util
import os, sys, shutil
import draft as fct
import vtk
import numpy as np
import scipy.integrate
import toolCase as toca 
import time
from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt


"############################  3D CHANNEL .vtu  #####################################"

# calculation of the main quantities to be PP, this function retrieve mean profiles from the .vtu file

def getChannel(workdir, chan, model): # workdir='string/path', chan='string.3Dfile.vtu'

	workdir = workdir + '/%s_cc/vtu/'%model
	resdir    = workdir + '/../../PP_results/profiles_plot/'
	check_dir = util.fct_dir_exist(resdir)        
	
	if check_dir == False:
		os.makedirs(resdir) # make sure the results directoru exists 

	file_d = workdir + chan
    
	ori  = [3,0,0]
	norm = [0,0,1]

	out_f     = fct.getOutputVTKwithPointDataFromFile(file_d)
	Slice_out = fct.getSlice(out_f, ori, norm)

	ori1  = [3,0,0]
	norm1 = [1,0,0]

	"############################  GET THE FIELDS  #####################################"

	gamma = 1.4
	R = 0.714286

	# scalar fields and bulk data calculation
	# [y_coo, P]      = util.getVarDataAtPositionNoWall(Slice_out, ori1, norm1, 'P')
	[y_coo, T]      = util.getVarDataAtPositionNoWall(Slice_out, ori1, norm1, 'T')
	[y_coo, Rho]    = util.getVarDataAtPositionNoWall(Slice_out, ori1, norm1, 'RHO')
	[y_coo, Mu]     = util.getVarDataAtPositionNoWall(Slice_out, ori1, norm1, 'MU_LAM') # not avg because it's almost constant
	[y_coo, Mu_sgs] = util.getVarDataAtPositionNoWall(Slice_out, ori1, norm1, 'MU_SGS_AVG')

	Rhoinf = scipy.integrate.simps(Rho[:], y_coo)/2
	Muinf  = scipy.integrate.simps(Mu[:], y_coo)/2
	Tinf  = scipy.integrate.simps(T[:], y_coo)/2
	a = np.sqrt(gamma*R*Tinf)

	# velocity fields and bulk velocity calculation
	[y_coo, Urms] = util.getVarDataAtPositionNoWall(Slice_out, ori1, norm1, 'U_RMS')
	[y_coo, Urey] = util.getVarDataAtPositionNoWall(Slice_out, ori1, norm1, 'U_REY')
	[y_coo, U]    = util.getVarDataAtPositionNoWall(Slice_out, ori1, norm1, 'U_AVG')

	"############################  NORMALIZATION  #####################################"

	Ut    = np.sum(U*norm1, axis=1) ; Ut[-1]=0 ; Ut[0]=0 # tang. velocity (more general than the 'along x' one)

	Uinf = scipy.integrate.simps(U[:,0], y_coo)/2 ; print 'Uinf classica : ', Uinf
	Ma = Uinf/a ; print 'calculated Mach number (%s) : Ma = '%model, Ma


	# wall density calculation to get tau_wall and normalize in wall unit
	ori2  = [3., -0.99999, 0.]
	norm2 = [0,1,0]
	[x_coo, Rhow] = toca.getVarDataAtWall(Slice_out, ori2, norm2, 'RHO')
	Rhowall       = scipy.integrate.simps(Rhow[:], x_coo)/6
	# print 'Rhowall = ', Rhowall, 'and the first Rho at y=0 : ', Rho[0] # to compare the density at the wall and the bulk one

	# wall shear stress to calculate the Reynolds wall number
	Tw = 0.5*(Ut[5]/y_coo[5]+Ut[-6]/(2.-y_coo[-6]))*Muinf # averaged to be more precise
	Utau = np.sqrt(Tw/Rhowall)
	# Cf = 2*Tw/(Rhoinf*Uinf**2) 			# useless in this case
	# errorCf = (Cf-0.00667)/0.00667*100 	# useless in this case

	""" NOTE : Retau and 1/Lw (reciprocal of the friction length) are hereafter confused 
	because of the channel half height h = 1, but they are not in general the same!!! """

	Retau = Utau*Rhowall/Muinf
	yp      = y_coo*Retau

	Up      = Ut/Utau 		 # wall units velocity
	Mu_sgsp = Mu_sgs.T/Muinf # normalized SGS viscosity
	Urms    = Urms.T/Utau 
	Urey    = Urey.T/Utau**2 # normalized with the square of the friction velocity

	#calculation of the grid size, which is not reliable. It is better to use the ruler on Paraview.
	"""
	"############################  GRID DIMENSIONS  ################################"
	# print 'Normalized grid dimensions : '
	# ori3  = [3,-0.99,0]
	# norm3 = [0,1,0]
	# mean_planeZ = fct.getSlice(out_f, ori1, norm1)
	# [z_coo, Rhowz] = toca.getVarDataAtWall_GridZ(mean_planeZ, ori3, norm3, 'RHO')
	# mean_planeX = fct.getSlice(out_f, ori, norm)
	# [x_coo, Rhowz] = toca.getVarDataAtWall(mean_planeX, ori3, norm2, 'RHO')

	# # print len(x_coo), len(z_coo)

	# # plt.figure()
	# # plt.plot(x_coo,z_coo)
	# # plt.show()

	# print 'deltaYp', np.ediff1d(y_coo)*Retau, '\n deltaXp', np.ediff1d(x_coo)*Retau, '\n deltaZp', np.ediff1d(z_coo)*Retau

	# a test to verify the validity of the linear rang near the wall
	"############################  BEGIN OF THE SLOPE TEST  #####################################"
	fig = plt.figure(figsize=(10, 6))
	ax = fig.add_subplot(111)
	ax.plot(yp,dudy,'b-.',marker='v',label=r'$slopetest$')#,posd[:,0],Cflog*1E3,'k--',pos[:,0],Cfpower*1E3,'b--')

	# ax.plot(yp_VRE[1:],(yp_VRE[1:])**3,ls='-.',color='0.5')
	ax.set_xlabel(r'$y^+$')
	ax.set_ylabel(r'$(du/dy)$')
	ax.set_xlim([min(y_coo),7]) # max(y_coo)])
	# ax.set_ylim([0,30])
	# ax.set_xscale('log')
	# ax.set_yscale('log')
	plt.legend(loc='best')
	plt.grid(which='both', axis='both',color='darkgrey')

	fig.savefig(resdir+'/slope_%s.eps'%model)

	fig = plt.figure(figsize=(10, 6))
	ax = fig.add_subplot(111)
	ax.plot(yp,Ut,'b-.',marker='v',label=r'$data$')
	ax.plot(yp,Tw/Muinf*y_coo,'r-.',label=r'$regression$')
	# ax.plot(yp_VRE[1:],(yp_VRE[1:])**3,ls='-.',color='0.5')
	ax.set_xlabel(r'$y^+$')
	ax.set_ylabel(r'$(U_t)$')
	ax.set_xlim([min(yp),5]) # max(y_coo)])
	ax.set_ylim([0,Ut[15]])
	# ax.set_xscale('log')
	# ax.set_yscale('log')
	plt.legend(loc='best')
	plt.grid(which='both', axis='both',color='darkgrey')

	fig.savefig(resdir+'/linear_range_%s.eps'%model)


	print 'check the first y = ', y_coo[:15]
	print 'check the velocity graient du/dy = ', dudy[:15]
	print 'check the first y+ = ', yp[:15]
	"############################  END OF THE SLOPE TEST  #####################################"
	"""

	# print of some of the results obtained
	yp_slice = np.array([5, 10, 98, 392]) # where we want to calculate the spectra (depending on th DNS data avalaibility)
	y_coo_slice = yp_slice / Retau - np.ones(4)
	y_coo_slice_sym = -y_coo_slice

	print '--> calculated Bulk density (%s): Rhoinf = '%model, Rhoinf
	print '--> calculated Bulk viscosity (%s): Muinf = '%model, Muinf
	print '--> calculated Bulk velocity (%s): Uinf = '%model, Uinf
	print ''
	print '--> calculated wall density (%s): Rhowall = '%model, Rhowall
	print '--> calculated wall velocity (%s): u_tau = '%model, Utau,
	# print '--> calculated viscosity: nu= ', Muinf/Rhoinf
	print '--> calculated wall Reynolds number (%s): Re_tau = '%model, Retau
	print '--> calculated friction length (%s): l_w = '%model, 1/Retau
	print '--> calculated bulk Reynolds number (%s): Re_b = '%model, Retau*Uinf/Utau
	print '--> calculated friction coefficient (%s): Cf = '%model, Cf, ' with an error respect to DNS of : ', errorCf, '%'
	print ''
	# print of the altitude where to slice, because the plane files are not 2D
	# print '--> calculated positions where to slice (%s): y = '%model, np.array_str(y_coo_slice)
	# print '--> the symmetric ones to average these (%s): y = '%model, np.array_str(y_coo_slice_sym)

	"############################  SAVE DATA  #####################################"

	np.savetxt(resdir+'/Profiles_%s.dat'%model,np.vstack([y_coo,yp,Up,Urms[0],Urms[1],Urms[2],Urey[2],Mu_sgsp]).T,fmt="%-13e",
        		header='   y             y+          u_mean+       u_rms+        v_rms+        w_rms+         -uv+          mu_sgs+', footer='', comments='# ')

	# np.savetxt(resdir+'/wall_distance_%s.dat'%model,np.vstack([y_coo,yp]).T,fmt="%-13e",
 #        		header='   y          y+', footer='', comments='# ')

	return [yp, Up, Mu_sgsp, y_coo, Ut, Urms, Urey, Uinf, Utau, Retau]

"############################  PROBES  #####################################"

# reading of the probes data

def	getProbes(workdir,model,Uinf):
	probe  ='/%s_cc/probes/'%model

	loc_probes, nb_probes = toca.read_probes_pos(workdir, probe)

	UX, time_UX = toca.read_probes_data(workdir,probe,'U_AVG-X',start=0)
	flowth_UX = Uinf*time_UX/6 #(Xmax-Xmin) to be set in this form

	MU, time_MU = toca.read_probes_data(workdir,probe,'MU_SGS_AVG',start=0)
	flowth_MU = Uinf*time_MU/6 #(Xmax-Xmin) to be set in this form

	return [flowth_UX, UX, flowth_MU, MU, nb_probes]


"############################  2D PLANES .vtu  #####################################"

# calculation of the spectra to be PP, this function uses instantaneous velocity fields from the .vtu file

def getSpectra(workdir, model):

	workdir = workdir + '/%s_cc/vtu/planes/'%model

	# this function calculates the spectra at the 4 different position, making averages of 
	# symmetric planes thaks to the homogeneous directions

	ID = ['1', '2', '3', '4']
	sliceA = [-0.98877096, -0.97754192, -0.77991079, -0.11964314]
	sliceB = [0.98877096, 0.97754192, 0.77991079, 0.11964314]
	
	coordX = [] ; coordZ = []  
	spectraX = [] ; spectraZ = []
	
	for plane in range(4):
		idp = ID[plane]

		resdir    = workdir + '/../../../PP_results/spectra_plot/plane_%s'%idp
		check_dir = util.fct_dir_exist(resdir)           
		if check_dir == False:
			# os.mkdir(resdir)
			os.makedirs(resdir)

		file_d = [] # ; cnt = 0 # what is that cnt ???...
		for name in os.listdir(workdir): # list entries in workdir in arbitrary order.
		    if os.path.isfile(os.path.join(workdir, name)): # it should not be an "is True:" ?
				if name.startswith("cc.%s_Plane%s"%(model,idp)) is True: 
					if name.endswith(".vtu") is True:
						file_d.append(name) # file_d is a list of planar .vtu files

		print '--> number of planes %sA and %sB founded: '%(idp,idp),  len(file_d) # number of files to be Post Processed? 

		norm = [0,0,1]
		ori = [0,0,1.5]

		"###############################################################################"
	    
		modX = [] ; modZ = []
		specX = [] ; specZ = []

		# y_tick = np.arange(0, 1.25, 0.25)
		# y_label = [r"$0$", r"$\frac{\pi}{4}$", r"$\frac{\pi}{2}$", r"$\frac{3\pi}{4}$", r"$\pi$"]

		# x_tick = np.arange(0, 2.5, 0.5)
		# x_label = [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"]

		time1 = time.time() ; cnt = 0

		for comp in range(3):

			SpanSpectrum = [] ; SpanK = []
			StreamSpectrum = [] ; StreamK = []

			for t in range(len(file_d)): # loop over all the planes (at y=fix) to be PP and then averaged

				out_f = fct.getOutputVTKwithPointDataFromFile(workdir+file_d[t])

				if file_d[t].startswith("cc.%s_Plane%sA"%(model,idp)) is True: # "A", 13) is True: 
					out_f = fct.getSlice(out_f, [0, sliceA[plane], 0], [0, 1, 0])
					[coord, U] = fct.getArrayFromPointData(out_f, 'U')

				elif file_d[t].startswith("cc.%s_Plane%sB"%(model,idp)) is True: # "B", 13) is True:
					out_f = fct.getSlice(out_f, [0, sliceB[plane], 0], [0, 1, 0])
					[coord, U] = fct.getArrayFromPointData(out_f, 'U')

				"# to be COMMENTED when the symmetric planes are available <---------------------  #"
				# out_f = fct.getSlice(out_f, [0, sliceA[plane], 0], [0, 1, 0])
				# [coord, U] = fct.getArrayFromPointData(out_f, 'U')
				"# to be COMMENTED when the symmetric planes are available <---------------------  #"


				x = set(coord.T[0])
				nx = int(6./(1*np.ediff1d(np.sort(coord.T[0])).max()))+1
				# print nx
				z = set(coord.T[2])
				nz = int(3./(1*np.ediff1d(np.sort(coord.T[2])).max()))+1
				# print nz
				line_x = np.linspace(0.,6.,nx)
				line_z = np.linspace(-1.5,1.5,nz)
				grid_x, grid_z = np.meshgrid(line_x,line_z) # coord of the 2D domain: array of x for each z and viceversa
				grid_data = griddata((coord[:,0],coord[:,2]), U[:,comp], (grid_x, grid_z), method='cubic') 

				index = []
				for i in range(nz):
					if np.isnan(grid_data[i]).any():
						index.append(i)

				grid_data = np.delete(grid_data, index, axis=0)
				grid_x = np.delete(grid_x, index, axis=0)
				grid_z = np.delete(grid_z, index, axis=0)

				nz = nz - len(index)

				""" # better if commented if not the calculation is much more slow # """
				# fig = plt.figure(figsize=(14,14))
				# ax1 = fig.add_subplot(211)
				# ax1.imshow(grid_data, extent=(0,6,-1.5,1.5), origin='lower')
				# ax1.plot(grid_x, grid_z, 'rs', ms=2)
				# ax1.plot(coord[:,0], coord[:,2], 'ko', ms=2)
				# ax1.set_title(r'$Original\/data$')
				# ax1.set_xlabel(r'$X$')
				# ax1.set_ylabel(r'$Z$')
				# # ax1.set_yticks(y_tick*np.pi)
				# # ax1.set_yticklabels(y_label, fontsize=18)
				# # ax1.set_xticks(x_tick*np.pi)
				# # ax1.set_xticklabels(x_label, fontsize=18)
				# ax1.set_ylim(min(z),max(z))
				# ax1.set_xlim(min(x),max(x))
				# ax2 = fig.add_subplot(212)
				# ax2.imshow(grid_data, extent=(0,6,-1.5,1.5), origin='lower')
				# ax2.set_title(r'$Interpolated\/-\/Cubic$')
				# ax2.set_xlabel(r'$X$')
				# ax2.set_ylabel(r'$Z$')
				# # ax2.set_yticks(y_tick*np.pi)
				# # ax2.set_yticklabels(y_label, fontsize=18)
				# # ax2.set_xticks(x_tick*np.pi)
				# # ax2.set_xticklabels(x_label, fontsize=18)
				# # plt.gcf().set_size_inches(6, 6)
				# # plt.show()
				# fig.savefig(resdir+'/U%d_plane%s_%s.eps'%(comp,idp,model))

				"##########################  SPECTRA  #################################"

				[KX, FFTX, KZ, FFTZ] = toca.computeSpectrum(nx,nz,grid_x,grid_z,grid_data)


				SpanK.append(np.mean(KZ,axis=0))
				SpanSpectrum.append(np.mean(FFTZ,axis=0))
				
				StreamK.append(np.mean(KX,axis=0))
				StreamSpectrum.append(np.mean(FFTX,axis=0))

				txt_min,seconds = util.time_stat(time1)
				percent = (cnt+1)*100/((len(file_d)-1)*3) # (len(file_d)-1) if we don't average erase -1
				cmd = '    --> Processing: ' + str(percent) + '%' + ' (elapsed time: %s %2.2f seconds)'%(txt_min,seconds)
				util.fct_progress_bar(cmd,percent/2,50)

				cnt += 1

			modX.append(np.mean(StreamK,axis=0))
			specX.append(np.mean(StreamSpectrum,axis=0))

			modZ.append(np.mean(SpanK,axis=0))
			specZ.append(np.mean(SpanSpectrum,axis=0))

		modX = np.array(modX).real ; modZ = np.array(modZ).real
		specX = np.array(specX).real ; specZ = np.array(specZ).real

		"############################  SAVE DATA  #####################################"

		np.savetxt(resdir+'/SpanSpectrum_plane%s_%s.dat'%(idp,model),np.vstack([modZ[0],specZ[0],specZ[1],specZ[2]]).T,fmt="%-13e",
		            header='    k            UU            VV            WW', footer='', comments='# ')
		np.savetxt(resdir+'/StreamSpectrum_plane%s_%s.dat'%(idp,model),np.vstack([modX[0],specX[0],specX[1],specX[2]]).T,fmt="%-13e",
		            header='    k            UU            VV            WW', footer='', comments='# ')

		coordX.append(modX) ; coordZ.append(modZ)
		spectraX.append(specX) ; spectraZ.append(specZ)

	return [coordX, spectraX, coordZ, spectraZ]
