# -*- coding: utf-8 -*-
"""
plotPP.py
"""

import vtk
import draft as fct
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import string, os
import modUtils as util
import toolCase as toca
import PPchan as ppch 
import scipy.integrate

"###############################################################################"

font = {'size'   : 20}

matplotlib.rc('font', **font)
matplotlib.rcParams['lines.linewidth'] = 3
matplotlib.rcParams['lines.markersize'] = 8
matplotlib.rcParams['axes.titlesize'] = 18
matplotlib.rcParams['axes.labelsize'] = 28

"#########################  PARAMETERS SETTING ################################"

# this is the main directory 
workdir = '/Users/andreagallegati/Documents/CharlesX/WALVREDSM'
# this are the 3D .vtk channel results to be PP
chan_VRE = '/cc.VRE.01.683108-0.vtu' 
chan_DSM = '/cc.DSM.01.697874-0.vtu'
chan_WAL = '/cc.01_WAL.697501-0.vtu'
chan_NON = '/cc.NON.01.474702-0.vtu'
# this is the directory where results will be saved
resdir = workdir + '/PP_results/profiles_plot/'
check_dir = util.fct_dir_exist(resdir)           

if check_dir == False:
	os.makedirs(resdir) # os.makedirs: Recursive directory creation function, like mkdir(), 
						# but makes all intermediate-level directories needed to contain the leaf directory. 

"#####################  GET THE 3 MODEL MEAN TIME X ITER #######################"
# calculation of the time per iter statistics: a file timeperiter.dat is needed (obtained by the run.log file,
# using the line command: cat run.log | grep "time-per-iter" | awk '{print $NF}' > timeperiter.dat)

[minimalTimeXiter_VRE, minimalTimeXiter_DSM, minimalTimeXiter_WAL, minimalTimeXiter_NON], [meanTimeXiter_VRE, meanTimeXiter_DSM, meanTimeXiter_WAL, meanTimeXiter_NON], [StD_VRE, StD_DSM, StD_WAL, StD_NON], [timeXiter_VRE, timeXiter_DSM, timeXiter_WAL, timeXiter_NON] = toca.readTimePerIter(workdir)
print 'minimal time per iter in model : (VRE) = ', minimalTimeXiter_VRE, ' (DSM) = ', minimalTimeXiter_DSM, ' (WAL) = ', minimalTimeXiter_WAL, '(NON) = ', minimalTimeXiter_NON
print 'mean time per iter in model : (VRE) = ', meanTimeXiter_VRE, ' (DSM) = ', meanTimeXiter_DSM, ' (WAL) = ', meanTimeXiter_WAL, '(NON) = ', meanTimeXiter_NON
print 'Standard deviation in model : (VRE) = ', StD_VRE, ' (DSM) = ', StD_DSM, ' (WAL) = ', StD_WAL, '(NON) = ', StD_NON
print'--------------------------------------------------------'
print''

"#########################  PP THE 3 CHANNELS ################################"
# calculation of the main quantities to be plotted in this PP script, for each model

print '------------------------------------ getChannel & getProbes in VREMAN case -----------------------------------'
[yp_VRE, Up_VRE, Mu_sgsp_VRE, y_coo_VRE, Ut_VRE, Urms_VRE, Urey_VRE, Uinf_VRE, Utau_VRE, Retau_VRE] = ppch.getChannel(workdir, chan_VRE, 'VRE')
print''
print '------------------------------- getChannel & getProbes in DYN SMAGORINSKY case -------------------------------'
[yp_DSM, Up_DSM, Mu_sgsp_DSM, y_coo_DSM, Ut_DSM, Urms_DSM, Urey_DSM, Uinf_DSM, Utau_DSM, Retau_DSM] = ppch.getChannel(workdir, chan_DSM, 'DSM')
print''
print '------------------------------------ getChannel & getProbes in WALE case -------------------------------------'
[yp_WAL, Up_WAL, Mu_sgsp_WAL, y_coo_WAL, Ut_WAL, Urms_WAL, Urey_WAL, Uinf_WAL, Utau_WAL, Retau_WAL] = ppch.getChannel(workdir, chan_WAL, 'WAL')
print''
print '------------------------------------ getChannel & getProbes in WALE case -------------------------------------'
[yp_NON, Up_NON, Mu_sgsp_NON, y_coo_NON, Ut_NON, Urms_NON, Urey_NON, Uinf_NON, Utau_NON, Retau_NON] = ppch.getChannel(workdir, chan_NON, 'NON')
print''

"#########################  GET THE DNS PROFILES ##############################"
# reading of the DNS database to be plotted as reference

database = '/chan395.means'
y_DNS, yp_DNS, Up_DNS = toca.readDNSdb_Umean(workdir,database,'Up')
database = '/chan395.reystress'
ReS_DNS = [Urms_DNS, Vrms_DNS, Wrms_DNS, UVrms_DNS] = toca.readDNSdb_ReStress(workdir,database)

"###########################  PLOTS PROFILES ##################################"

######################### plot log-law #########################

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

ax.plot(yp_VRE[1:],Up_VRE[1:],'b-.',label=r'$VRE$') # ,marker='v',markevery=5 #,posd[:,0],Cflog*1E3,'k--',pos[:,0],Cfpower*1E3,'b--')
ax.plot(yp_DSM[1:],Up_DSM[1:],'r--',label=r'$DSM$') # ,marker='o',markevery=5 #,posd[:,0],Cflog*1E3,'k--',pos[:,0],Cfpower*1E3,'b--')
ax.plot(yp_WAL[1:],Up_WAL[1:],'k-',label=r'$WAL$') # ,marker='s',markevery=5 #,posd[:,0],Cflog*1E3,'k--',pos[:,0],Cfpower*1E3,'b--')
ax.plot(yp_DNS,Up_DNS,'k',ls='None',marker='o',mfc='none',mew='2',label=r'$DNS$')#,posd[:,0],Cflog*1E3,'k--',pos[:,0],Cfpower*1E3,'b--')

ax.plot(yp_WAL[1:],yp_WAL[1:],ls='-.',color='0.5')
ax.plot(yp_WAL[1:],(1./0.415)*np.log(yp_WAL[1:])+5,ls='-.',color='0.5')
ax.set_xlabel(r'$y^{+}$')
ax.set_ylabel(r'$U^{+}$')
ax.set_xlim([1,400])
ax.set_ylim([0,25])
ax.set_xscale('log')
plt.legend(loc='best')
plt.grid(which='both', axis='both',color='darkgrey')

fig.savefig(resdir+'/Up.eps')

######################## plot viscosity ########################

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

ax.plot(yp_VRE[1:],Mu_sgsp_VRE[1:],'b-.',label=r'$VRE$') # ,marker='v',markevery=5 #,posd[:,0],Cflog*1E3,'k--',pos[:,0],Cfpower*1E3,'b--')
ax.plot(yp_DSM[1:],Mu_sgsp_DSM[1:],'r--',label=r'$DSM$') # ,marker='o',markevery=5 #,posd[:,0],Cflog*1E3,'k--',pos[:,0],Cfpower*1E3,'b--')
ax.plot(yp_WAL[1:],Mu_sgsp_WAL[1:],'k-',label=r'$WAL$') # ,marker='s',markevery=5 #,posd[:,0],Cflog*1E3,'k--',pos[:,0],Cfpower*1E3,'b--')
ax.plot(yp_WAL[1:],(yp_WAL[1:]**3)/2E4,ls='-.',color='0.5')

# ax.plot(yp_VRE[1:],(yp_VRE[1:])**3,ls='-.',color='0.5')
ax.set_xlabel(r'$y^{+}$')
ax.set_ylabel(r'$\mu_{sgs}/\mu$')
ax.set_xlim([min(yp_DSM[1:]),100])
ax.set_ylim([min(Mu_sgsp_WAL),1.4])#,max(Mu_sgsp_WAL)])
ax.set_xscale('log')
# ax.set_yscale('log')
plt.legend(loc='best')
plt.grid(which='both', axis='both',color='darkgrey')

fig.savefig(resdir+'/Mu_sgs.eps')

######################## plot velocity ########################

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

ax.plot(y_coo_VRE,Ut_VRE/Uinf_VRE,'b-.',label=r'$VRE$') # ,marker='v',markevery=5
ax.plot(y_coo_DSM,Ut_DSM/Uinf_DSM,'r--',label=r'$DSM$') # ,marker='o',markevery=5
ax.plot(y_coo_WAL,Ut_WAL/Uinf_WAL,'k-',label=r'$WAL$') # ,marker='s',markevery=5

ax.set_xlabel(r'$y/h$')
ax.set_ylabel(r'$U/U_b$') # ^{+}$')
ax.set_xlim([0,1])
plt.legend(loc='lower right')
plt.grid(which='both', axis='both',color='darkgrey')

# ax.set_xlim([min(y_coo_VRE),0.00001]) # to be commented
# ax.set_ylim([min(Ut_VRE),0.0001]) # to be commented


fig.savefig(resdir+'/Ut.eps')

#################### plot Re stress tensor ####################

# Four subplots, for the reynolds stress tensor components, the axes array is 1-d
# these plots share the x axes fig, ax = plt.subplots(4, sharex=True, figsize=(10, 20))

for sbf in range(0, 4):
	axarr[sbf].grid(which='both', axis='both',color='darkgrey')
	axarr[sbf].locator_params(axis='y',nbins=5)
	axarr[sbf].set_xlim([min(yp_WAL),max(yp_WAL)/2])

	# plot over y_coo instead of yp
	if sbf != 3:
		axarr[sbf].plot(yp_VRE,Urms_VRE[sbf],'b-.',label=r"$VRE$") # marker='v',markevery=5,
		axarr[sbf].plot(yp_DSM,Urms_DSM[sbf],'r--',label=r"$DSM$") # marker='o',markevery=5,
		axarr[sbf].plot(yp_WAL,Urms_WAL[sbf],'k-',label=r"$WAL$") # marker='s',markevery=5,
		axarr[sbf].plot(yp_DNS,np.sqrt(ReS_DNS[sbf]),'k',ls='None',marker='o',mfc='none',mew='2',label=r'$DNS$')

	elif sbf == 3 :
		axarr[3].plot(yp_VRE,Urey_VRE[2],color='blue',linestyle='-.',label=r"$VRE$") # marker='v',markevery=5,
		axarr[3].plot(yp_DSM,Urey_DSM[2],color='red',linestyle='--',label=r"$DSM$") # marker='o',markevery=5,
		axarr[3].plot(yp_WAL,Urey_WAL[2],color='black',linestyle='-',label=r"$WAL$") # marker='s',markevery=5,
		axarr[3].plot(yp_DNS,ReS_DNS[sbf],'k',ls='None',marker='o',mfc='none',mew='2',label=r'$DNS$')

axarr[0].set_ylabel(r"$<u'\,^2>^{1/2}/u_\tau$")
axarr[1].set_ylabel(r"$<v'\,^2>^{1/2}/u_\tau$")
axarr[2].set_ylabel(r"$<w'\,^2>^{1/2}/u_\tau$")
axarr[3].set_ylabel(r"$<-u'v'>/u_\tau^2$")
axarr[3].set_xlabel(r'$y^+$') # ^{+}$')
axarr[0].legend(loc='best',fontsize=16)
axarr[3].legend(loc='best',fontsize=16)


axarr[1].set_yticks(np.arange(0, 2, 0.5))
axarr[3].set_yticks(np.arange(-1, 1, 0.5))
axarr[3].set_ylim([-1,0])
# [min([min(Urey_VRE[2]), min(Urey_DSM[2]), min(Urey_WAL[2])]),0])


fig.savefig(resdir+'/Re_tens.eps')


"############################  PLOTS PROBES #################################"

# reading of the probes data form the working directory

[flowth_UX_VRE, UX_VRE, flowth_MU_VRE, MU_VRE, nb_probes] = ppch.getProbes(workdir, 'VRE', Uinf_VRE)
[flowth_UX_DSM, UX_DSM, flowth_MU_DSM, MU_DSM, nb_probes] = ppch.getProbes(workdir, 'DSM', Uinf_DSM)
[flowth_UX_WAL, UX_WAL, flowth_MU_WAL, MU_WAL, nb_probes] = ppch.getProbes(workdir, 'WAL', Uinf_WAL)

# probes data of the average velocity along X direction

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111)

for sbf in range(nb_probes):
	ax.plot(flowth_UX_VRE.T,UX_VRE.T[sbf],'b-.')#,marker='v',markevery=100,label=r'$VRE$')
	ax.plot(flowth_UX_DSM.T,UX_DSM.T[sbf],'r--')#,marker='o',markevery=100,label=r'$DSM$')
	ax.plot(flowth_UX_WAL.T,UX_WAL.T[sbf],'k-')#,marker='s',markevery=100,label=r'$DSM$')

ax.set_xlabel(r'$flow \quad through$')
ax.set_ylabel(r'$<U_{x}>$')
# plt.legend(loc='best')
plt.grid(which='both', axis='both',color='darkgrey')

fig.savefig(resdir+'/probes_U_AVG.eps')

# probes data of the average sgs viscosity 

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111)

for sbf in range(nb_probes):
	ax.plot(flowth_MU_VRE.T,MU_VRE.T[sbf],'b-.')#,marker='v',markevery=100,label=r'$VRE$')
	ax.plot(flowth_MU_DSM.T,MU_DSM.T[sbf],'r--')#,marker='o',markevery=100,label=r'$DSM$')
	ax.plot(flowth_MU_WAL.T,MU_WAL.T[sbf],'k-')#,marker='s',markevery=100,label=r'$DSM$')

ax.set_xlabel(r'$flow \quad through$')
ax.set_ylabel(r'$<\mu_{sgs}>$')
# plt.legend(loc='best')
plt.grid(which='both', axis='both',color='darkgrey')

fig.savefig(resdir+'/probes_MUsgs_AVG.eps')

print ''
print '--> hurray!!! check your plots at: %s' %resdir

"###########################  PP THE 12 PLANES #################################"

"""	# the following lines are to be decommented if the spectra are not yet calculated. After the first time this script is 
	# launched the calculation of the spectra is saved in the directory so it is possible to read it the following times. # """

# calculation of the spectra for each different model 
# print '------------------------------------ getSpectra in VREMAN case -----------------------------------'
# [coordX_VRE, spectraX_VRE, coordZ_VRE, spectraZ_VRE] = ppch.getSpectra(workdir, 'VRE')
# print''
# print '------------------------------- getSpectra in DYN SMAGORINSKY case -------------------------------'
# [coordX_DSM, spectraX_DSM, coordZ_DSM, spectraZ_DSM] = ppch.getSpectra(workdir, 'DSM')
# print''
# print '------------------------------- getSpectra in WALE case -------------------------------'
# [coordX_WAL, spectraX_WAL, coordZ_WAL, spectraZ_WAL] = ppch.getSpectra(workdir, 'WAL')

# reading of the calculated spectra 
print '------------------------------------ readSpectraDatas in VREMAN case -----------------------------------'
[coordX_VRE, spectraX_VRE, coordZ_VRE, spectraZ_VRE] = toca.readSpectraDatas(workdir, 'VRE')
print''
print '------------------------------- readSpectraDatas in DYN SMAGORINSKY case -------------------------------'
[coordX_DSM, spectraX_DSM, coordZ_DSM, spectraZ_DSM] = toca.readSpectraDatas(workdir, 'DSM')
print''
print '------------------------------- readSpectraDatas in WALE case -------------------------------'
[coordX_WAL, spectraX_WAL, coordZ_WAL, spectraZ_WAL] = toca.readSpectraDatas(workdir, 'WAL')

"#########################  GET THE DNS SPECTRA ################################"
# reading of the DNS database to be plotted as reference

Retau_DNS = 392.24
database_X = '/chan395.xspec.'
database_Z = '/chan395.zspec.'

kDNS_X, SpectraDNS_X = toca.readDNSdb_Spectra(workdir,database_X)
kDNS_Z, SpectraDNS_Z = toca.readDNSdb_Spectra(workdir,database_Z)

"###########################  PLOTS SPECTRA ####################################"

ID = ['1', '2', '3', '4'] # counter for each plane
yPlab = [5, 10, 98, 392]  # wall coordinates for each plane

for plane in range(4):
	idp = ID[plane]

	resdir    = workdir + '/PP_results/spectra_plot/plane_%s'%idp
	check_dir = util.fct_dir_exist(resdir)           
	if check_dir == False:
		os.makedirs(resdir)

	fig = plt.figure(figsize=(10, 10))
	ax = fig.add_subplot(111)


	"########################### BEGIN NORMALIZATIONS ####################################"
	# following lines to be decommented if you want to plot spectra over wavelength instead of wavenumber
	spectraZ_VRE[plane][:,1:] = spectraZ_VRE[plane][:,1:]/(Utau_VRE**2)/31
	# coordZ_VRE[plane][:,1:] = np.divide(2*np.pi*Retau_VRE, coordZ_VRE[plane][:,1:])
	spectraZ_DSM[plane][:,1:] = spectraZ_DSM[plane][:,1:]/(Utau_DSM**2)/31
	# coordZ_DSM[plane][:,1:] = np.divide(2*np.pi*Retau_DSM, coordZ_DSM[plane][:,1:])
	spectraZ_WAL[plane][:,1:] = spectraZ_WAL[plane][:,1:]/(Utau_WAL**2)/31
	# coordZ_WAL[plane][:,1:] = np.divide(2*np.pi*Retau_WAL, coordZ_WAL[plane][:,1:])
	SpectraDNS_Z[plane][:,1:] = SpectraDNS_Z[plane][:,1:]
	# kDNS_Z[plane][1:] = np.divide(2*np.pi*Retau_DNS, kDNS_Z[plane][1:])	

	# check to compare the misalignment of results and DNS reference
	# print 'comparison on the plane%s LES/DNS, UU: '%idp, spectraZ_VRE[plane][0,1]/SpectraDNS_Z[plane][0,1]
	# print 'comparison on the plane%s LES/DNS, VV: '%idp, spectraZ_VRE[plane][1,1]/SpectraDNS_Z[plane][1,1]
	# print 'comparison on the plane%s LES/DNS, WW: '%idp, spectraZ_VRE[plane][2,1]/SpectraDNS_Z[plane][2,1]
	"########################### END NORMALIZATIONS ####################################"

	# plot span-wise spectra
	ax.plot(coordZ_VRE[plane][0,1:], spectraZ_VRE[plane][0,1:], 'r-.')#,label=r'$UU_{VRE}$') # ,marker='v',markevery=5
	ax.plot(coordZ_VRE[plane][1,1:], spectraZ_VRE[plane][1,1:], 'b-.')#,label=r'$VV_{VRE}$') # ,marker='v',markevery=5
	ax.plot(coordZ_VRE[plane][2,1:], spectraZ_VRE[plane][2,1:], 'k-.',label=r'$VRE$')#,label=r'$WW_{VRE}$') # ,marker='v',markevery=5
	ax.plot(coordZ_DSM[plane][0,1:], spectraZ_DSM[plane][0,1:], 'r--')#,label=r'$UU_{DSM}$') # ,marker='o',markevery=5
	ax.plot(coordZ_DSM[plane][1,1:], spectraZ_DSM[plane][1,1:], 'b--')#,label=r'$VV_{DSM}$') # ,marker='o',markevery=5
	ax.plot(coordZ_DSM[plane][2,1:], spectraZ_DSM[plane][2,1:], 'k--',label=r'$DSM$')#,label=r'$WW_{DSM}$') # ,marker='o',markevery=5
	ax.plot(coordZ_WAL[plane][0,1:], spectraZ_WAL[plane][0,1:], 'r-')#,label=r'$UU_{WAL}$') # ,marker='s',markevery=5
	ax.plot(coordZ_WAL[plane][1,1:], spectraZ_WAL[plane][1,1:], 'b-')#,label=r'$VV_{WAL}$') # ,marker='s',markevery=5
	ax.plot(coordZ_WAL[plane][2,1:], spectraZ_WAL[plane][2,1:], 'k-',label=r'$WAL$')#,label=r'$WW_{WAL}$') # ,marker='s',markevery=5

	ax.plot(kDNS_Z[plane][1:], SpectraDNS_Z[plane][0,1:],'k',ls='None',marker='o',mfc='none',mew='2',label=r'$DNS$')
	ax.plot(kDNS_Z[plane][1:], SpectraDNS_Z[plane][1,1:],'k',ls='None',marker='o',mfc='none',mew='2')#,label=r'$DNS$')
	ax.plot(kDNS_Z[plane][1:], SpectraDNS_Z[plane][2,1:],'k',ls='None',marker='o',mfc='none',mew='2')#,label=r'$DNS$')

	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlabel(r'$k_z$') # ax.set_xlabel(r'$\lambda_z^+$') if plot over wavelength
	ax.set_ylabel(r'$E/(u_\tau^2)$') # ax.set_ylabel(r'$k_zE/u_\tau^2$') if normalized in this manner
	ax.set_title(r'Spanwise spectra at $y^+=%d$'%yPlab[plane])
	plt.grid(which='both', axis='both',color='darkgrey')

	if plane==0: ax.set_xlim([9,2000])
	if plane==1: ax.set_xlim([9,2000])
	if plane==3: 
		ax.plot(kDNS_Z[plane][1:],(kDNS_Z[plane][1:]**(-5/3))*1E0,ls='-.',color='0.5')
		ax.set_xlim([2,300])

	plt.legend(loc='lower center')
	fig.savefig(resdir+'/SpanSpectrum_%s.eps'%idp)


	fig = plt.figure(figsize=(10, 10))
	ax = fig.add_subplot(111)
	"########################### BEGIN NORMALIZATIONS ####################################"
	# following lines to be decommented if you want to plot spectra over wavelength instead of wavenumber
	spectraX_VRE[plane][:,1:] = spectraX_VRE[plane][:,1:]/(Utau_VRE**2)/31
	# coordX_VRE[plane][:,1:] = np.divide(2*np.pi*Retau_VRE, coordX_VRE[plane][:,1:])
	spectraX_DSM[plane][:,1:] = spectraX_DSM[plane][:,1:]/(Utau_DSM**2)/31
	# coordX_DSM[plane][:,1:] = np.divide(2*np.pi*Retau_DSM, coordX_DSM[plane][:,1:])
	spectraX_WAL[plane][:,1:] = spectraX_WAL[plane][:,1:]/(Utau_WAL**2)/31
	# coordX_WAL[plane][:,1:] = np.divide(2*np.pi*Retau_WAL, coordX_WAL[plane][:,1:])
	SpectraDNS_X[plane][:,1:] = SpectraDNS_X[plane][:,1:]
	# kDNS_X[plane][1:] = np.divide(2*np.pi*Retau_DNS, kDNS_X[plane][1:])	
	"########################### END NORMALIZATIONS ####################################"
	
	# plot stream-wise spectra
	ax.plot(coordX_VRE[plane][0,1:], spectraX_VRE[plane][0,1:], 'r-.')#,label=r'$UU_{VRE}$') # ,marker='v',markevery=5
	ax.plot(coordX_VRE[plane][1,1:], spectraX_VRE[plane][1,1:], 'b-.')#,label=r'$VV_{VRE}$') # ,marker='v',markevery=5
	ax.plot(coordX_VRE[plane][2,1:], spectraX_VRE[plane][2,1:], 'k-.',label=r'$VRE$')#,label=r'$WW_{VRE}$') # ,marker='v',markevery=5
	ax.plot(coordX_DSM[plane][0,1:], spectraX_DSM[plane][0,1:], 'r--')#,label=r'$UU_{DSM}$') # ,marker='o',markevery=5
	ax.plot(coordX_DSM[plane][1,1:], spectraX_DSM[plane][1,1:], 'b--')#,label=r'$VV_{DSM}$') # ,marker='o',markevery=5
	ax.plot(coordX_DSM[plane][2,1:], spectraX_DSM[plane][2,1:], 'k--',label=r'$DSM$')#,label=r'$WW_{DSM}$') # ,marker='o',markevery=5
	ax.plot(coordX_WAL[plane][0,1:], spectraX_WAL[plane][0,1:], 'r-')#,label=r'$UU_{WAL}$') # ,marker='s',markevery=5
	ax.plot(coordX_WAL[plane][1,1:], spectraX_WAL[plane][1,1:], 'b-')#,label=r'$VV_{WAL}$') # ,marker='s',markevery=5
	ax.plot(coordX_WAL[plane][2,1:], spectraX_WAL[plane][2,1:], 'k-',label=r'$WAL$')#,label=r'$WW_{WAL}$') # ,marker='s',markevery=5

	ax.plot(kDNS_X[plane][1:], SpectraDNS_X[plane][0,1:],'k',ls='None',marker='o',mfc='none',mew='2',label=r'$DNS$')
	ax.plot(kDNS_X[plane][1:], SpectraDNS_X[plane][1,1:],'k',ls='None',marker='o',mfc='none',mew='2',label=r'$DNS$')
	ax.plot(kDNS_X[plane][1:], SpectraDNS_X[plane][2,1:],'k',ls='None',marker='o',mfc='none',mew='2',label=r'$DNS$')

	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlabel(r'$k_x$') # ax.set_xlabel(r'$\lambda_x^+$')
	ax.set_ylabel(r'$E/(h\,u_\tau^2)$') # ax.set_ylabel(r'$k_xE/u_\tau^2$')
	ax.set_title(r'$Streamwise$')
	plt.grid(which='both', axis='both',color='darkgrey')
	plt.legend(loc='upper right')
	fig.savefig(resdir+'/StreamSpectrum_%s.eps'%idp)

	print ''
	print '--> hurray!!! check your plots at: %s'%resdir

# to be decommented if you want a plot of the time per iter statistics
"""
######################## plot time per iter ########################

fig = plt.figure(figsize=(24, 8))
ax = fig.add_subplot(111)

iter_VRE = len(timeXiter_VRE)
iter_DSM = len(timeXiter_DSM)
iter_WAL = len(timeXiter_WAL)

ax.plot(np.arange(iter_VRE),timeXiter_VRE,'b-.',ls='None',marker='o',mfc='b',markevery=750,label=r'$VRE$') # ,ms=0.1) # ,marker='v',markevery=5
ax.plot(np.arange(iter_DSM),timeXiter_DSM,'r--',ls='None',marker='o',mfc='r',markevery=750,label=r'$DSM$') # ,ms=0.1) # ,marker='o',markevery=5
ax.plot(np.arange(iter_WAL),timeXiter_WAL,'k-',ls='None',marker='o',mfc='k',markevery=750,label=r'$WAL$') # ,ms=0.1) # ,marker='s',markevery=5

ax.set_xlabel(r'$\#\, iterations$')
ax.set_ylabel(r'$\Delta t_{iter}\,[\mu s]$') # ^{+}$')
ax.set_xlim([0,max(iter_VRE, iter_DSM, iter_WAL)])
ax.set_ylim([6,8])
plt.legend(loc='best')
plt.grid(which='both', axis='both',color='darkgrey')

# ax.set_xlim([min(y_coo_VRE),0.00001]) # to be commented
# ax.set_ylim([min(Ut_VRE),0.0001]) # to be commented


fig.savefig(resdir+'/../timeXiter.eps')
"""

# to be decommented if you want to verify the validity of the Parseval Thm.
"""
"############################# recovering the PARSEVAL THM for the VRE data #############################"

URMS_5yPlus = (Urms_VRE[0].T[11]+Urms_VRE[0].T[10])/2 ; print 'URMS_5yPlus = ', URMS_5yPlus
VRMS_5yPlus = (Urms_VRE[1].T[11]+Urms_VRE[1].T[10])/2 ; print 'VRMS_5yPlus = ', VRMS_5yPlus
WRMS_5yPlus = (Urms_VRE[2].T[11]+Urms_VRE[2].T[10])/2 ; print 'WRMS_5yPlus = ', WRMS_5yPlus
print 'u_rmsPlus energy : ', (URMS_5yPlus**2)
print 'v_rmsPlus energy : ', (VRMS_5yPlus**2)
print 'w_rmsPlus energy : ', (WRMS_5yPlus**2)

Kuu = np.real((spectraZ_VRE[0][0,0] + np.sum(spectraZ_VRE[0][0,1:])))
Kvv = np.real((spectraZ_VRE[0][1,0] + np.sum(spectraZ_VRE[0][1,1:])))
Kww = np.real((spectraZ_VRE[0][2,0] + np.sum(spectraZ_VRE[0][2,1:])))

print 'u_rms spectra energy : ', Kuu, '\n v_rms spectra energy : ', Kvv, '\n w_rms spectra energy : ', Kww
print 'LES -----> parseval theorem RATIO integ.Euu/Urms^2 = ', Kuu/((URMS_5yPlus**2))
print 'LES -----> parseval theorem RATIO integ.Evv/Vrms^2 = ', Kvv/((VRMS_5yPlus**2))
print 'LES -----> parseval theorem RATIO integ.Eww/Wrms^2 = ', Kww/((WRMS_5yPlus**2))

print 'ratio between myEuu/DNS_Euu first components', spectraZ_VRE[0][0,1]/SpectraDNS_Z[plane][0,1]
print 'ratio between myEvv/DNS_Evv first components', spectraZ_VRE[0][1,1]/SpectraDNS_Z[plane][1,1]
print 'ratio between myEww/DNS_Eww first components', spectraZ_VRE[0][2,1]/SpectraDNS_Z[plane][2,1]



"############################# recovering the PARSEVAL THM for the DNS data #############################"

URMS_5yP_DNS = 1.8 ; print 'URMS_5yP_DNS = ', URMS_5yP_DNS
VRMS_5yP_DNS = 0.144 ; print 'VRMS_5yP_DNS = ', VRMS_5yP_DNS
WRMS_5yP_DNS = 0.727 ; print 'WRMS_5yP_DNS = ', WRMS_5yP_DNS
print 'u_rmsPlus energy : ', (URMS_5yP_DNS**2)
print 'v_rmsPlus energy : ', (VRMS_5yP_DNS**2)
print 'w_rmsPlus energy : ', (WRMS_5yP_DNS**2)

Kuu = np.real((SpectraDNS_Z[0][0,0] + np.sum(SpectraDNS_Z[0][0,1:])))
Kvv = np.real((SpectraDNS_Z[0][1,0] + np.sum(SpectraDNS_Z[0][1,1:])))
Kww = np.real((SpectraDNS_Z[0][2,0] + np.sum(SpectraDNS_Z[0][2,1:])))



print 'u_rms spectra energy : ', Kuu, '\n v_rms spectra energy : ', Kvv, '\n w_rms spectra energy : ', Kww
print 'DNS -----> parseval theorem RATIO integ.Euu/Urms^2 = ', Kuu/((URMS_5yP_DNS**2))
print 'DNS -----> parseval theorem RATIO integ.Evv/Vrms^2 = ', Kvv/((VRMS_5yP_DNS**2))
print 'DNS -----> parseval theorem RATIO integ.Eww/Wrms^2 = ', Kww/((WRMS_5yP_DNS**2))
"""

























