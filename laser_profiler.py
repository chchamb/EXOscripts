# Module to do laser profile analysis 
date=20170120
#-----------------------------------------------------------------------------------------

# Script to read and average files in a directory (using fringe_temp.vi)
def pow_ave(filedir,dx,units,**kwargs):

	# Import numpy and asciitable
	import asciitable
	import numpy as np
	import os
	
	# Keywords
	zvals=kwargs.get('zvals','default')
	
	# Get filenames from filedir
	filenames=os.listdir(filedir)
	filenames.sort()
	avpow=[]
	after_shot=False
	for filename in filenames:
		if filename=='ave.txt':
			os.remove(filedir+'/'+filename)
			continue
		if filename[-4:]!='.txt':
			continue
		datafile= open(filedir+'/'+filename,'r')
		if zvals=='default':
			z=filename[:3]
		if zvals!='default':
			z=1
		gnarpow=[]	
		for line in datafile:
			line=line.strip()
			columns=line.split()
			gnarpow.append(float(columns[2]))
		if filename[-5:]=='a.txt':
			after_shot=True
			after=np.mean(gnarpow)
		if filename[-5:]!='a.txt':
			avpow.append(np.mean(gnarpow))
	datafile.close()
	pows=np.array(avpow)
	# Make position array
	pos=np.zeros(np.shape(pows),dtype=float)
	if units=='standard':
		conv=25.4
	if units=='mm':
		conv=1
	if units=='micron':
		conv=1E-3
	for i in range(1,len(pos)):
		pos[i]=pos[i-1]+(dx*conv)
	# Make file to write data to
	writefile=filedir+'/ave.txt'
	newfile=open(writefile,'a')
	# Convert to numpy array and write to ascii 
	asciitable.write({'Position':pos,'Average Power':pows},newfile,names=['Position','Average Power'])
	newfile.close()
	
	if after_shot==False:
		return pos,pows,z
	if after_shot==True:
		return pos,pows,z,after
	
#-----------------------------------------------------------------------------------------

# Script to read profile data and fit erf from fringe_temp type files

def profile_fit(filedir,dx,**kwargs):

	# Define kwargs
	imsave=kwargs.get('savefile','none')
	fitrange=kwargs.get('fit_range','full')
	
	# Import modules
	import numpy as np
	import matplotlib.pyplot as plt
	import math
	from scipy.special import erf
	from scipy.optimize import curve_fit
	import pylab
	import os
	
	# Define fit function (erf for low to high power)
	def fit_erf(x,a,x_0,w):
		return (a/2.)*(1+erf((np.sqrt(2)*(x-x_0))/w))
		
	# Read profile data
	filenames=os.listdir(filedir)
	filenames.sort()
	beam_width=np.zeros((len(filenames),),dtype=float)
	var=np.zeros((len(filenames),),dtype=float)
	pows=[]
	for filename in filenames:
		# Get the data
		datafile=open(filedir+'/'+filename,'r')
		pows2=[]
		for line in datafile:
			line=line.strip()
			columns=line.split()
			pows2.append(float(columns[2]))
		pows.append(np.mean(pows2))
		datafile.close()
		
	# Make position array and power array (normalized low to high power)
	pos=np.arange(0,len(pows)*dx,dx)*1e-3
	if pows[0]>pows[-1]:
		pows_norm=np.flipud(np.array(pows))/np.amax(np.array(pows))
	else:
		pows_norm=np.array(pows)/np.amax(np.array(pows))
	# Get initial parameter guesses
	param_guess=[1,0,0]
	# Guess for w
	for i in range(len(pows)):
		if pows_norm[i]<=.1:
			w_low=i
			continue
		else:
			break
	for i in range(len(pows)):
		if pows_norm[i]<=.9:
			w_high=i
			continue
		else:
			break
	#param_guess[2]=pos[w_high]-pos[w_low]
	param_guess[2]=w_high-w_low
	# Guess for x_0
	for i in range(len(pows)):
		if pows_norm[i]<=.5:
			param_guess[1]=pos[i]
			continue
		else:
			if pows_norm[i]-.5>=.5-pows_norm[i-1]:
				param_guess[1]=pos[i-1]
				break
			else:
				break
	# Do the fit
	if fitrange!='full':
		param,covar=curve_fit(fit_erf,pos[fitrange[0]:fitrange[1]],pows_norm[fitrange[0]:fitrange[1]],p0=param_guess)
	if fitrange=='full':
		param,covar=curve_fit(fit_erf,pos,pows_norm,p0=param_guess)
	beam_width[filenames.index(filename)]=param[2]
	var[filenames.index(filename)]=covar[2,2]
		
	# Plot the fits
	x=np.linspace(0,np.amax(pos),100)
	fit=fit_erf(x,param[0],param[1],param[2])
	fig=plt.figure('laser_prof')
	plt.clf()
	ax=fig.add_subplot(111)
	plt.plot(pos,pows_norm,'bo',label='Data')
	plt.plot(x,fit,'k--',label='Fit')
	plt.axis([0,1.1*np.amax(pos),0,1.25])
	plt.xlabel('Razor Blade Position (mm)')
	plt.ylabel('Normalized Laser Power')
	plt.legend()
	bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
	fig_text='Spot Size: '+repr(1000*round(param[2],4))+' +/- '+repr(round(1000*np.sqrt(covar[2,2]),2))+' micron'
	#fig_text='Amplitude: '+repr(round(param[0],2))+'\nOffset: '+repr(round(param[1],2))+' mm'+'\nSpot Size: '+repr(1000*round(param[2]),3)+' micron'
	plt.text(.02,.97,fig_text,fontsize=10,bbox=bbox_props,va='top',ha='left',transform=ax.transAxes)
	if imsave!='none':
		plt.savefig(imsave+'beam_profile.png')
	plt.show()
	
	return beam_width
		
#-----------------------------------------------------------------------------------------
	
# Script to do laser profile
def beam_profile(maindir,dx_in,imsave,**kwargs):
	
	# Define Keyword Arguments
	z_units=kwargs.get('z_units','standard')
	data_units=kwargs.get('data_units','standard')
	afters=kwargs.get('after',False)
	param_guess=kwargs.get('p0','find')
	param_guess2=kwargs.get('p0_waist','find')
	zvals=kwargs.get('zvals','default')
	fitlim=kwargs.get('fitlim','full')
	stages=kwargs.get('stage_type','auto')
	
	# Import Modules
	import numpy as np
	import matplotlib.pyplot as plt
	import math
	from scipy.special import erf
	from scipy.optimize import curve_fit
	import pylab
	import os
	import utility as ut
	import gen_reader as gr
	
	# Create directory for imsave
	ut.create_dir(imsave)
	
	# Define fit function (erf for low to high power)
	def fit_erf(x,a,x_0,w):
		return (a/2.)*(1+erf((np.sqrt(2)*(x-x_0))/w))
	
	# Get the files in the directory and remove '.DS_Store' if it exists
	filedirs=os.listdir(maindir)
	filedirs.sort()
	
	# If all dx are the same
	if len(dx_in)==len(filedirs):
		dx=dx_in
	else:
		dx=np.zeros((len(filedirs),),dtype=float)
		dx[:]=dx_in
		
	# Fit the data to get spot size for each z-position	
	beam_width=np.zeros((len(filedirs),),dtype=float)
	weights=np.zeros((len(filedirs),),dtype=float)
	z=np.zeros((len(filedirs),),dtype=int)
	after=np.zeros((len(filedirs),),dtype=int)
	for filedir in filedirs:
		if stages=='manual':
			if afters==True:
				pos,pows,z[filedirs.index(filedir)],after[filedirs.index(filedir)]=pow_ave(maindir+'/'+filedir,dx[filedirs.index(filedir)],data_units,zvals=zvals)
			if afters==False:
				pos,pows,z[filedirs.index(filedir)]=pow_ave(maindir+'/'+filedir,dx[filedirs.index(filedir)],data_units,zvals=zvals)
		if stages=='auto':
			z[filedirs.index(filedir)]=np.float(filedir[-9:-4])
			matrix=gr.reader(maindir+filedir,header=False)
			pos,pows=(.166E-3/25.)*matrix[:,0],matrix[:,1]
		if stages=='profman':
			xstep,z[filedirs.index(filedir)]=filedir[-9:-4],filedir[-29:-24]
			matrix=gr.reader(maindir+filedir,header=False)
			pows=matrix[:,0]
			pos=1E-3*np.arange(0,float(xstep)*len(pows),float(xstep))
			
		# Normalize and reorder the data for low to high power
		pows_norm=pows/np.amax(pows)
		if pows_norm[0]>pows_norm[-1]:
			pows_norm=np.flipud(pows_norm)
		
		# Get initial parameter guesses
		if param_guess=='find':
			param_guess=[1,0,0]
			# Guess for w
			bnds=ut.bound_finder(pows_norm,[.1,.9])
			param_guess[2]=pos[bnds[1]]-pos[bnds[0]]
			for i in range(len(pows)):
				if pows_norm[i]<=.5:
					param_guess[1]=pos[i]
					continue
				else:
					if pows_norm[i]-.5>=.5-pows_norm[i-1]:
						param_guess[1]=pos[i-1]
						break
					else:
						break
	
		# Do the fit
		if fitlim=='full':
			posf=pos
			pows_normf=pows_norm
		if fitlim !='full':
			posf=pos[fitlim[0]:fitlim[1]]
			pows_normf=pows_norm[fitlim[0]:fitlim[1]]
		param,covar=curve_fit(fit_erf,posf,pows_normf,p0=param_guess)
		beam_width[filedirs.index(filedir)]=param[2]
		weights[filedirs.index(filedir)]=1/np.sqrt(covar[2,2])
		
		# Plot Fit
		x=np.linspace(0,pos[-1],100)
		fit=fit_erf(x,param[0],param[1],param[2])
		fig=plt.figure('profile')
		plt.clf()
		ax=fig.add_subplot(111)
		plt.plot(pos,pows_norm,'bo',label='Data')
		plt.plot(x,fit,'k--',label='Fit')
		plt.axis([0,1.1*np.amax(pos),0,1.25])
		plt.xlabel('Razor Blade Position (mm)')
		plt.ylabel('Normalized Laser Power')
		plt.legend()
		bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
		fig_text='Spot Size: '+repr(round(1000.*param[2],2))+' +/- '+repr(round(1000*np.sqrt(covar[2,2]),2))+' micron'
		#fig_text='Amplitude: '+repr(round(param[0],2))+'\nOffset: '+repr(round(param[1],2))+' mm'+'\nSpot Size: '+repr(1000*round(param[2]),3)+' micron'
		plt.text(.02,.97,fig_text,fontsize=10,bbox=bbox_props,va='top',ha='left',transform=ax.transAxes)
		if zvals=='default':
			pylab.savefig(imsave+'/'+repr(z[filedirs.index(filedir)])+'_fit.png')
		if zvals!='default':
			pylab.savefig(imsave+'/'+repr(zvals[filedirs.index(filedir)])+'_fit.png')
		plt.show()
	
	# Fit the beam widths to find beam waist
	# Define the fit function
	wlength=572 # nm
	def fit_waist(z,w_0,f):
		z_0=(math.pi*np.square(w_0))/(wlength*1E-6)
		return w_0*np.sqrt(1+np.square((z-f)/z_0))
	
	# Switch to manual z vals
	if zvals!='default':
		z=np.array(zvals)
	
	# Convert z-positions to array with 0 = lowest position and units are mm
	if z_units=='standard':
		conv=25.4
	if z_units=='mm':
		conv=1
	if z_units=='micron':
		conv=.001
	z=(z-z[0])*conv
	
	# Set initial parameters for beam fit
	if param_guess2=='find':
		# Get initial parameter guesses
		param_guess2=[]
		slope=(beam_width[1]-beam_width[0])/(z[1])
		theta=2*np.arctan(np.absolute(slope))
		param_guess2.append(2*(wlength*1E-6)/(math.pi*theta))
		param_guess2.append(-beam_width[0]/slope)
	
	# Normalize the weights
	weights_n=weights/np.amax(weights)
	
	# Do the fit
	param2,covar2=curve_fit(fit_waist,z,beam_width,p0=param_guess2,sigma=weights_n)
	
	# Plot the fit and data
	z_cont=np.linspace(0,np.amax(z),100)
	fit2=fit_waist(z_cont,param2[0],param2[1])
	fig=plt.figure('Beam Fit')
	plt.clf()
	ax=fig.add_subplot(111)
	plt.plot(z,beam_width,'bo',label='Data')
	plt.plot(z_cont,fit2,'k--',label='Fit')
	plt.axis([0,np.amax(z),0,1.25*np.amax(beam_width)])
	plt.xlabel('z (mm)')
	plt.ylabel('Spot Size (mm)')
	plt.legend()
	plt.errorbar(z,beam_width,xerr=None,yerr=1/weights,fmt=None)
	bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
	fig_text='Beam Waist: '+repr(round(1000*param2[0],2))+' +/- '+repr(round(1000*np.sqrt(covar2[0,0]),2))+' micron \nWaist Pos: '+repr(round(param2[1],2))+' +/- '+repr(round(np.sqrt(covar2[1,1]),2))+' mm'
	plt.text(.02,.97,fig_text,fontsize=10,bbox=bbox_props,va='top',ha='left',transform=ax.transAxes)
	pylab.savefig(imsave+'/beam_fit.png')
	plt.show()
	
#-----------------------------------------------------------------------------------------
	
# Fit the beam widths to find beam waist
def waist_fit(z,beam_width,weight,**kwargs):

	imsave=kwargs.get('savepath',False)
	disp=kwargs.get('print',False)
	
	# Import Modules
	import numpy as np
	import matplotlib.pyplot as plt
	import math
	from scipy.special import erf
	from scipy.optimize import curve_fit
	import utility as ut
	
	# Convert z-positions to array with 0 = lowest position and units are mm
	z=np.array(z)
	z=(z-z[0])*(25.4/1000)

	# Define the fit function
	
	def fit_waist(z,w_0,f):
		z_0=(math.pi*np.square(w_0))/(572*1E-6)
		return w_0*np.sqrt(1+np.square((z-f)/z_0))
	param_guess=[4E-3,np.mean(z)]
	
	# Normalize the weights
	weight=weight/np.amax(weight)
	
	# Do the fit
	param,covar=curve_fit(fit_waist,z,beam_width,p0=param_guess,sigma=weight)
	
	# Plot the fit and data
	z_cont=np.linspace(0,np.amax(z),100)
	if flength=='fit':
		fit=fit_waist(z_cont,param[0],param[1])
	if flength!='fit':
		fit=fit_waist(z_cont,param[0])
	fig=plt.figure('Beam Fit')
	plt.clf()
	ax=fig.add_subplot(111)
	plt.plot(z,beam_width,'bo',label='Data')
	plt.plot(z_cont,fit,'k--',label='Fit')
	plt.axis([-.1,1.1*np.amax(z),0,1.25*np.amax(beam_width)])
	plt.xlabel('z (mm)')
	plt.ylabel('Spot Size (mm)')
	plt.legend()
	bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
	fig_text='Beam Waist: '+ut.rnd(1000*param[0],2)+' +/- '+ut.rnd(1000*np.sqrt(covar[0,0]),2)+' micron'
	plt.text(.02,.97,fig_text,fontsize=10,bbox=bbox_props,va='top',ha='left',transform=ax.transAxes)
	if imsave!=False:
		ut.create_dir(imsave)
		plt.savefig(imsave+'/beam_fit.png')
	if disp==True:
		print('w0: '+ut.rnd(1000*param[0],2)+' um\nz0: '+ut.rnd((math.pi*np.square(1000*param[0]))/(572*1E-3),2)+' um\nf: '+ut.rnd(param[1],2)+' mm')
	plt.show()
	
#-----------------------------------------------------------------------------------------

# Script to fit single profile
def solo_prof(filedir,dx,imsave,**kwargs):	
	
	#Define Keyword Arguments
	data_units=kwargs.get('data_units','standard')
	calc_units=kwargs.get('calc_units','micron')
	
	# Import Modules
	import numpy as np
	import matplotlib.pyplot as plt
	import math
	from scipy.special import erf
	from scipy.optimize import curve_fit
	import pylab
	import os

	# Define fit function (erf for low to high power)
	def fit_erf(x,a,x_0,w):
		return (a/2.)*(1+erf((np.sqrt(2)*(x-x_0))/w))
		
	pos,pows,z,after=pow_ave(filedir,dx,data_units)
	# Normalize and reorder the data for low to high power
	pows_norm=pows/np.amax(pows)
	pows_norm=pows_norm-np.amin(pows_norm)
	if pows_norm[0]>pows_norm[-1]:
		pows_norm=np.flipud(pows_norm)
	# Get initial parameter guesses
	param_guess=[1,0,0]
	for i in range(len(pows)):
		if pows_norm[i]<=.1:
			w_low=i
			continue
		else:
			break
	for i in range(len(pows)):
		if pows_norm[i]<=.9:
			w_high=i
			continue
		else:
			break
	param_guess[2]=pos[w_high]-pos[w_low]
	for i in range(len(pows)):
		if pows_norm[i]<=.5:
			param_guess[1]=pos[i]
			continue
		else:
			if pows_norm[i]-.5>=.5-pows_norm[i-1]:
				param_guess[1]=pos[i-1]
				break
			else:
				break
	# Do the fit
	param,covar=curve_fit(fit_erf,pos,pows_norm,p0=param_guess)
		
	# Plot Fit
	x=np.linspace(0,pos[-1],100)
	fit=fit_erf(x,param[0],param[1],param[2])
	fig=plt.figure(z)
	plt.clf()
	ax=fig.add_subplot(111)
	plt.plot(pos,pows_norm,'bo',label='Data')
	plt.plot(x,fit,'k--',label='Fit')
	plt.axis([0,1.1*np.amax(pos),0,1.25])
	plt.xlabel('Razor Blade Position (mm)')
	plt.ylabel('Normalized Laser Power')
	plt.legend()
	bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
	if calc_units=='micron':
		fig_text='Spot Size: '+repr(round(1000*param[2],1))+' +/- '+repr(round(1000*np.sqrt(covar[2,2]),2))+' micron'
	if calc_units=='mm':
		fig_text='Spot Size: '+repr(round(param[2],3))+' +/- '+repr(round(np.sqrt(covar[2,2]),2))+' mm'
	#fig_text='Amplitude: '+repr(round(param[0],2))+'\nOffset: '+repr(round(param[1],2))+' mm'+'\nSpot Size: '+repr(1000*round(param[2]),3)+' micron'
	plt.text(.02,.97,fig_text,fontsize=10,bbox=bbox_props,va='top',ha='left',transform=ax.transAxes)
	pylab.savefig(imsave)
	plt.show()

#-----------------------------------------------------------------------------------------

# Script to read profile data and fit erf from auto stage

def prof_solo(profnum,z,**kwargs):

	# Define kwargs
	imsave=kwargs.get('savefile','none')
	fitrange=kwargs.get('fit_range','full')
	xstep=kwargs.get('xstep','auto')
	units=kwargs.get('units','microns')
	
	# Import modules
	import numpy as np
	import matplotlib.pyplot as plt
	import math
	from scipy.special import erf
	from scipy.optimize import curve_fit
	import pylab
	import os
	import gen_reader as gr
	import utility as ut
	
	# Define fit function (erf for low to high power)
	def fit_erf(x,a,x_0,w,offset):
		return (a/2.)*(1+erf((np.sqrt(2)*(x-x_0))/w))+offset
		
	# Read in data from file
	if xstep=='auto':
		matrix=gr.reader('/home/chris/anaconda/data/'+str(date)+'/lv/prof'+str(profnum)+'/prof'+str(profnum)+'_amp25_zpos_micron_'+z+'.txt',header=False)
		pos,pows=(.166E-3/25.)*matrix[:,0],matrix[:,1]
	if xstep!='auto':
		matrix=gr.reader('/home/chris/anaconda/data/'+str(date)+'/lv/profman'+str.zfill(profnum,2)+'/profman'+str.zfill(profnum,2)+'_z-pos-'+units+str.zfill(z,5)+'_x-step-'+units+str.zfill(xstep,5)+'.txt',header=False)
		pows=matrix[:,0]
		pos=1E-3*np.arange(0,float(xstep)*len(pows),float(xstep))
		if units=='minches':
			pos=pos*25.4

	# Make position array and power array (normalized low to high power)
	if pows[0]>pows[-1]:
		pows_norm=np.flipud(np.array(pows))/np.amax(np.array(pows))
	else:
		pows_norm=np.array(pows)/np.amax(np.array(pows))
	# Get initial parameter guesses
	param_guess=[1,0,0,0]
	# Guess for w
	bnds=ut.bound_finder(pows_norm,[.1,.9])
	param_guess[2]=pos[bnds[1]]-pos[bnds[0]]
	# Guess for x_0
	for i in range(len(pows)):
		if pows_norm[i]<=.5:
			param_guess[1]=pos[i]
			continue
		else:
			if pows_norm[i]-.5>=.5-pows_norm[i-1]:
				param_guess[1]=pos[i-1]
				break
			else:
				break
	# Do the fit
	if fitrange!='full':
		param,covar=curve_fit(fit_erf,pos[fitrange[0]:fitrange[1]],pows_norm[fitrange[0]:fitrange[1]],p0=param_guess)
	if fitrange=='full':
		param,covar=curve_fit(fit_erf,pos,pows_norm,p0=param_guess)
		
	# Plot the fits
	x=np.linspace(0,np.amax(pos),100)
	fit=fit_erf(x,param[0],param[1],param[2],param[3])
	fig=plt.figure('laser_prof')
	plt.clf()
	ax=fig.add_subplot(111)
	plt.plot(pos,pows_norm,'bo',label='Data')
	plt.plot(x,fit,'k--',label='Fit')
	plt.axis([0,1.1*np.amax(pos),0,1.25])
	plt.xlabel('Razor Blade Position (mm)')
	plt.ylabel('Normalized Laser Power')
	plt.legend()
	bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
	fig_text='Spot Size: '+ut.rnd(1000*param[2],2)+' +/- '+ut.rnd(1000*np.sqrt(covar[2,2]),2)+' micron'
	#fig_text='Amplitude: '+repr(round(param[0],2))+'\nOffset: '+repr(round(param[1],2))+' mm'+'\nSpot Size: '+repr(1000*round(param[2]),3)+' micron'
	plt.text(.02,.97,fig_text,fontsize=10,bbox=bbox_props,va='top',ha='left',transform=ax.transAxes)
	if imsave!='none':
		ut.create_dir(imsave)
		plt.savefig(imsave+'beam_profile_'+str.zfill(z,5)+'.png')
	plt.show()
	
	return param[2],np.sqrt(covar[2,2]),pos,pows_norm
