# Script for handling image analysis

import numpy as np
import matplotlib.pyplot as plt
import utility as ut
import gen_reader as gr
import lv_analysis as lv
import scope_cal as scp
import funk
import math as m

# Define global variable for date and darkfile
date=20170720
darkfiles=gr.reader('/home/chris/anaconda/EXOscripts/darkframes.txt',header=False,dtype=int)
nullpulses=gr.reader('/home/chris/anaconda/EXOscripts/nullpulses.txt',header=False,dtype=int)
darkfile=int(darkfiles[np.where(darkfiles==date)[0],1])
if date in nullpulses[:,0]:
	null_pulse=int(nullpulses[np.where(nullpulses==date)[0],1])

#=========================================================================================

# Function to remove cosmic rays from images
def im_remove_cosmics_old(filenum,tag,output):

	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	import utility as ut
	
	# Import file to be read and it's information
	expnum=os.listdir('/home/chris/anaconda/data/'+str(date)+'/raw/ascii')[1][:7]
	filename='/home/chris/anaconda/data/'+str(date)+'/raw/ascii/'+expnum+str(filenum)+'.txt'
	datafile= open(filename,'r')
	
	# Get run number from data file name
	run=''
	j=1                                                        # j is the incrementor for characters in the run name/number
	for i in range(1,len(filename)+1):
		if filename[-i]=='/':                                  # Steps backwards in path until it finds a '/'
			while filename[-i+j]!='.':
				run += filename[-i+j]                          # Adds on characters after the '/' until it finds a '.'
				j+=1
			break
		else:
			continue

	# Initialize list for unorganized data. Going to make one long list of floats 
	# corresponding to the pixel values (i.e. column 4 of the datafile)
	column1=[]
	column2=[]
	column4=[]
	for line in datafile:
		line=line.strip()
		values=line.split()
		column1.append(float(values[0]))
		column4.append(float(values[3]))
		y_max=int(values[1])
		N= int(values[2])
	# Close the datafile, what were you, raised in MATLAB?
	datafile.close()
	
	# Get the number of x and y points and the arrays for each (usually 1340 x 400, but could change)
	x_max= int(len(column1)/(N*y_max))
	x=np.array(range(x_max))
	y=range(1,y_max+1)
	
	# Organize data into a list of lists of lists (frame,y,x).
	pix=[]
	# Big loop makes a list of the frame sub lists
	for i in range(N):
		frame_pix=[]
		# Small loop makes a list of the row sub lists
		for j in range(y_max*i,y_max*(i+1)):  
			# This grabs a chunk of length x_max from the unorganized list (one row of pixels)
			frame_pix.append(column4[x_max*(j):x_max*(j+1)]) 
		pix.append(frame_pix)
	
	# Make copy array for removal of cosmics
	pix_raw=np.array(pix)
	pix_cf=np.array(pix)
	
	# Remove Cosmic Rays
	for k in range(N):
		reps= 0
		rep_waves=[]
		for j in range(y_max):
			if j<=2:
				continue
			if j>= (y_max-3):
				continue
			else:
				
				for i in range(x_max): 
					if i <= 2:
						continue
					if i >= (x_max-3):
						continue
					else:
						bin= np.ndarray.flatten(pix_raw[k,j-2:j+3,i-2:i+3])
						bin=np.delete(bin,7)
					bin_ave= np.mean(bin)
					if pix_raw[k,j,i]>=(3*bin_ave):
						rep_vals= np.concatenate([bin,[bin_ave]])
						pix_cf[k,j,i-1:i+2]= np.amin(rep_vals)
						reps= reps+1
						rep_waves.append([round(x[i],2),round(y[j+1],2)])
		if reps > 0 and output==True:
			print('Cosmics removed from frame %d'%(k+1))
				
	return x,y,pix_cf,x_max,y_max,N,run

#-----------------------------------------------------------------------------------------

# Create dark frame average:
if darkfile!=0:
	dark_pix_g=np.mean(im_remove_cosmics_old(darkfile,'bg',False)[2],axis=0)
	#print(np.mean(dark_pix_g))
	print('dark run: '+str(darkfile))
if darkfile==0:
	dark_pix_g=0.

#------------------------------------------------------------------------------------------
# Function to import data from ascii files
def im_import(filenum,tag,output):

	# Import modules for array and plotting
	import numpy as np
	import os
	import pandas as pd
	
	# Import file to be read and it's information
	expnum=os.listdir('/home/chris/anaconda/data/'+str(date)+'/raw/ascii')[1][:7]
	filename='/home/chris/anaconda/data/'+str(date)+'/raw/ascii/'+expnum+str(filenum)+'.txt'
	
	# Get run number from data file name
	run=''
	j=1                                                        # j is the incrementor for characters in the run name/number
	for i in range(1,len(filename)+1):
		if filename[-i]=='/':                                  # Steps backwards in path until it finds a '/'
			while filename[-i+j]!='.':
				run += filename[-i+j]                          # Adds on characters after the '/' until it finds a '.'
				j+=1
			break
		else:
			continue

	# Import data using pandas
	dat=pd.read_csv(filename,header=None,sep=' ')
	dat.columns=['x','y','f','px']
	# Reshape into [f,y,x] ndarray
	y_max,N=np.amax(dat.y),np.amax(dat.f)
	x_max=int(len(dat.px)/(y_max*N))
	x,y=np.arange(x_max),np.arange(y_max)
	pix=np.reshape(dat.px,(N,y_max,x_max))
				
	return x,y,pix,x_max,y_max,N,run

#-----------------------------------------------------------------------------------------

# Function to background subtract images
def im_sub(sigfile,bgfile,**kwargs):
	
	# Define keyword arguments
	scale=kwargs.get('bg_scale',1.)
	shift=kwargs.get('bg_shift',[0,0])
	frame=kwargs.get('bg_frame','ave')
	sigframe=kwargs.get('frames','all')
	
	import numpy as np
	import utility as ut
	
	# Import data and remove cosmics for signal and background with image_remove_cosmics
	sig_x,sig_y,sig_pix,sig_x_max,sig_y_max,sig_N,sig_run=im_sub_dark(sigfile)
	bg_x,bg_y,bg_pix,bg_x_max,bg_y_max,bg_N,bg_run=im_sub_dark(bgfile)
	
	if sigframe!= 'all':
		sig_pix=sig_pix[sigframe,:,:]
		sig_N=len(sigframe)
	if frame=='ave':
		# Average the background frames together
		bg_pix=np.mean(bg_pix,axis=0)
	if frame=='each':
		bg_pix=bg_pix
	else:
		# Choose the frame to sub
		bg_pix=bg_pix[frame-1,:,:]
		
	# Shift if necessary
	bg_final=ut.shifter(bg_pix,shift[0],shift[1])
			
	# Do the subtraction
	pix_sub=sig_pix-(scale*bg_final)
	
	return sig_x,sig_y,pix_sub,sig_x_max,sig_y_max,sig_N,sig_run
	
#-----------------------------------------------------------------------------------------

# Function to background subtract images
def im_sub_dark(sigfile,speed,**kwargs):
	
	dark=kwargs.get('dark',False)
	
	import numpy as np
	import utility as ut
	
	# Import data and remove cosmics for signal and background with image_remove_cosmics
	sig_x,sig_y,sig_pix,sig_x_max,sig_y_max,sig_N,sig_run=im_import(sigfile,'signal',False)
	if speed=='slow':
		if dark==0.:
			dark_pix=np.zeros(np.shape(sig_pix))
			dark_pix[:,:]=0.
		if dark==False:
			dark_pix=dark_pix_g
		if dark not in [False,0.]:
			dark_x,dark_y,dark_pix,dark_x_max,dark_y_max,dark_N,dark_run=im_remove_cosmics_old(dark,'background',False)
			dark_pix=np.mean(dark_pix,axis=0)
			print('manual dark frame: '+str(dark))
			
	if speed=='fast':
		sig_pix=sig_pix[:,:,2:]
		dark_pix=np.zeros(np.shape(sig_pix))
	
	# Do the subtraction
	pix_sub=[]
	for i in range(sig_N):
		pix_sub.append(sig_pix[i,:,:]-dark_pix)
	pix_sub=np.array(pix_sub)
	if dark==0:
		pix_sub=sig_pix
	
	return sig_x,sig_y,pix_sub,sig_x_max,sig_y_max,sig_N,sig_run

#-----------------------------------------------------------------------------------------

# Function to plot 0th order images
def implot_0(filename,**kwargs):

	# Define keyword arguments
	colmap=kwargs.get('cmap','jet')
	frames=kwargs.get('frames','all')
	bgfile=kwargs.get('bgfile',False)
	store=kwargs.get('store',False)
	savepath=kwargs.get('savepath',False)
	xylim=kwargs.get('axis','default')
	zscale=kwargs.get('z_scale','default')
	tag=kwargs.get('tag','')
	bgscale=kwargs.get('bg_scale',1.)
	shift=kwargs.get('bg_shift',[0,0])
	bgframe=kwargs.get('bg_frame','ave')
	overlay=kwargs.get('overlay',False)
	rect=kwargs.get('rect','none')
	lvscale=kwargs.get('lv',False)
	speed=kwargs.get('speed','slow')
	scale=kwargs.get('scale',1.)
	dark=kwargs.get('dark',False)
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	import utility as ut
	import lv_analysis as lv
	
	if bgfile==False:
		# Import data and remove cosmics with image_remove_cosmics
		x,y,pix,x_max,y_max,N,run=im_sub_dark(filename,speed,dark=dark)
	else:
		# Import data and remove cosmics with image_sub
		if frames!='all':
			frames=np.array(frames)-1
		x,y,pix,x_max,y_max,N,run=im_sub(filename,bgfile,bg_scale=bgscale,bg_shift=shift,bg_frame=bgframe,frames=frames)
	
	if frames=='all':
		frames=range(1,N+1)
	
	# Scale, if neccesary
	pix=pix*scale
	
	if lvscale!=False:
		pows=lv.lv_energy(lvscale)*1000
		pix=pix/pows
					
	# Plot the frames desired
	if frames != 'none':
		for frame in frames:
			if overlay==False:
				plt.figure(run)
			if overlay==True:
				plt.figure('im_plot')	
			plt.clf()
			im=plt.imshow(pix[(frame-1),:,:],cmap=colmap,interpolation='nearest')
			plt.colorbar(im,orientation='vertical')
			plt.title(run+': frame '+repr(frame))
			plt.xlabel('x Pixel')
			plt.ylabel('y Pixel')
			# want to start at beginning of first lim, and end of last lim
			plt.axis([x[0]-.5,x[-1]-.5,y[0]-.5,y[-1]-.5])
			if xylim!='default':
				plt.axis(xylim)
			if zscale!='default':
				plt.clim(zscale)
			if rect!='none':
				ut.rect(rect)	
			plt.show()
			if savepath!=False:
				ut.create_dir(savepath)
				if frame<=9:
					framestr='00'+str(frame)
				if 10<=frame<=99:
					framestr='0'+str(frame)
				if 100<=frame:
					framestr=str(frame)
				plt.savefig(savepath+'/'+str(date)+'_run'+str(filename)+'_'+framestr+tag+'.png')
	 
	# Return data if desired
	if store==True:
		return pix
	
#-----------------------------------------------------------------------------------------

# Function to plot 0th order images
def im_read(filename,**kwargs):

	# Define keyword arguments
	bgfile=kwargs.get('bgfile',False)
	bgscale=kwargs.get('bg_scale',1.)
	shift=kwargs.get('bg_shift',[0,0])
	bgframe=kwargs.get('bg_frame','ave')
	lvscale=kwargs.get('lv',False)
	speed=kwargs.get('speed','slow')
	scale=kwargs.get('scale',1.)
	dark=kwargs.get('dark',False)
	T_R=kwargs.get('T_R',1.)
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	import utility as ut
	import lv_analysis as lv
	
	if bgfile==False:
		# Import data and remove cosmics with image_remove_cosmics
		x,y,pix,x_max,y_max,N,run=im_sub_dark(filename,speed,dark=dark)
	else:
		# Import data and remove cosmics with image_sub
		if frames!='all':
			frames=np.array(frames)-1
		x,y,pix,x_max,y_max,N,run=im_sub(filename,bgfile,bg_scale=bgscale,bg_shift=shift,bg_frame=bgframe,frames='all')
	
	# Scale, if neccesary
	pix=pix*scale
	if lvscale!=False:
		pows=lv.lv_pow(lvscale,T_R=T_R)*1000
		for i in range(len(pows)):
			pix[i,:,:]=pix[i,:,:]/pows[i]
					
	return pix

#-----------------------------------------------------------------------------------------

# Function to integrate area in 0-order image (can do fractional pixels)
def im_int(filename,bounds,**kwargs):

	# Define keyword arguments
	colmap=kwargs.get('cmap','jet')
	frames=kwargs.get('frames','all')
	bgfile=kwargs.get('bgfile',False)
	plt_rect=kwargs.get('plt_rect',True)
	savepath=kwargs.get('savepath',False)
	store=kwargs.get('store',True)
	overlay=kwargs.get('overlay',False)
	lab=kwargs.get('label','default')
	scale=kwargs.get('bg_scale',1.)
	shift=kwargs.get('bg_shift',[0,0])
	bgframe=kwargs.get('bg_frame','ave')
	xylim=kwargs.get('axis','default')
	zscale=kwargs.get('z_scale','default')
	submin=kwargs.get('sub_min',False)
	dt=kwargs.get('frametime',False)
	lvscale=kwargs.get('lv',False)
	timescale=kwargs.get('timescale',False)
	T_R=kwargs.get('T_R',1.)
	bndtype=kwargs.get('boundtype','edge')
	norm=kwargs.get('normalize',False)
	dark=kwargs.get('dark',False)
	powscale=kwargs.get('pow','frame')
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	import os
	import lv_analysis as lv
	
	if frames!='all':
			frames=np.array(frames)
				
	if bgfile==False:
		# Import data and remove cosmics with image_remove_cosmics
		x,y,pix,x_max,y_max,N,run=im_sub_dark(filename,'slow',dark=dark)
	else:
		# Import data and remove cosmics with image_sub
		x,y,pix,x_max,y_max,N,run=im_sub(filename,bgfile,bg_scale=scale,bg_shift=shift,bg_frame=bgframe,frames=frames)
	
	if frames=='all':
		frames=np.arange(N)+1	
	scaled=''
	
	if lvscale!=False:
		pows=lv.lv_energy(lvscale,T_R=T_R)
		if powscale=='ave':
			pix=pix/(1000*np.mean(pows))
		if powscale=='frame':
			for i in range(len(pows)):
				pix[i,:,:]=pix[i,:,:]/(1000*pows[i])
		scaled=' per mW'
		if timescale==True:
			times=lv.get_avtime(lvscale)
			pix=pix/(times*1E-3)
			scaled=' per mW*s'
		
	# Integrate over the bounds
	if bndtype=='edge':
		[l,r]=np.array(bounds[:2])
		[t,b]=np.array(bounds[-2:])
		intsize=(bounds[1]-bounds[0])*(bounds[3]-bounds[2])
	if bndtype=='center':
		[l,r]=np.array([bounds[0]-(bounds[2]/2.),bounds[0]+(bounds[2]/2.)])
		[t,b]=np.array([bounds[1]-(bounds[3]/2.),bounds[1]+(bounds[3]/2.)])
		intsize=bounds[2]*bounds[3]
	lf=np.ceil(l)-l
	rf=r-np.floor(r)
	tf=np.ceil(t)-t
	bf=b-np.floor(b)
	integral=[]
	for j in frames:
		i=j-1
		fullpix=np.sum(pix[i,m.ceil(t):m.floor(b),m.ceil(l):m.floor(r)])
		sidepixl=lf*np.sum(pix[i,m.ceil(t):m.floor(b),m.floor(l)])
		sidepixr=rf*np.sum(pix[i,m.ceil(t):m.floor(b),m.floor(r)])
		sidepixt=tf*np.sum(pix[i,m.floor(t),m.ceil(l):m.floor(r)])
		sidepixb=bf*np.sum(pix[i,m.floor(b),m.ceil(l):m.floor(r)])
		cornerbl=bf*lf*pix[i,m.floor(b),m.floor(l)]
		cornertl=tf*lf*pix[i,m.floor(t),m.floor(l)]
		cornertr=tf*rf*pix[i,m.floor(t),m.floor(r)]
		cornerbr=tf*lf*pix[i,m.floor(b),m.floor(r)]
		integral.append(fullpix+sidepixl+sidepixr+sidepixt+sidepixb+cornerbl+cornertl+cornerbr+cornertr)
	integral=np.array(integral,dtype=float)
	
	if norm==True:
		integral=integral/float(intsize)
	
	# Plot the frames desired
	if plt_rect==True:
		for frame in frames:
			plt.figure('int_area')
			plt.clf()
			im=plt.imshow(pix[(frame-1),:,:],cmap=colmap,interpolation='nearest')
			plt.colorbar(im,orientation='vertical')
			plt.title(run+': frame '+repr(frame))
			plt.xlabel('x Pixel')
			plt.ylabel('y Pixel')
			# want to start at beginning of first bound, and end of last bound
			plt.axis([x[0]-.5,x[-1]+.5,y[0]-.5,y[-1]+.5])
			ut.rect([l,r,t,b])
			# put pixel number at far pixel edge, instead of middle
			plt.xticks(np.arange(np.shape(pix)[2]),np.arange(np.shape(pix)[2])+1)
			plt.yticks(np.arange(np.shape(pix)[1]),np.arange(np.shape(pix)[1])+1)
			if xylim!='default':
				plt.axis(xylim)
			if zscale!='default':
				plt.clim(zscale)
			plt.show()
			if frame<=9:
				framestr='00'+str(frame)
			if 10<=frame<=99:
				framestr='0'+str(frame)
			if 100<=frame:
				framestr=str(frame)
			if savepath!=False:
				ut.create_dir(savepath)
				plt.savefig(savepath+'/'+str(date)+'_run'+str(filename)+'_'+framestr+'_int_area.png')
				
	
	# Plot integral
	if submin==True:
		integral=integral-np.amin(integral)
	if lab=='default':
		lab=run
	if overlay==False:
		plt.figure(run+'_int')
	if overlay==True:
		plt.figure('int')
	plt.clf()
	ax=plt.gca()
	if dt==False:
		plt.plot(frames,integral,'-o',label=lab)
		plt.axis([0,N+1,0,1.1*np.amax(integral)])
		plt.xlabel('Frame')
		xdata=frames
	if dt!=False:
		time=np.arange(dt,(np.shape(integral)[0]+1)*dt,dt)
		plt.plot(time,integral,'-o',label=lab)
		plt.axis([0,time[-1]+dt,0,1.1*np.amax(integral)])
		plt.xlabel('Time (s)')
		xdata=time	
	plt.ylabel('Integrated Counts'+scaled)
	plt.title(run+' Integral')
	bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
	fig_text='X-bounds: '+repr([l,r])+'\nY-bounds: '+repr([t,b])+'\nTotal Pixels: '+repr(np.size(pix[i,t:b,l:r]))
	plt.text(.05,.05,fig_text,fontsize=10,bbox=bbox_props,va='bottom',ha='left',transform=ax.transAxes)
	plt.show()
	if savepath!=False:
		plt.savefig(savepath+'/'+str(date)+'_run'+str(filename)+'_int_plot.png')
	
	if store==True:
		return integral	
#-----------------------------------------------------------------------------------------

def improf_0(filename,ax,bounds,**kwargs):
	
	lvscale=kwargs.get('lvscale',True)
	T_R=kwargs.get('T_R',1.)
	plot=kwargs.get('plot',False)
	
	import numpy as np
	import matplotlib.pyplot as plt
	
	# Read in the file
	x,y,pix,x_max,y_max,N,run=im_sub_dark(filename,speed='slow')
	
	# Scale by energy
	if lvscale==True:
		for i in range(np.shape(pix)[0]):
			pix[i,:,:]=pix[i,:,:]/(1e3*lv.lv_energy(filename,T_R=T_R)[i])
	
	# Integrate in one axis to form profile
	if ax=='x':
		prof=np.sum(pix[:,bounds[0]-1:bounds[1],:],axis=1)
	if ax=='y':
		prof=np.sum(pix[:,:,bounds[0]-1:bounds[1]],axis=2)
	
	# plot the profile
	if plot==True:
		plt.figure('profile')
		plt.clf()
		for i in range(np.shape(prof)[0]):
			plt.plot(prof[i,:],label='frame '+str(i+1))
		plt.legend(fontsize=11)
		plt.show()
		
	return prof
	
#-----------------------------------------------------------------------------------------

# Function to integrate each frame in the scan and make a 2d plot
def scan_plot(filename,lvnum,center,width,stepsize,numsteps,**kwargs):
	
	# Define keyword arguments
	colmap=kwargs.get('cmap','jet')
	frames=kwargs.get('frames','all')
	bgfile=kwargs.get('bgfile',False)
	plt_rect=kwargs.get('plt_rect',True)
	savepath=kwargs.get('savepath',False)
	store=kwargs.get('store',True)
	overlay=kwargs.get('overlay',False)
	lab=kwargs.get('label','default')
	scale=kwargs.get('bg_scale',1.)
	shift=kwargs.get('bg_shift',[0,0])
	bgframe=kwargs.get('bg_frame','ave')
	xylim=kwargs.get('axis','default')
	zscale=kwargs.get('z_scale','default')
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	import os
	import lv_analysis as lv
		
	if bgfile==False:
		# Import data and remove cosmics with image_remove_cosmics
		x,y,pix,x_max,y_max,N,run=im_sub_dark(filename)
	else:
		# Import data and remove cosmics with image_sub
		if frames!='all':
			frames=np.array(frames)-1
		x,y,pix,x_max,y_max,N,run=im_sub(filename,bgfile,bg_scale=scale,bg_shift=shift,bg_frame=bgframe,frames=frames)
				
	# Make a box centered at some pixel to integrate
	bounds=np.zeros((N,4))
	bounds[0,:]=[center[0]-.5*width[0],center[0]+.5*width[0],center[1]-.5*width[1],center[1]+.5*width[1]]
	xinc,xdir,ymoved=0,1,0
	stepsizepix=np.array(stepsize)/5. 
	for i in range(N):
		if i==0: 
			continue
		if xinc==numsteps[0]-1:
			bounds[i,:]=[bounds[i-1,0],bounds[i-1,1],bounds[i-1,2]+stepsizepix[1],bounds[i-1,3]+stepsizepix[1]]
			ymoved=1
			xdir=xdir*-1
		if xinc<numsteps[0]-1:
			bounds[i,:]=[bounds[i-1,0]+xdir*stepsizepix[0],bounds[i-1,1]+xdir*stepsizepix[0],bounds[i-1,2],bounds[i-1,3]]
			xinc+=1
		if ymoved==1:
			ymoved,xinc=0,0
	
	# Integrate over the bounds
	xbounds=bounds[:,:2]-1
	ybounds=bounds[:,-2:]-1
	xbounds_i=np.ceil(xbounds)
	ybounds_i=np.ceil(ybounds)
	xbounds_f=xbounds-np.floor(xbounds)
	ybounds_f=ybounds-np.floor(ybounds)
	integral=np.zeros ((N,),dtype=float)
	for i in range(N):
		fullpix=np.sum(pix[i,ybounds_i[i,0]+1:ybounds_i[i,1],xbounds_i[i,0]+1:xbounds_i[i,1]])
		sidepixl=(1-xbounds_f[i,0])*np.sum(pix[i,ybounds_i[i,0]+1:ybounds_i[i,1],xbounds_i[i,0]])
		sidepixr=xbounds_f[i,1]*np.sum(pix[i,ybounds_i[i,0]+1:ybounds_i[i,1],xbounds_i[i,1]])
		sidepixb=(1-ybounds_f[i,0])*np.sum(pix[i,ybounds_i[i,0],xbounds_i[i,0]+1:xbounds_i[i,1]])
		sidepixt=ybounds_f[i,1]*np.sum(pix[i,ybounds_i[i,1],xbounds_i[i,0]+1:xbounds_i[i,1]])
		cornerbl=(1-xbounds_f[i,0])*(1-ybounds_f[i,0])*pix[i,ybounds_i[i,0],xbounds_i[i,0]]
		cornertl=(1-xbounds_f[i,0])*ybounds_f[i,1]*pix[i,ybounds_i[i,1],xbounds_i[i,0]]
		cornertr=xbounds_f[i,1]*ybounds_f[i,1]*pix[i,ybounds_i[i,1],xbounds_i[i,1]]
		cornerbr=xbounds_f[i,1]*(1-ybounds_f[i,0])*pix[i,ybounds_i[i,0],xbounds_i[i,1]]
		integral[i]=fullpix+sidepixl+sidepixr+sidepixt+sidepixb+cornerbl+cornertl+cornerbr+cornertr
	# Power scale
	integral=integral/lv.get_pow(lvnum)
	
	# Plot the frames desired
	for frame in range(1,N+1):
		plt.figure('int_area_scan')
		plt.clf()
		im=plt.imshow(pix[(frame-1),:,:],cmap=colmap,interpolation='nearest')
		plt.colorbar(im,orientation='vertical')
		plt.title(run+': frame '+repr(frame))
		plt.xlabel('x Pixel')
		plt.ylabel('y Pixel')
		# want to start at beginning of first bound, and end of last bound
		plt.axis([x[0]-.5,x[-1]+.5,y[0]-.5,y[-1]+.5])
		ut.rect(bounds[frame-1,:])
		# put pixel number at far pixel edge, instead of middle
		plt.xticks(np.arange(np.shape(pix)[2])+.5,np.arange(np.shape(pix)[2]))
		plt.yticks(np.arange(np.shape(pix)[1])+.5,np.arange(np.shape(pix)[1]))
		if xylim!='default':
			plt.axis(xylim)
		plt.clim(0,1.05*np.amax(pix))
		if zscale!='default':
			plt.clim(zscale)
		plt.show()
		if frame<=9:
			framestr='00'+str(frame)
		if 10<=frame<=99:
			framestr='0'+str(frame)
		if 100<=frame:
			framestr=str(frame)
		if savepath!=False:
			ut.create_dir(savepath)
			plt.savefig(savepath+'/'+run+'_'+framestr+'_int_area.png')
	
	# Create integral map in 2d
	integral2d=integral.reshape(np.flipud(numsteps))
	# Flip every other row since scan is rasterized
	for i in range(numsteps[1]):
		if i%2==1:
			integral2d[i,:]=np.flipud(integral2d[i,:])
	
	# Plot integral map
	if lab=='default':
		lab=run
	if overlay==False:
		plt.figure(run+'_int2d')
	if overlay==True:
		plt.figure('int2d')
	plt.clf()
	im=plt.imshow(integral2d,cmap=colmap,interpolation='nearest',extent=[0,numsteps[0]*stepsize[0],0,numsteps[1]*stepsize[1]])
	plt.colorbar(im,orientation='vertical')
	plt.title(run+'Integral Scan')
	plt.xlabel('x location (micron)')
	plt.ylabel('y location (micron)')
	if zscale!='default':
		plt.clim(zscale)
	plt.show()
	if savepath!=False:
		plt.savefig(savepath+'/'+run+'_intscan.png')
	
	if store==True:
		return integral2d,bounds	

#------------------------------------------------------------------------------------------

# Function to integrate each frame in the scan and make a 2d plot
def scanplot_static(filename,bounds,numsteps,**kwargs):
	
	# Define keyword arguments
	colmap=kwargs.get('cmap','jet')
	plt_rect=kwargs.get('plt_rect',True)
	savepath=kwargs.get('savepath',False)
	store=kwargs.get('store',True)
	xylim=kwargs.get('axis','default')
	zscale=kwargs.get('z_scale','default')
	manrun=kwargs.get('run','expXXX_YYY')
	lvnum=kwargs.get('lv','none')
	T_R=kwargs.get('T_R',1.)
	drop=kwargs.get('drop',[])
	sap=kwargs.get('sapphire',False)
	dark=kwargs.get('dark',False)
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	import os
	import gen_reader as gr
	from lmfit import minimize,Parameters,Parameter,Model
	import lv_analysis as lv
	
	if sap==True:
		savepath=savepath[:-1]
		savepath+='_sapphire/'
	ut.create_dir(savepath)
		
	# Import data and remove cosmics with image_remove_cosmics
	if isinstance(filename,int)==True:
		x,y,pix,x_max,y_max,N,run=im_sub_dark(filename,'slow',dark=dark)
	# If already imported, just use array
	if isinstance(filename,np.ndarray)==True:
		pix=filename
		N=np.shape(pix)[0]
		run=manrun
		
	
	# Integrate the data in the bounds and plot the data and box
	integral=[]
	for i in range(N):
		integral.append(ut.integrate(pix[i,:,:],bounds,boundtype='edge'))
		if i+1 in drop:
			integral.pop(-1)
		if plt_rect==True and savepath!=False:
			plt.figure('int_area')
			plt.clf()
			im=plt.imshow(pix[i,:,:],cmap=colmap,interpolation='nearest')
			plt.colorbar(im,orientation='vertical')
			plt.title(run+': frame '+repr(i+1))
			plt.xlabel('x Pixel')
			plt.ylabel('y Pixel')
			if xylim!='default':
				plt.axis(xylim)
			ut.rect(bounds)
			if xylim!='default':
				plt.axis(xylim)
			plt.show()
			if savepath!=False:
				if (i+1)<=9:
					framestr='00'+str(i+1)
				if 10<=(i+1)<=99:
					framestr='0'+str(i+1)
				if 100<=(i+1)																																																																																																																																																																																			:
					framestr=str(i+1)
				plt.savefig(savepath+'/'+str(date)+'_run'+str(filename)+'_'+framestr+'_int_area.png')
	
	for i in drop:
		integral=np.append(integral,50)
		
	if lvnum!='none':
		pows=lv.lv_pow(lvnum,T_R=T_R)
		for i in drop:
			pows=np.insert(pows,i,pows[i])
		integral=integral/(np.absolute(pows)*1E3)
	if len(drop)!=0:
		intmin,intminf=np.amin(integral[:-len(drop)]),np.argmin(integral[:-len(drop)])+1
	if len(drop)==0:
		intmin,intminf=np.amin(integral),np.argmin(integral)+1
	print(drop)
	print(intmin,intminf)

	# Create integral map in 2d
	integral2d=integral.reshape(np.flipud(numsteps))
	# Flip every other row since scan is rasterized
	for i in range(numsteps[1]):
		if i%2==1:
			integral2d[i,:]=np.flipud(integral2d[i,:])
	
	# Plot integral map
	plt.figure('int2d')
	plt.clf()
	im=plt.imshow(integral2d,cmap=colmap,interpolation='nearest',origin='lower')
	cb=plt.colorbar(im,orientation='vertical')
	cb.formatter.set_powerlimits((0, 0))
	cb.update_ticks()
	plt.title(run+' Integral Scan')
	plt.xlabel('x step')
	plt.ylabel('y step')
	ut.textbox('Minimum: '+str(int(intmin))+', frame: '+str(intminf),[.05,.95])
	if zscale!='default':
		plt.clim(zscale)
	plt.show()
	if savepath!=False:
		plt.savefig(savepath+'/'+str(date)+'_run'+str(filename)+'_intscan_static.png')
		
	if store==True:
		return integral2d

#------------------------------------------------------------------------------------------

# Function to integrate each frame in the scan and make a 2d plot
def scanplot_fixed(filename,startpos,stepsize,numsteps,**kwargs):
	
	# Define keyword arguments
	colmap=kwargs.get('cmap','jet')
	plt_rect=kwargs.get('plt_rect',True)
	savepath=kwargs.get('savepath',False)
	store=kwargs.get('store',True)
	xylim=kwargs.get('axis','default')
	zscale=kwargs.get('z_scale','default')
	intsize=kwargs.get('intsize',[5,5])
	manrun=kwargs.get('run','expXXX_YYY')
	lvnum=kwargs.get('lv','none')
	T_R=kwargs.get('T_R',1.)
	drop=kwargs.get('drop',[])
	sap=kwargs.get('sapphire',False)
	dark=kwargs.get('dark',False)
	pixscale=kwargs.get('pixscale','default')
	cs=kwargs.get('checkstart',False)
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	import os
	import gen_reader as gr
	from lmfit import minimize,Parameters,Parameter,Model
	import lv_analysis as lv
	
	if sap==True:
		savepath=savepath[:-1]
		savepath+='_sapphire/'
	ut.create_dir(savepath)
	ut.create_dir(savepath+'framedata')
		
	# Import data and remove cosmics with image_remove_cosmics
	if isinstance(filename,int)==True:
		x,y,pix,x_max,y_max,N,run=im_sub_dark(filename,'slow',dark=dark)
	# If already imported, just use array
	if isinstance(filename,np.ndarray)==True:
		pix=filename
		N=np.shape(pix)[0]
		run=manrun
		
	# Make a box centered at some pixel to integrate
	bounds=np.zeros((N,4))
	bounds[0,:]=[startpos[0]-.5*intsize[0],startpos[0]+.5*intsize[0],startpos[1]-.5*intsize[1],startpos[1]+.5*intsize[1]]
	if sap==True:
		bounds[0,2:]=bounds[0,2:]-12
	xinc,xdir,ymoved=0,1,0
	stepsizepix=np.array(stepsize)/5. 
	for i in range(N):
		if i==0: 
			continue
		if i+1 in drop:
			bounds[i,:]=bounds[i-1,:]
			continue
		if xinc==numsteps[0]-1:
			bounds[i,:]=[bounds[i-1,0],bounds[i-1,1],bounds[i-1,2]+stepsizepix[1],bounds[i-1,3]+stepsizepix[1]]
			ymoved=1
			xdir=xdir*-1
		if xinc<numsteps[0]-1:
			bounds[i,:]=[bounds[i-1,0]+xdir*stepsizepix[0],bounds[i-1,1]+xdir*stepsizepix[0],bounds[i-1,2],bounds[i-1,3]]
			xinc+=1
		if ymoved==1:
			ymoved,xinc=0,0
	
	# Integrate the data in the bounds and plot the data and box
	integral=[]
	for i in range(N):
		integral.append(ut.integrate(pix[i,:,:],bounds[i,:],boundtype='edge'))
		if i+1 in drop:
			integral.pop(-1)
		if plt_rect==True and savepath!=False:
			plt.figure('int_area')
			plt.clf()
			im=plt.imshow(pix[i,:,:],cmap=colmap,interpolation='nearest')
			plt.colorbar(im,orientation='vertical')
			plt.title(run+': frame '+repr(i+1))
			plt.xlabel('x Pixel')
			plt.ylabel('y Pixel')
			if xylim!='default':
				plt.axis(xylim)
			ut.rect(bounds[i,:])
			if pixscale!='default':
				plt.clim(pixscale)
			plt.show()
			if savepath!=False:
				if (i+1)<=9:
					framestr='00'+str(i+1)
				if 10<=(i+1)<=99:
					framestr='0'+str(i+1)
				if 100<=(i+1)																																																																																																																																																																																			:
					framestr=str(i+1)
				plt.savefig(savepath+'/'+str(date)+'_run'+str(filename)+'_'+framestr+'_int_area.png')
				ut.save(np.ravel(pix[i,xylim[2]:xylim[3],xylim[0]:xylim[1]]),['data'],savepath+'framedata/'+str(date)+'_run'+str(filename)+'_'+framestr+'.txt')
	
	if cs==True:
		plt.figure('start')
		plt.clf()
		plt.imshow(pix[0,:,:],cmap=colmap,interpolation='nearest')
		ut.rect(bounds[0,:])
		if xylim!='default':
				plt.axis(xylim)
		plt.show()
	
	integral=np.array(integral)
	
	for i in drop:
		integral=np.append(integral,50)
		
	if lvnum!='none':
		pows=lv.lv_energy(lvnum,T_R=T_R)
		for i in drop:
			pows=np.insert(pows,i,pows[i])
		integral=integral/(np.absolute(pows)*1E3)
	if len(drop)!=0:
		intmin,intminf=np.amin(integral[:-len(drop)]),np.argmin(integral[:-len(drop)])+1
	if len(drop)==0:
		intmin,intminf=np.amin(integral),np.argmin(integral)+1
	print(drop)
	print(intmin,intminf)

	# Create integral map in 2d
	integral2d=integral.reshape(np.flipud(numsteps))
	# Flip every other row since scan is rasterized
	for i in range(numsteps[1]):
		if i%2==1:
			integral2d[i,:]=np.flipud(integral2d[i,:])
	
	# Plot integral map
	plt.figure('int2d')
	plt.clf()
	im=plt.imshow(integral2d,cmap=colmap,interpolation='nearest',origin='lower')
	cb=plt.colorbar(im,orientation='vertical')
	cb.formatter.set_powerlimits((0, 0))
	cb.update_ticks()
	plt.title(run+' Integral Scan '+str(intsize[0])+'x'+str(intsize[1]))
	plt.xlabel('x step')
	plt.ylabel('y step')
	#ut.textbox('Minimum: '+str(int(intmin))+', frame: '+str(intminf),[.05,.95])
	if zscale!='default':
		plt.clim(zscale)
	plt.show()
	if savepath!=False:
		plt.savefig(savepath+'/'+str(date)+'_run'+str(filename)+'_intscan_'+str(intsize[0])+'x'+str(intsize[1])+'.png')
		ut.save(np.array([startpos[0],startpos[1],stepsize[0],stepsize[1],intmin,intminf]),['xstart','ystart','xstep','ystep','minimum','min frame'],savepath+'/'+str(date)+'_run'+str(run)+'_scandata.txt')
	
	if store==True:
		return integral2d,bounds,integral

#------------------------------------------------------------------------------------------

# Function to integrate each frame in the scan and make a 2d plot
def scanplot_fit(filename,fitbnds,stepsize,numsteps,**kwargs):
	
	# Define keyword arguments
	colmap=kwargs.get('cmap','jet')
	plt_rect=kwargs.get('plt_rect',True)
	savepath=kwargs.get('savepath','/home/chris/anaconda/plots/'+str(date))
	store=kwargs.get('store',True)
	xylim=kwargs.get('axis','default')
	zscale=kwargs.get('z_scale','default')
	intsize=kwargs.get('intsize',3.)
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	import os
	from lmfit import minimize,Parameters,Parameter,Model
		
	# Import data and remove cosmics with image_remove_cosmics
	x,y,pix,x_max,y_max,N,run=im_sub_dark(filename,'slow')
	
	# Slice data for just fit bounds
	pixslice=pix[:,fitbnds[2]:fitbnds[3],fitbnds[0]:fitbnds[1]]
	x,y=np.arange(fitbnds[0],fitbnds[1]),np.arange(fitbnds[2],fitbnds[3])
	
	# Define 2d gaussian
	def gauss2d_flat(x,y,amp,xcenter,ycenter,xwidth,ywidth):  # This ends up as f(y,x)
		my,mx=np.meshgrid(x,y)
		xwidth,ywidth=xwidth/2.,ywidth/2.	
		return np.ravel(amp*np.exp(-2*(np.square(mx-ycenter)/(np.square(ywidth))))*np.exp(-2*(np.square(my-xcenter))/(np.square(xwidth))))
	
	# Find the peak and integrate region around it
	integral,param_all,cov=[],[],[]
	p0=[100,np.mean([fitbnds[0],fitbnds[1]]),np.mean([fitbnds[2],fitbnds[3]]),4,4]
	p_0s=['amp','xcenter','ycenter','xwidth','ywidth']
	func_model=Model(gauss2d_flat,independent_vars=['x','y'],param_names=p_0s)
	for i in range(N):
		
		result=func_model.fit(np.ravel(pixsclice[i,:,:]),x=x,y=y,amp=p0[0],xcenter=p0[1],ycenter=p0[2],xwidth=p0[3],ywidth=p0[4],verbose=False)
		center=[result.best_values['xcenter'],result.best_values['ycenter']]
		param_all.append(center)
		integral.append(ut.integrate(data,[center[0]-fitbnds[0],center[1]-fitbnds[2],intsize,intsize]))
		if plt_rect==True:
			plt.figure('int_area')
			plt.clf()
			im=plt.imshow(data,cmap=colmap,interpolation='nearest')
			plt.colorbar(im,orientation='vertical')
			plt.title(run+': frame '+repr(i+1))
			plt.xlabel('x Pixel')
			plt.ylabel('y Pixel')
			ut.rect([center[0]-intsize/2.,center[0]+intsize/2.,center[1]-intsize/2.,center[1]+intsize/2.])
			if xylim!='default':
				plt.axis(xylim)
			if zscale!='default':
				plt.clim(zscale)
			plt.show()
			if (i+1)<=9:
				framestr='00'+str(i+1)
			if 10<=(i+1)<=99:
				framestr='0'+str(i+1)
			if 100<=(i+1)																																																																																																																																																																																			:
				framestr=str(i+1)
			plt.savefig(savepath+'/'+str(date)+'_run'+str(filename)+'_'+framestr+'_int_area.png')
	param_all=np.array(param_all)
	cov=np.array(cov)
	integral=np.array(integral)

	# Create integral map in 2d
	integral2d=integral.reshape(np.flipud(numsteps))
	# Flip every other row since scan is rasterized
	for i in range(numsteps[1]):
		if i%2==1:
			integral2d[i,:]=np.flipud(integral2d[i,:])
	
	# Plot integral map
	plt.figure('int2d')
	plt.clf()
	im=plt.imshow(integral2d,cmap=colmap,interpolation='nearest')
	plt.colorbar(im,orientation='vertical')
	plt.title(run+'Integral Scan')
	plt.xlabel('x step ('+str(stepsize)+' micron)')
	plt.ylabel('y step ('+str(stepsize)+' micron)')
	if zscale!='default':
		plt.clim(zscale)
	plt.show()
	plt.savefig(savepath+'/'+str(date)+'_run'+str(filename)+'_intscan.png')
	# Plot center positions
	plt.figure('center_pos')
	plt.clf()
	plt.plot(param_all[:,0],param_all[:,1],'b-o')
	plt.title(run+'Center Fit Positions')
	plt.xlabel('x position (pixel)')
	plt.ylabel('y position (pixel)')
	plt.show()
	plt.savefig(savepath+'/'+str(date)+'_run'+str(filename)+'_centerpos.png')
	
	if store==True:
		return integral2d,param_all,integral	

#------------------------------------------------------------------------------------------
# Function to plot 1st order images
def image_plot_1(filename,**kwargs):

	# Define keyword arguments
	colmap=kwargs.get('cmap','jet')
	dc=kwargs.get('disp_cosmics',False)
	frames=kwargs.get('frames','all')
	df=kwargs.get('diff_frames',False)
	sub=kwargs.get('bgfile',False)
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as u
	
	if bgfile==False:
		# Import data and remove cosmics with image_remove_cosmics
		wavelength,y,pix,wave_max,y_max,N,run=im_remove_cosmics(filename,'image',dc)
	else:
		# Import data and remove cosmics with image_sub
		wavelength,y,pix,wave_max,y_max,N,run=im_sub(filename,bgfile,dc)
	
	# Calculate the aspect ratio based on the size and number of pixels
	extent=[wavelength[0],wavelength[-1],1,y_max]
	# Calculate the length/pixel in wavelength and y directions
	dx=(wavelength[-1]-wavelength[0])/float(wave_max)
	dy=(y[-1]-y[0])/float(y_max)
	aspect_ratio=dx/dy
	
	if frames=='all':
		frames=range(1,N+1)
	
	# Plot the frames desired
	if df==True:
		for frame in frames:
			fig_new=plt.figure(aspect_ratio,1)
			im=plt.imshow(pix[(frame-1),:,:],extent=extent,aspect=aspect_ratio,cmap=colmap,interpolation='nearest')
			plt.colorbar(im,orientation='vertical')
			plt.title(run+': frame '+repr(frame))
			plt.xlabel('Wavelength (nm)')
			plt.ylabel('y Pixel')
			plt.axis([wavelength[0],wavelength[-1],y[0],y[-1]])
			plt.show()
	if df==False:
		fig=plt.figure(3/aspect_ratio,3)
		for frame in frames:
			im=plt.imshow(pix[(frame-1),:,:],extent=extent,aspect=aspect_ratio,cmap=colmap,interpolation='nearest',label='Frame '+repr(frame))
		plt.colorbar(im,orientation='vertical')
		plt.title(run)
		plt.xlabel('Wavelength (nm)')
		plt.ylabel('y Pixel')
		plt.axis([wavelength[0],wavelength[-1],y[0],y[-1]])
		plt.legend()
		plt.show()
	
#-----------------------------------------------------------------------------------------

# Function to plot integrated 1st order images (turn them into spectra)
def image_int_1(filename,bounds,**kwargs):
	
	# Define keyword arguments (define frames later)
	dc=kwargs.get('disp_cosmics',False)
	frames=kwargs.get('frames','all')
	df=kwargs.get('diff_frames',False)
	sub=kwargs.get('bgfile',False)
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as u
	
	if bgfile==False:
		# Import data and remove cosmics with image_remove_cosmics
		wavelength,y,pix,wave_max,y_max,N,run=image_remove_cosmics(filename,'image',dc)
	else:
		# Import data, remove cosmics, and background subtract with image_sub
		wavelength,y,pix,wave_max,y_max,N,run=image_sub(filename,bgfile,disp_cosmics=dc)
	
	# Integrate over bounds
	pix_int=np.sum(pix[:,bounds[0]:bounds[1],:],1)
	
	if frames=='all':
		frames=range(1,N+1)
	
	# Plot the frames desired
	if df==True:
		for frame in frames:
			fig_new=plt.figure()
			plt.plot(wavelength,pix_int[(frame-1),:])
			plt.title(run+': frame '+repr(frame)+' ROI: '+repr(bounds))
			plt.xlabel('Wavelength (nm)')
			plt.ylabel('Counts in ROI')
			plt.axis([wavelength[0],wavelength[-1],0,1.25*np.amax(pix_int[(frame-1),:])])
			plt.show()
	if df==False:
		fig=plt.figure()
		max_counts=[]
		for frame in frames:
			plt.plot(wavelength,pix_int[(frame-1),:],label='Frame '+repr(frame))
			max_counts.append(np.amax(pix_int[(frame-1),:]))
		plt.title(run+' ROI: '+repr(bounds))
		plt.xlabel('Wavelength (nm)')
		plt.ylabel('Counts in ROI')
		plt.axis([wavelength[0],wavelength[-1],0,1.25*np.amax(np.array(max_counts))])
		plt.legend()
		plt.show()
		
#-----------------------------------------------------------------------------------------

# Function to get laser profile from images (zinfo=[start,step])
def focus_find(filelist,bounds,zinfo,**kwargs):
	
	savepath=kwargs.get('savepath',False)
	fig=kwargs.get('fig','')
	speed=kwargs.get('speed','slow')
	store=kwargs.get('store',False)
	
	# Import modules
	import numpy as np
	import matplotlib.pyplot as plt
	import math
	from scipy.optimize import curve_fit
	import os
	import utility as ut
	
	if savepath!=False:
		ut.create_dir(savepath)
		
	# Get data from files
	allpix=[]
	for filename in filelist:
		x,y,pix,x_max,y_max,N,run=im_sub_dark(filename,speed=speed)
		allpix.append(pix)
	allpix=np.array(allpix)[:,0,bounds[2]:bounds[3],bounds[0]:bounds[1]]
	x,y=x[bounds[0]:bounds[1]],y[bounds[2]:bounds[3]]

	# Do a 2D gauss fit to each file
	widths=np.zeros((2,np.shape(allpix)[0]))
	for i in range(np.shape(allpix)[0]):
		params=ut.gaussfit2d(allpix[i,:,:],x,y,[100,np.median(x),np.median(y),4,4])
		widths[:,i]=[params['xwidth'],params['ywidth']]
	
	# Plot width vs micrometer setting
	z=np.arange(zinfo[0],zinfo[0]+len(filelist)*zinfo[1],zinfo[1])
	plt.figure('beam_profs')
	plt.clf()
	plt.plot(z,5*widths[0,:],'bo',label='x profile')
	plt.plot(z,5*widths[1,:],'ro',label='y profile')
	plt.xlabel('micrometer position (mm)')
	plt.ylabel('width (micron)')
	plt.legend()
	plt.title(str(date)+' Laser Focus Finding')
	if savepath!=False:
		plt.savefig(savepath+str(date)+'_focus_finder.png')
	
	if store==True:
		return z,widths

#------------------------------------------------------------------------------------------

# Function to integrate area in 0-order image
def centroid(filename,bounds,**kwargs):

	# Define keyword arguments
	colmap=kwargs.get('cmap','jet')
	frames=kwargs.get('frames','all')
	bgfile=kwargs.get('bgfile',False)
	plt_rect=kwargs.get('plt_rect',True)
	savepath=kwargs.get('savepath',False)
	overlay=kwargs.get('overlay',False)
	scale=kwargs.get('bg_scale',1.)
	shift=kwargs.get('bg_shift',[0,0])
	bgframe=kwargs.get('bg_frame','ave')
	xylim=kwargs.get('axis','default')
	zscale=kwargs.get('z_scale','default')
	stats=kwargs.get('stats',False)
	annotate=kwargs.get('annotate',False) 
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	import os
	
	if savepath!=False:
		ut.create_dir(savepath)
		
	if bgfile==False:
		# Import data and remove cosmics with image_remove_cosmics
		x,y,pix,x_max,y_max,N,run=im_sub_dark(filename,speed='slow')
	
	# Import data and remove cosmics with image_sub
	else:
		if frames!='all':
			frames=np.array(frames)-1
		x,y,pix,x_max,y_max,N,run=im_sub(filename,bgfile,bg_scale=scale,bg_shift=shift,bg_frame=bgframe,frames=frames)
	
	if frames=='all':
		frames=range(1,N+1)
	frames=np.array(frames)
					
	# Find the Centroid
	bounds=np.array(bounds)-1
	cent=[]
	for i in frames-1:
		centx,centy=[],[]
		for j in range(bounds[2],bounds[3]+1):
			for k in range(bounds[0],bounds[1]+1):
				centx.append(pix[i,j,k]*(k+1))
				centy.append(pix[i,j,k]*(j+1))
		cent.append(((1./np.sum(pix[i,bounds[2]:bounds[3]+1,bounds[0]:bounds[1]+1]))*np.sum(centx),(1./np.sum(pix[i,bounds[2]:bounds[3]+1,bounds[0]:bounds[1]+1]))*np.sum(centy)))
	cent=np.array(cent)
	
	# Plot the frames desired
	if plt_rect==True:
		j=0
		for frame in frames:
			plt.figure('cent_area')
			plt.clf()
			im=plt.imshow(pix[(frame-1),:,:],cmap=colmap,interpolation='nearest')
			plt.plot(cent[j,0],cent[j,1],'ko',markersize=5)
			plt.colorbar(im,orientation='vertical')
			plt.title(run+': frame '+repr(frame))
			plt.xlabel('x Pixel')
			plt.ylabel('y Pixel')
			plt.axis([bounds[0]-2,bounds[1]+2,bounds[2]-2,bounds[3]+2])
			ut.rect(bounds)
			if xylim!='default':
				plt.axis(xylim)
			if zscale!='default':
				plt.clim(zscale)
			plt.show()
			if savepath!=False:
				ut.create_dir(savepath)
				framestr=str(frame).zfill(4)
				plt.savefig(savepath+str(date)+'_run'+str(filename)+'_'+framestr+'_cent_area.png')
			j+=1
	plt.figure('cent_all')
	plt.clf()
	
	# Plot each centroid position
	plt.figure('cent_all')
	plt.plot(cent[:,0],cent[:,1],'bo')
	plt.title('Centroid for All Frames in '+run)
	plt.xlabel('pixel')
	plt.ylabel('pixel')
	if annotate==True:
		for i in frames-1:
			plt.annotate(str(i+1),xy=(cent[i,0],cent[i,1]),xytext=(5,5),textcoords='offset points')
	if savepath!=False:
		plt.savefig(savepath+str(date)+'_run'+str(filename)+'centroid_all.png')
		
	if stats==True:
		
		# Calculate Average and Standard Deviation
		cent_ave=np.mean(cent,axis=0)
		cent_std=np.std(cent,axis=0)

		# Plot average position
		plt.figure('cent_all')
		plt.plot(cent_ave[0],cent_ave[1],'ro',markersize=8)
		ut.textbox('Average Centroid: ('+str(round(cent_ave[0],2))+','+str(round(cent_ave[1],2))+')\nCentroid StDev: ('+str(round(cent_std[0],2))+','+str(round(cent_std[1],2))+')',[.05,.95])
		if savepath!=False:
			plt.savefig(savepath+'centroid/'+str(date)+'_run'+str(filename)+'centroid_all.png')
		ut.rem('cent_all','text',0)
		plt.axis([cent_ave[0]-3*cent_std[0],cent_ave[0]+3*cent_std[0],cent_ave[1]-3*cent_std[1],cent_ave[1]+3*cent_std[1]])
		ut.textbox('Average Centroid: ('+str(round(cent_ave[0],2))+','+str(round(cent_ave[1],2))+')\nCentroid StDev: ('+str(round(cent_std[0],2))+','+str(round(cent_std[1],2))+')',[.05,.95])
		if savepath!=False:
			plt.savefig(savepath+str(date)+'_run'+str(filename)+'centroid_all_zoom.png')
		return cent,cent_ave,cent_std
	
	return cent	
#------------------------------------------------------------------------------------------

# Function to integrate area in 0-order image
def centroid_scan(filename,center,width,skip,numsteps,**kwargs):

	# Define keyword arguments
	colmap=kwargs.get('cmap','jet')
	frames=kwargs.get('frames','all')
	bgfile=kwargs.get('bgfile',False)
	plt_rect=kwargs.get('plt_rect',True)
	savepath=kwargs.get('savepath',False)
	overlay=kwargs.get('overlay',False)
	scale=kwargs.get('bg_scale',1.)
	shift=kwargs.get('bg_shift',[0,0])
	bgframe=kwargs.get('bg_frame','ave')
	xylim=kwargs.get('axis','default')
	zscale=kwargs.get('z_scale','default')
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	import os
	
	if savepath!=False:
		ut.create_dir(savepath+'centroid/')
		ut.create_dir(savepath+'centroid/im/')
		
	if bgfile==False:
		# Import data and remove cosmics with image_remove_cosmics
		x,y,pix,x_max,y_max,N,run=im_sub_dark(filename)
	
	# Import data and remove cosmics with image_sub
	else:
		if frames!='all':
			frames=np.array(frames)-1
		x,y,pix,x_max,y_max,N,run=im_sub(filename,bgfile,bg_scale=scale,bg_shift=shift,bg_frame=bgframe,frames=frames)
	
	if frames=='all':
		frames=range(1,N+1)
	frames=np.array(frames)
	
	# Make a box centered at some pixel to integrate
	bounds=np.zeros((N,4),dtype=int)
	bounds[0,:]=[center[0]-width[0],center[0]+width[0],center[1]-width[1],center[1]+width[1]]
	stepsize=[1,1]
	xinc,xdir,ymoved=0,1,0
	for i in range(N):
		if i==0: 
			continue
		if xinc==numsteps[0]-1:
			bounds[i,:]=[bounds[i-1,0],bounds[i-1,1],bounds[i-1,2]+stepsize[1],bounds[i-1,3]+stepsize[1]]
			ymoved=1
			xdir=xdir*-1
		if xinc<numsteps[0]-1:
			if xinc % skip[0] == 0:
				bounds[i,:]=[bounds[i-1,0]+xdir*stepsize[0],bounds[i-1,1]+xdir*stepsize[0],bounds[i-1,2],bounds[i-1,3]]
			else:
				bounds[i,:]=bounds[i-1,:]
			xinc+=1
		if ymoved==1:
			ymoved,xinc=0,0
							
	# Find the Centroid
	cent=[]
	for i in frames-1:
		aoi=pix[i,bounds[i,2]:bounds[i,3]+1,bounds[i,0]:bounds[i,1]+1]
		centx,centy=[],[]
		for j in range(bounds[i,2],bounds[i,3]+1):
			for k in range(bounds[i,0],bounds[i,1]+1):
				centx.append(aoi[j-bounds[i,2],k-bounds[i,0]]*(k+1))
				centy.append(aoi[j-bounds[i,2],k-bounds[i,0]]*(j+1))
		cent.append(((1./np.sum(aoi))*np.sum(centx),(1./np.sum(aoi))*np.sum(centy)))
	cent=np.array(cent)
	
	# Plot the frames desired
	if plt_rect==True:
		for frame in frames:
			plt.figure('int_area')
			plt.clf()
			im=plt.imshow(pix[(frame-1),:,:],cmap=colmap,interpolation='nearest')
			plt.plot(cent[frame-1][0]-1,cent[frame-1][1]-1,'ko',markersize=5)
			plt.colorbar(im,orientation='vertical')
			plt.title(run+': frame '+repr(frame))
			plt.xlabel('x Pixel')
			plt.ylabel('y Pixel')
			plt.axis([x[0],x[-1],y[0],y[-1]])
			ut.rect(bounds[frame-1,:])
			if xylim!='default':
				plt.axis(xylim)
			if zscale!='default':
				plt.clim(zscale)
			plt.show()
			if savepath!=False:
				ut.create_dir(savepath)
				if frame<=9:
					framestr='00'+str(frame)
				if 10<=frame<=99:
					framestr='0'+str(frame)
				if 100<=frame:
					framestr=str(frame)
				plt.savefig(savepath+'/centroid/im/'+str(date)+'_run'+str(filename)+'_'+framestr+'_cent_area.png')
	
	plt.figure('cent_all')
	plt.clf()
	
	# Plot each centroid position
	centxlist,centylist=[],[]
	for i in range(len(cent)):
		plt.figure('cent_solo')
		plt.clf()
		plt.plot(cent[i][0],cent[i][1],'ro')
		plt.xlabel('x pixel')
		plt.ylabel('y pixel')
		if i+1<=9:
			framestr='00'+str(i+1)
		if 10<=i+1<=99:
			framestr='0'+str(i+1)
		if 100<=i+1:
			framestr=str(i+1)
		plt.title('Centroid Movement '+run+'_fr'+framestr)
		plt.axis([1,np.shape(pix)[2]+1,1,np.shape(pix)[1]+1])
		plt.ylabel('y pixel')
		plt.xlabel('x pixel')
		if savepath!=False:
			plt.savefig(savepath+'centroid/centroid'+str(date)+'_run'+str(filename)+'_fr'+framestr+'.png')
		# Plot all centroid positions
		plt.figure('cent_all')
		plt.plot(cent[i][0],cent[i][1],'ro')
		plt.ylabel('y pixel')
		plt.xlabel('x pixel')
	if savepath!=False:
		plt.savefig(savepath+'centroid/'+str(date)+'_run'+str(filename)+'centroid_all.png')
	
	return cent	
	
#------------------------------------------------------------------------------------------
	
# Display the scan pattern and contour for bleaching
def scan_setup(shape_s,spacing_s,width_s,shape_b,spacing_b,width_b,wtype,**kwargs):
	
	# Define keyword arguments
	contours=kwargs.get('contours',[.9])
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	import os
	import funk
	
	# Make scan and bleach arrays array
	scan_array=funk.gaussarray(shape_s,spacing_s,width_s[0],width_s[1])
	bleach_array=funk.gaussarray(shape_b,spacing_b,width_b[0],width_b[1])
	
	# Plot the scan array and the contour for bleaching
	plt.figure('scan_setup')
	plt.clf()
	extent_s=[-(shape_s[0]-1)*spacing_s,(shape_s[0]-1)*spacing_s,-(shape_s[1]-1)*spacing_s,(shape_s[1]-1)*spacing_s]
	extent_b=[-(shape_b[0]-1)*spacing_b,(shape_b[0]-1)*spacing_b,-(shape_b[1]-1)*spacing_b,(shape_b[1]-1)*spacing_b]
	aspect_ratio_s=float(shape_s[0])/shape_s[1]
	aspect_ratio_b=float(shape_b[0])/shape_b[1]
	plt.imshow(scan_array,extent=extent_s,aspect=aspect_ratio_s)
	plt.colorbar(orientation='vertical')
	plt.xlabel('x position (micron)')
	plt.ylabel('y position (micron)')
	# Plot contour(s)
	contours_real=np.array(contours)*np.amax(bleach_array)
	cs=plt.contour(bleach_array,levels=contours_real,colors='firebrick',linewidths=3,linestyles=['dashed'],extent=extent_b,aspect=aspect_ratio_b)
	fmt={}
	for i,j in zip(cs.levels,contours):
		fmt[i]=str(j)
	plt.clabel(cs,inline=1,fmt=fmt,fontsize=10)
	plt.axis(extent_s)
	#offset=start_b-start_s
	#ut.textbox('Bleach Start Offset: ('+str(offset[0])+','+str(offset[1])+') micron',[.05,.98])

#------------------------------------------------------------------------------------------
	
# Display the bleaching expected
def bleach_setup(shape,spacing,width,**kwargs):
	
	# Define keyword arguments
	contours=kwargs.get('contours',[.9])
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	import os
	import funk
	
	# Make scan and bleach arrays array
	bleach_array,start=funk.gaussarray(shape,spacing,width[0],width[1],start=True)
	
	# Plot the scan array and the contour for bleaching
	plt.figure('scan_setup')
	plt.clf()
	extent_b=[-(shape[0]-1)*spacing,(shape[0]-1)*spacing,-(shape[1]-1)*spacing,(shape[1]-1)*spacing]
	aspect_ratio_b=float(shape[0])/shape[1]
	plt.imshow(bleach_array,extent=extent_b,aspect=aspect_ratio_b)
	plt.colorbar(orientation='vertical')
	plt.xlabel('x position (micron)')
	plt.ylabel('y position (micron)')
	# Plot contour(s)
	contours_real=np.array(contours)*np.amax(bleach_array)
	cs=plt.contour(bleach_array,levels=contours_real,colors='firebrick',linewidths=3,linestyles=['dashed'],extent=extent_b,aspect=aspect_ratio_b)
	fmt={}
	for i,j in zip(cs.levels,contours):
		fmt[i]=str(j)
	plt.clabel(cs,inline=1,fmt=fmt,fontsize=10)
	plt.axis(extent_b)

#------------------------------------------------------------------------------------------

def gauss_center(data,xfit,yfit,xint,yint,p0,**kwargs):
	
	savepath=kwargs.get('savepath',False)
	run=kwargs.get('run',0)
	fig=kwargs.get('fig','')
	highlander=kwargs.get('only','both')
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	import os

	if savepath!=False:
		ut.create_dir(savepath+'gauss_prof/')
		
	def gauss(x,amp,cen,width):
		return amp*np.exp(-2*np.square(x-cen)/(np.square(width)))
		
	# Integrate the xstrip and ystrip to get 1D arrays
	yprof=np.vstack((np.arange(yfit[0],yfit[1]+1),np.sum(data[yfit[0]-1:yfit[1],xint[0]-1:xint[1]],axis=1)))
	xprof=np.vstack((np.arange(xfit[0],xfit[1]+1),np.sum(data[yint[0]-1:yint[1],xfit[0]-1:xfit[1]],axis=0)))
	
	# Fit 1D gaussian to each cross section
	if highlander!='y':
		xparam=ut.fitter1d(gauss,xprof,[p0[0],p0[1],p0[3]],fig='xprof'+fig)
		ut.textbox('Center: '+str(round(xparam[1],2))+'\nWidth: '+str(round(xparam[2],2)),[.05,.95])
		if savepath!=False:
			plt.savefig(savepath+'gauss_prof/run'+str(run)+'xprof.png')
	if highlander!='x':
		yparam=ut.fitter1d(gauss,yprof,[p0[0],p0[2],p0[4]],fig='yprof'+fig)
		ut.textbox('Center: '+str(round(yparam[1],2))+'\nWidth: '+str(round(yparam[2],2)),[.05,.95])
		if savepath!=False:
			plt.savefig(savepath+'gauss_prof/run'+str(run)+'yprof.png')
	if highlander=='x':
		return xparam[1]
	if highlander=='y':
		return yparam[1]
		
	return xparam[1],yparam[1]
#------------------------------------------------------------------------------------------
def center_find(filename,fitbnds,**kwargs):

	# Define keyword arguments
	colmap=kwargs.get('cmap','jet')
	frames=kwargs.get('frames','all')
	store=kwargs.get('store',False)
	savepath=kwargs.get('savepath',False)
	units=kwargs.get('units','pixel')
	speed=kwargs.get('speed','slow')
	fit=kwargs.get('fit','gauss')
	plot_each=kwargs.get('plot_each',True)
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	import utility as ut
	import lv_analysis as lv
	from lmfit import minimize,Parameters,Parameter,Model
	
	subdir=fit+'_center/'
	
	if savepath!=False:
		ut.rm_dir(savepath+subdir)
		ut.create_dir(savepath+subdir)
		
	# Import data and remove cosmics with image_remove_cosmics
	x,y,pix,x_max,y_max,N,run=im_sub_dark(filename,speed=speed)
	# Slice data for just fit bounds
	pixslice=pix[:,fitbnds[2]:fitbnds[3]+1,fitbnds[0]:fitbnds[1]+1]
	x,y=np.arange(fitbnds[0],fitbnds[1]+1),np.arange(fitbnds[2],fitbnds[3]+1)
	
	if frames=='all':
		frames=range(1,N+1)
	
	# Define 2d gaussian
	def gauss2d_flat(x,y,amp,xcenter,ycenter,xwidth,ywidth):  # This ends up as f(y,x)
		my,mx=np.meshgrid(x,y)
		xwidth,ywidth=xwidth/2.,ywidth/2.	
		return np.ravel(amp*np.exp(-2*(np.square(mx-ycenter)/(np.square(ywidth))))*np.exp(-2*(np.square(my-xcenter))/(np.square(xwidth))))
	
	param_all,frames_fit=[],[]
	p0=[500,np.mean([fitbnds[0],fitbnds[1]]),np.mean([fitbnds[2],fitbnds[3]]),2,2]
	p_0s=['amp','xcenter','ycenter','xwidth','ywidth']
	for frame in frames:
		f=frame-1
		data=pixslice[f,:,:]
		if np.amax(data)>=3.5*np.mean(data):
			if fit=='gauss':
				# Define x,y grid and p_0 names
				func_model=Model(gauss2d_flat,independent_vars=['x','y'],param_names=p_0s)
				result=func_model.fit(np.ravel(data),x=x,y=y,amp=p0[0],xcenter=p0[1],ycenter=p0[2],xwidth=p0[3],ywidth=p0[4],verbose=False)
				param_all.append([result.best_values['xcenter'],result.best_values['ycenter']+1])
				frames_fit.append(frame)
			if fit=='centroid':
				# Find the Centroid
				bounds=np.array(fitbnds)-1
				centx,centy=[],[]
				for j in range(bounds[2],bounds[3]+1):
					for k in range(bounds[0],bounds[1]+1):
						centx.append(pix[f,j,k]*(k+1))
						centy.append(pix[f,j,k]*(j+1))
				param_all.append(((1./np.sum(pix[f,bounds[2]:bounds[3]+1,bounds[0]:bounds[1]+1]))*np.sum(centx),(1./np.sum(pix[f,bounds[2]:bounds[3]+1,bounds[0]:bounds[1]+1]))*np.sum(centy)))
				#centx,centy=[],[]
				#for l in range(fitbnds[2]-1,fitbnds[3]):
					#for k in range(fitbnds[0]-1,fitbnds[1]):
						#centx.append(pix[f,l,k]*(k+1))
						#centy.append(pix[f,l,k]*(l+1))
				#param_all.append([(1./np.sum(data))*np.sum(centx),(1./np.sum(data))*np.sum(centy)])
				frames_fit.append(frame)
	param_all=np.array(param_all)
	meanpos=np.mean(param_all,axis=0)
	if units=='micron':
		param_all=5*param_all
			
	# Plot the frames desired
	if plot_each==True:
		for j in range(np.shape(param_all)[0]):
			plt.figure(fit+'_fit')	
			plt.clf()
			im=plt.imshow(pix[frames_fit[j]-1,:,:],cmap=colmap,interpolation='nearest')
			plt.plot(param_all[j,0]-1,param_all[j,1]-1,'ko',markersize=4)
			plt.colorbar(im,orientation='vertical')
			plt.title(run+': frame '+repr(frames_fit[j]))
			if units=='pixel':
				plt.xlabel('x Pixel')
				plt.ylabel('y Pixel')
				plt.axis([fitbnds[0],fitbnds[1],fitbnds[2],fitbnds[3]])
			if units=='micron':
				plt.xlabel('x (micron)')
				plt.ylabel('y (micron)')
				plt.axis([fitbnds[0]*5,fitbnds[1]*5,fitbnds[2]*5,fitbnds[3]*5])
			ut.textbox('('+str(round(param_all[j,0],1))+','+str(round(param_all[j,1],1))+')',[.05,.95])
			plt.show()
			if savepath!=False:
				plt.savefig(savepath+subdir+str(date)+'_run'+str(filename)+'_'+str(frames_fit[j]).zfill(3)+'_'+fit+'_center.png')
	
	# Plot just the centers
	plt.figure('all_'+fit+'_centers')
	plt.clf()
	plt.plot(param_all[:,0],param_all[:,1],'bo')
	plt.title(str(date)+' run'+str(filename)+' Center for All Frames')
	if units=='pixel':
		plt.xlabel('x Pixel')
		plt.ylabel('y Pixel')
		#plt.axis([meanpos[0]-3,meanpos[0]+3,meanpos[1]-3,meanpos[1]-3])
	if units=='micron':
		plt.xlabel('x (micron)')
		plt.ylabel('y (micron)')
		plt.axis([fitbnds[0]*5,fitbnds[1]*5,fitbnds[2]*5,fitbnds[3]*5])
	#ut.textbox('Fit Bounds: '+repr(fitbnds),[.05,.95])
	if savepath!=False:
		plt.savefig(savepath+subdir+str(date)+'_run'+str(filename)+'_'+fit+'_center_all.png')
		ut.save(np.transpose(param_all),['x pixel','ypixel'],savepath+subdir+str(date)+'_run'+str(filename)+'_'+fit+'_center_all.txt')
		
	return param_all,frames_fit

#------------------------------------------------------------------------------------------
# create linearity plots for ion deposits
def linearity(l_area,T_R,**kwargs):
	
	xwidth=kwargs.get('xwidth',3)
	ywidth=kwargs.get('ywidth',3)
	savepath=kwargs.get('savepath',False)
	A_obs=kwargs.get('A_obs',8.6)
	axis=kwargs.get('axis','default')
	bndfind=kwargs.get('chgbnds','auto')
	vib=kwargs.get('vib_shutter',True)
	int_area=kwargs.get('int_area','manual')
	cutoff=kwargs.get('cutoff',5.)
	bndtype=kwargs.get('boundtype','edge')
	fit0=kwargs.get('fit0',True)
	charge=kwargs.get('charge','cup')
	drops=kwargs.get('runs2skip',[])
	ionlim=kwargs.get('ion_limit','all')
	
	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	import utility as ut
	import lv_analysis as lv
	import gen_reader as gr
	import scope_cal as scp
	import scipy.special as sp
	
	if savepath!=False:
		ut.create_dir(savepath)
		savescope=savepath+'scope/'
		ut.create_dir(savescope)
	else:
		savescope=savepath
	
	# make l_area the edges of integral area
	if bndtype=='center':
		l_area=[l_area[0]-.5*l_area[2],l_area[0]+.5*l_area[2],l_area[1]-.5*l_area[3],l_area[1]+.5*l_area[3]]	
		
	# import rundata.csv to set up lists
	runlist,lvlist,scplist,pulselist=ut.data_lists(date)
	
	# be able to ignore runs
	drops_i=[]
	for i in drops:
		drops_i.append(np.where(runlist==i)[0][0])
	runlist,lvlist,scplist,pulselist=np.delete(runlist,drops_i),np.delete(lvlist,drops_i),np.delete(scplist,drops_i,axis=1),np.delete(pulselist,drops_i)
	
	# calculate charge correction for each pulse
	pulsecorr=.802*sp.erf(.007*pulselist)
	
	# calculate charge/pulse for cup3 with correction
	fC=pulsecorr*scp.pulse_charge_lin(scplist,null_pulse,savepath=savescope,cutoff=7.,chan=3,bounds=bndfind)
	
	# integrate signal,scale by mW*s and calculate number of ions
	k=0
	ions,ints,labels,xe,zp,energy,pos=[],[],[],[],[],[],[]
	for i in runlist:
		labels.append('run'+str(i))
		# Get laser energy (power * exposure time)
		if vib==True:
			energy.append(1000*lv.lv_energy(i,T_R=T_R)[0])
		if vib=='old':
			energy.append(np.mean(lv.lv_pow(lvlist[np.where(runlist==i)[0][0]],T_R=T_R))*np.mean(lv.lv_time(lvlist[np.where(runlist==i)[0][0]])))
			# old way of calculating energy
		if vib==False:
			energy.append(1E3)
		# Do the integral
		if int_area=='auto':
			# Find center of laser
			param,cov,fr=gauss2d_fit(i,l_area,savepath=savepath)
			pos.append(param[0,:]+1)
			if savepath!=False:
				integral=im_int_0(i,[param[0,0],param[0,1],xwidth,ywidth],boundtype='center',savepath=savepath+'ints/',overlay=True,axis=l_area)
			if savepath==False:
				integral=im_int_0(i,[param[0,0],param[0,1],xwidth,ywidth],boundtype='center',overlay=True,axis=l_area)
		if int_area=='manual':
			pos.append([np.mean(l_area[:2]),np.mean(l_area[2:])])
			if savepath!=False:
				integral=im_int_0(i,l_area,boundtype='edge',savepath=savepath+'ints/',overlay=True,axis=[l_area[0]-5,l_area[1]+5,l_area[2]-5,l_area[3]+5])
			if savepath==False:
				integral=im_int_0(i,l_area,boundtype='edge',overlay=True,axis=[l_area[0]-5,l_area[1]+5,l_area[2]-5,l_area[3]+5])
			xwidth,ywidth=l_area[1]-l_area[0],l_area[3]-l_area[2]
		ints.append(integral[0])
		# Find number of ions from pulse data
		if pulselist[k]==-1:
			ions.append(0.)
			xe.append(k)
		if pulselist[k]==0:
			ions.append(0.)
			zp.append(k)
		if pulselist[k] not in [0,-1]:
			ions.append(.8*pulselist[k]*(fC[1,k]/1.602E-4)*(A_obs/(np.pi*np.square(1.4E3)))*.5)
		k+=1
	ions=np.array(ions)
	ints=np.array(ints)/np.array(energy)
	pos=np.array(pos)
	
	# option to fit only lower ion runs
	low=[]
	if ionlim!='all':
		for i in range(len(ions)):
			if ions[i]<=ionlim:
				low.append(i)
		ions,ints,pos,runlist,pulselist,fC=np.take(ions,low),np.take(ints,low),np.take(pos,low,axis=0),np.take(runlist,low),np.take(pulselist,low),np.take(fC,low)
		# Re-assign indicies to xe and zero-pulse runs
		xe,zp=[],[]
		for k in range(len(pulselist)):
			if pulselist[k]==-1:
				xe.append(k)
			if pulselist[k]==0:
				zp.append(k)
			
	# Subtract average of surrounding Xe BG runs from pulse runs 
	ints_sub=np.zeros(np.shape(ints))
	for i in range(len(runlist)):
		if pulselist[i]!=-1:
			for j in range(len(runlist)-(i+1)):
				if pulselist[i+(j+1)]==-1:
					xe_after=ints[i+(j+1)]
					break
			for k in range(i):
				if pulselist[i-(k+1)]==-1:
					xe_before=ints[i-(k+1)]
					break
			ints_sub[i]=ints[i]-np.mean([xe_before,xe_after])
	ints_sub_0=np.take(ints_sub,zp)
	ions_sub_0=np.take(ions,zp)
	
	# Save the data in ascii
	if savepath!=False:
		ut.save(np.vstack((runlist,ions,ints,fC,pos[:,0],pos[:,1])),['run','ions','counts','fC/pulse','xpos','ypos'],savepath+'rundata.txt')
		ut.save(np.delete(np.vstack((runlist,ions,ints_sub,fC,pos[:,0],pos[:,1])),xe,axis=1),['run','ions','counts(sub)','fC/pulse','xpos','ypos'],savepath+'rundata_sub.txt')
		ut.save(np.take(np.vstack((runlist,ions,ints,fC,pos[:,0],pos[:,1])),xe,axis=1),['run','ions','counts','fC/pulse','xpos','ypos'],savepath+'bg_data.txt')
		
	#Fit line to data
	if fit0==False:
		params_line=ut.fitter1d(funk.poly,np.vstack((np.delete(ions,xe),np.delete(ints_sub,xe))),[0,200])
		fitx=np.linspace(-1,1.1*np.amax(np.delete(ions,xe)),20)
		fit=funk.poly(np.linspace(-1,1.1*np.amax(np.delete(ions,xe)),20),params_line[0],params_line[1])
	if fit0==True:
		params_line=ut.fitter1d(funk.line0,np.vstack((np.delete(ions,xe),np.delete(ints_sub,xe))),[200])
		fitx=np.linspace(-1,1.1*np.amax(np.delete(ions,xe)),20)
		fit=funk.line0(np.linspace(-1,1.1*np.amax(np.delete(ions,xe)),20),params_line[0])
		
	# Make the Linearity Plot
	plt.figure('linearity')
	plt.clf()
	plt.plot(np.delete(ions,xe),np.delete(ints,xe),'bo')
	plt.plot(np.take(ions,xe),np.take(ints,xe),'ro')
	plt.plot(np.take(ions,zp),np.take(ints,zp),'go')
	plt.plot(fitx,fit+np.mean(np.take(ints,xe)),'r--')
	plt.xlabel('Number of Ions in Laser Region')
	plt.title(str(date)+ ' Linearity Plot',fontsize=12)
	plt.ylabel('Integrated Counts per mW in First Frame ('+str(xwidth)+'x'+str(ywidth)+' Region)')
	plt.legend(['Ba Runs','Xe Runs','0-pulse Runs'],numpoints=1,fontsize=12,loc=2)
	for label,x,y in zip(labels,ions,ints):
		plt.annotate(label,xy=(x,y),xytext=(0,5),textcoords='offset points',ha = 'left',va ='bottom')
	plt.xlim(-1,1.25*np.amax(ions))
	if np.amax(ints)>=1E3:
		ut.sci()
	plt.show()
	if savepath!=False:
		plt.savefig(savepath+'linearity_plot.png')
	if ionlim=='all':
		plt.axis([-1,10,0,1.25*fit[ut.bound_finder(fitx,[10])[0]]+np.mean(np.take(ints,xe))])
		if savepath!=False:
			plt.savefig(savepath+'linearity_plot_small.png')
	
	# Make the Linearity Plot for Subtracted Data
	plt.figure('linearity sub')
	plt.clf()
	plt.plot(np.delete(ions,xe),np.delete(ints_sub,xe),'bo')
	plt.plot(ions_sub_0,ints_sub_0,'go')
	plt.plot(fitx,fit,'r--')
	plt.xlabel('Number of Ions in Laser Region')
	plt.title(str(date)+ ' Linearity Plot Subtracted',fontsize=12)
	plt.ylabel('Integrated Counts per mW in First Frame ('+str(xwidth)+'x'+str(ywidth)+' Region)')
	ut.textbox('Counts/Ion: '+str(round(params_line[0],2)),[.05,.95])
	for label,x,y in zip(np.delete(labels,xe),np.delete(ions,xe),np.delete(ints_sub,xe)):
		plt.annotate(label,xy=(x,y),xytext=(5,5),textcoords='offset points',ha = 'right',va ='bottom')
	ut.zero(bounds=[-1,1.2*np.amax(ions)])
	if np.amax(ints_sub)>=1E3:
		ut.sci()	
	plt.show()
	if ionlim=='all':
		plt.axis([-1,10,-200,1.25*fit[ut.bound_finder(fitx,[10])[0]]])
		if savepath!=False:
			plt.savefig(savepath+'linearity_plot_sub_small.png')
	plt.axis([-1,1.05*np.amax(ions),-200,1.25*np.amax(ints_sub)])
	if savepath!=False:
		plt.savefig(savepath+'linearity_plot_sub.png')
	
	# Plot the Xe runs and statistics
	xe_ave,xe_std=np.mean(np.take(ints,xe)),np.std(np.take(ints,xe))
	plt.figure('xe runs')
	plt.clf()
	plt.plot(range(1,len(xe)+1),np.take(ints,xe),'bo')
	plt.xlim(0,len(xe)+2)
	ut.hline(xe_ave,legend=False)
	plt.xlabel('Xe Run')
	plt.ylabel('Integrated Counts per mW in First Frame ('+str(xwidth)+'x'+str(ywidth)+' Region)')
	plt.title(str(date)+' Counts per mW of Xe-only runs')
	for label,x,y in zip(np.take(labels,xe),range(1,len(xe)+1),np.take(ints,xe)):
		plt.annotate(label,xy=(x,y),xytext=(5,5),textcoords='offset points',ha = 'right',va ='bottom')
	ut.textbox('Average: '+str(round(xe_ave,1))+'\nStDev: '+str(round(xe_std,1)),[.05,.95])	
	if savepath!=False:
		plt.savefig(savepath+'xe_only_runs.png')
		
	return ions,ints,ints_sub,xe,zp,runlist
#------------------------------------------------------------------------------------------

def laser_overlay(positions,w,**kwargs):
	
	store=kwargs.get('store',True)
	fit_gauss=kwargs.get('fit_gauss',False)
	savepath=kwargs.get('savepath',False)
	savefile=kwargs.get('savefile',False)
	
	import matplotlib as mpl
	import numpy as np
	import utility as ut
	import funk
	
	mean_pos=[np.mean(positions[:,0]),np.mean(positions[:,1])]
	range_pos=[np.absolute(np.amax(positions[:,0])-np.amin(positions[:,0])),np.absolute(np.amax(positions[:,1])-np.amin(positions[:,1]))]
	x=np.linspace(np.amin(positions[:,0])-1.5*w[0],np.amax(positions[:,0])+1.5*w[0],500)
	y=np.linspace(np.amin(positions[:,1])-1.5*w[1],np.amax(positions[:,1])+1.5*w[1],500)
	laser_cov=np.zeros((500,500))
	for i in range(np.shape(positions)[0]):
		laser_cov+=funk.gauss2dw(x,y,1.,positions[i,0],positions[i,1],w[0],w[1])
	
	if fit_gauss==True:
		params=ut.gaussfit2d(laser_cov,x,y,[10,mean_pos[0],mean_pos[1],w[0],w[1]])
		
	area=ut.e_fold_area(laser_cov,x,y)
	
	plt.figure('laser')
	plt.clf()
	plt.imshow(np.flipud(laser_cov),extent=[x[0],x[-1],y[0],y[-1]])
	plt.colorbar(orientation='vertical')
	plt.xlabel('x (micron)')
	plt.ylabel('y (micron)')
	ut.textbox('w_x: '+str(w[0])+'micron\nw_y: '+str(w[1])+'micron\nArea: '+str(round(area,2))+'sq.micron',[.05,.25],fontsize=10)
	plt.show()
	if savepath!=False:
		ut.create_dir(savepath)
		plt.savefig(savepath+'laser_area.png')
	if savefile!=False:
		plt.savefig(savefile)
		
	if store==True:
		if fit_gauss==True:
			return laser_cov,area,params
		if fit_gauss==False:
			return laser_cov,area
	
#------------------------------------------------------------------------------------------

def gauss2d_fit(filename,fitbnds,**kwargs):
	
	frames=kwargs.get('frames','all')
	speed=kwargs.get('speed','slow')
	
	import matplotlib as mpl
	import numpy as np
	import utility as ut
	import funk
	from lmfit import minimize,Parameters,Parameter,Model
	
	# Import data and remove cosmics with image_remove_cosmics
	x,y,pix,x_max,y_max,N,run=im_sub_dark(filename,speed=speed)
	# Slice data for just fit bounds
	data=pix[:,fitbnds[2]:fitbnds[3]+1,fitbnds[0]:fitbnds[1]+1]
	x,y=np.arange(fitbnds[0],fitbnds[1]+1),np.arange(fitbnds[2],fitbnds[3]+1)
	
	if frames=='all':
		frames=range(1,N+1)
	
	# Define 2d gaussian
	def gauss2d_flat(x,y,amp,xcenter,ycenter,xwidth,ywidth):  # This ends up as f(y,x)
		my,mx=np.meshgrid(x,y)
		xwidth,ywidth=xwidth/2.,ywidth/2.	
		return np.ravel(amp*np.exp(-2*(np.square(mx-ycenter)/(np.square(ywidth))))*np.exp(-2*(np.square(my-xcenter))/(np.square(xwidth))))
	
	param_all,frames_fit=[],[]
	p0=[500,np.mean([fitbnds[0],fitbnds[1]]),np.mean([fitbnds[2],fitbnds[3]]),2,2]
	p_0s=['amp','xcenter','ycenter','xwidth','ywidth']
	# Define x,y grid and p_0 names
	for i in frames:
		func_model=Model(gauss2d_flat,independent_vars=['x','y'],param_names=p_0s)
		result=func_model.fit(np.ravel(data),x=x,y=y,amp=p0[0],xcenter=p0[1],ycenter=p0[2],xwidth=p0[3],ywidth=p0[4],verbose=False)
		param_all.append([result.best_values['amp'],result.best_values['xcenter']+1,result.best_values['ycenter']+1,result.best_values['xwidth'],result.best_values['ywidth']])
	
	param_all=np.array(param_all)
	
	return param_all
	
#------------------------------------------------------------------------------------------

def scan_corr(scan1,scan2,rng):
	
	import numpy as np
	import matplotlib as mpl
	import scipy.optimize as opt
	
	corr,shift=[],[]
	for i in np.arange(-rng[0],rng[0]+1):
		for j in np.arange(-rng[1],rng[1]+1):
			corr.append(np.sum((scan1-np.mean(scan1))*(ut.int_shift(scan2,i,j,val=0)-np.mean(ut.int_shift(scan2,i,j,val=0))))/((np.shape(scan1)[0]-abs(i))*(np.shape(scan1)[1]-abs(j))))
			shift.append([i,j])
	corr=np.transpose(np.array(corr).reshape(2*rng[1]+1,2*rng[0]+1))
	plt.figure('correlation')
	plt.imshow(corr,interpolation='nearest',origin='lower',extent=[-rng[0],rng[0],-rng[1],rng[1]])
	return corr,shift
			
#------------------------------------------------------------------------------------------

def scan_chi2(scan1,scan2,rng):
	
	import numpy as np
	import matplotlib as mpl
	import scipy.optimize as opt
	
	corr,shift=[],[]
	for i in np.arange(-rng[0],rng[0]+1):
		for j in np.arange(-rng[1],rng[1]+1):
			corr.append(np.sum(np.square(scan1-ut.int_shift(scan2,i,j,val=0)))/((np.shape(scan1)[0]-abs(i))*(np.shape(scan1)[1]-abs(j))))
			shift.append([i,j])
	corr=np.transpose(np.array(corr).reshape(2*rng[1]+1,2*rng[0]+1))
	plt.figure('correlation')
	plt.imshow(corr,interpolation='nearest',origin='lower',extent=[-rng[0],rng[0],-rng[1],rng[1]])
	return corr,shift
			
