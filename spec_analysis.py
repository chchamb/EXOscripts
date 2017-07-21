# Script For handling spectra analysis

date=20160603
darkrun=0
#=========================================================================================

# Function to import data from ascii files
def spec_import(filenum):

	# Import modules for array and plotting
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	import utility as ut
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
	dat.columns=['x','ROI','f','px']
	# Reshape into [f,ROI,x] ndarray
	y_max,N=np.amax(dat.ROI),np.amax(dat.f)
	x_max=int(len(dat.px)/(y_max*N))
	x,ROI=np.arange(x_max),np.arange(y_max)
	pix=np.reshape(dat.px,(N,y_max,x_max))
				
	return x,y,pix,x_max,y_max,N,run

#-----------------------------------------------------------------------------------------
# Function to remove and display location of cosmics
def spec_cosmic_remove(filenum,ROI,f_rem):
	
	# Import numpy and pyplot
	import numpy as np
	import os

	# Import data from file
	expnum=os.listdir('/home/chris/anaconda/data/'+str(date)+'/raw/ascii')[1][:6]
	filename='/home/chris/anaconda/data/'+str(date)+'/raw/ascii/'+expnum+str(filenum)+'.txt'
	datafile= open(filename,'r')
	
	# Get run number from data file name
	run=''
	j=1                                                        # j is the incrementor for characters in the run name/number
	for i in range(1,len(filename)+1):
		if filename[-i]=='/':                                      # Steps backwards in path until it finds a '/'
			while filename[-i+j]!='.':
				run += filename[-i+j]                              # Adds on characters after the '/' until it finds a '.'
				j+=1
			break
		else:
			continue
			
	# Initialize lists
	column1=[]
	column4=[]
	for line in datafile:
		line=line.strip()
		values=line.split()
		column1.append(float(values[0]))
		column4.append(float(values[3]))
		M= int(values[1])
		N= int(values[2])
	# Close the datafile, what were you, raised in MATLAB?
	datafile.close()
	
	# Get the number of Wavelengths and ROI's pixels (usually 1340 x 1, but could change)
	wave_max= int(len(column1)/(N*M))
	wavelength=column1[:wave_max]
	wavelength= np.array(wavelength)
	
	# Calibrate the wavelength
	if date==20140618:
		wavelength_s=wavelength-500
		def dwavelength(x):
			return 2.22E-3*np.square(x)-.303*x+8.21
		wave_corr=dwavelength(wavelength_s)	
		wavelength=wavelength_s+wave_corr+500
	
	if date==20140916 and filenum>=103:
		def dwavelength(x):
			return 6.4432E-5*np.square(x)-.063483*x+18.481
		wavelength=wavelength+dwavelength(wavelength)

	
	# Organize data into a numpy array (frame,y,x).
	counts=np.zeros((N,M,wave_max),dtype=float)
	for i in range(N):
		for j in range(M):  
			# This grabs a chunk of length wave_max from the unorganized list (one ROI)
			counts[i,j,:]=column4[(i*M+j)*wave_max:(i*M+j+1)*wave_max]
	
	# Remove Frames
	if f_rem!=False:
		if len(f_rem)>1:
			f_rem=np.flipud(np.sort(np.array(f_rem)))
		for f in f_rem:
			counts=np.delete(counts,(f-1),0)
		N=N-len(f_rem)	
	# Convert to numpy arrays
	counts_cf= np.array(counts)
	
	# Remove cosmic rays by comparing points to average of surrounding points
	for k in range(N):
		for j in range(M):
			reps= 0
			rep_waves=[]
			for i in range(wave_max): 
				if i <= 2:
					continue
				if i >= (wave_max-3):
					continue
				else:
					bin= np.concatenate([counts[k,j,i-2:i],counts[k,j,i+1:i+3]])
				bin_ave= np.mean(bin)
				if counts[k,j,i]>=(1.5*bin_ave):
					rep_vals= np.concatenate([bin,[bin_ave]])
					counts_cf[k,j,i-1:i+2]= np.amin(rep_vals)
					reps= reps+1
					rep_waves.append(round(wavelength[i],2))
	
	return wavelength,counts_cf[:,(ROI-1),:],N,M,run
	
#-----------------------------------------------------------------------------------------

# Function to average and subtract background from signals
def spec_sub(sigpath,bgpath,ROI,output,darkframe,bgscale,f_rem1,f_rem2):
	 
	# Import Modules
	import numpy as np
	import matplotlib.pyplot as plt
	 
	# Import data and remove cosmics from background and signal with spec_cosmic_remove
	wavelength,sig_counts_cf,sig_N,sig_M,sig_run=spec_cosmic_remove(sigpath,'signal',output,darkframe,ROI,f_rem1)
	wavelength,bg_counts_cf,bg_N,bg_M,bg_run=spec_cosmic_remove(bgpath,'background',output,darkframe,ROI,f_rem2)
	 
	# Average the backgrounds
	bg_ave=bgscale*np.mean(bg_counts_cf,axis=0)
		
	# Subtract averaged background from each signal frame
	counts_sub=[]
	for i in range(sig_N):
		counts_sub.append(sig_counts_cf[i,:]-bg_ave)
	counts_sub=np.array(counts_sub)
	
	return wavelength,counts_sub,sig_N,sig_M,sig_run,bg_run
	
#-----------------------------------------------------------------------------------------

# Function to average and subtract background from signals
def spec_sub_dark(sigpath):
	 
	# Import Modules
	import numpy as np
	import matplotlib.pyplot as plt
	 
	# Import data and remove cosmics from background and signal with spec_cosmic_remove
	wavelength,sig_counts_cf,sig_N,sig_M,sig_run=spec_import(sigpath)
	if darkrun!=0:
		wavelength,dk_counts_cf,dk_N,dk_M,dk_run=spec_import(darkrun)
		# Average the backgrounds
		dk_ave=np.mean(dk_counts_cf,axis=0)
	if darkrun==0:
		dk_ave=np.zeros(np.shape(sig_counts_cf))
		
	# Subtract averaged background from each signal frame
	counts_sub=[]
	for i in range(sig_N):
		counts_sub.append(sig_counts_cf[i,:]-dk_ave)
	counts_sub=np.array(counts_sub)
	
	return wavelength,counts_sub,sig_N,sig_M,sig_run
	
#-----------------------------------------------------------------------------------------
	
# Function to average and subtract background from signals
def spec_ave(filepath,col,**kwargs):

	# Define keyword arguments
	ROI=kwargs.get('ROI',1)
	output=kwargs.get('disp_cosmics',False)
	darkframe=kwargs.get('dark_frame',False)
	scale=kwargs.get('scale',1.)
	bgscale=kwargs.get('bg_scale',1)
	store=kwargs.get('store',False)
	f_rem=kwargs.get('remove_frames',False)
	bgfile=kwargs.get('bgfile',False)
	
	 
	# Import Modules
	import numpy as np
	import matplotlib.pyplot as plt
	 
	if bgfile==False:
		# Remove cosmic rays from data with spec_cosmic_remove
		wavelength,counts,N,M,run=spec_cosmic_remove(filepath,'spectrum',False,darkframe,ROI,f_rem)
		lab=run+' -- ave'
	
	else:
		# Import data and remove cosmics from signal and subtract background with spec_sub
		if f_rem==False:
			f_rem=[False,False]
		wavelength,counts,N,M,run,b_run=spec_sub(filepath,bgfile,ROI,False,darkframe,bgscale,f_rem[0],f_rem[1])
		lab=run+' - '+repr(bgscale)+'*'+b_run+' -- aves'
	
	# Average the counts
	counts_ave=np.mean(counts,axis=0)

	# Plot on whatever the current axis is
	counts_ave=scale*counts_ave
	plt.plot(wavelength,counts_ave,color=col,label=lab)	
	plt.show()
	
	if store==True:
		return wavelength,counts_ave
	
#-----------------------------------------------------------------------------------------
	
# Function to plot spectra
def spec_overlay(filepath,frame,col,**kwargs):
	
	# Define keyword arguments
	ROI=kwargs.get('ROI',1)
	dc=kwargs.get('disp_cosmics',False)
	bgfile=kwargs.get('bg_file',False)
	darkframe=kwargs.get('dark_frame',False)
	scale=kwargs.get('scale',1)
	bgscale=kwargs.get('bg_scale',1.)
	store=kwargs.get('store',False)
	tag=kwargs.get('tag','')
	f_rem=kwargs.get('remove_frames',False) # [[sig frames to remove],[bg frames to remove]]
	
	# Import numpy and pyplot
	import numpy as np
	import matplotlib.pyplot as plt

	if bgfile==False:
		# Remove cosmic rays from data with spec_cosmic_remove
		wavelength,counts,N,M,run=spec_cosmic_remove(filepath,'spectrum',dc,darkframe,ROI,f_rem)
		if frame==0:
			frame=N
		lab=run+'.'+repr(frame)
	
	else:
		# Import data and remove cosmics from signal and subtract background with spec_sub
		if f_rem==False:
			f_rem=[False,False]
		wavelength,counts,N,M,run,b_run=spec_sub(filepath,bgfile,ROI,dc,darkframe,bgscale,f_rem[0],f_rem[1])
		if frame==0:
			frame=N
		lab=run+'.'+repr(frame)+' - '+b_run
		
	# Plot on whatever the current axis is
	counts=counts*scale
	ax=plt.gca()
	lines=ax.plot(wavelength,counts[(frame-1),:],color=col,label=lab+' '+tag)
	plt.draw()
	plt.show()
	
	if store==True:
		return wavelength,counts[(frame-1),:]
	
#-----------------------------------------------------------------------------------------

# Function to plot spectra in 2-D
def spec_2D(filepath,**kwargs):
	
	# Define keyword arguments
	ROI=kwargs.get('ROI',1)
	dc=kwargs.get('disp_cosmics',False)
	bgfile=kwargs.get('bgfile',False)
	df=kwargs.get('diff_frames',False)
	darkframe=kwargs.get('dark_frame',False)
	bgscale=kwargs.get('bg_scale',1)
	col=kwargs.get('cmap','jet')
	f_rem=kwargs.get('remove_frames',False)
	
	# Import numpy and pyplot
	import numpy as np
	import matplotlib.pyplot as plt

	if bgfile==False:
		# Remove cosmic rays from data with spec_cosmic_remove
		wavelength,counts,N,M,run=spec_cosmic_remove(filepath,'spectrum',dc,darkframe,ROI,f_rem)
	
	else:
		# Import data and remove cosmics from signal and subtract background with spec_sub
		if f_rem==False:
			f_rem=[False,False]
		wavelength,counts,N,M,run,b_run=spec_sub(filepath,bgfile,ROI,dc,darkframe,bgscale,f_rem[0],f_rem[1])
		
	# Set the axes of the plot
	extent=[wavelength[0],wavelength[-1],1,N]
	#a_ratio=(wavelength[-1]-wavelength[0])/float(len(wavelength))
	# Find indices for bounds
	bound_index=[]
	bounds=[560,610]
	for k in range(2):
		j=-1
		for i in wavelength:
			j=j+1
			if i<bounds[k]:
				last_wave=i
			else:
				a=abs(bounds[k]-last_wave)
				b=abs(i-bounds[k])
				if a>b:
					bound_index.append(j)
					break
				else:
					bound_index.append(j-1)
					break
		
	# Plot on whatever the current axis is
	counts=np.flipud(counts*scale)
	fig=plt.figure(run)
	fig.clf()
	im=plt.imshow(counts,extent=extent,cmap=col,interpolation='nearest',vmin=0,vmax=np.amax(counts[:,bound_index[0]:bound_index[1]]))
	plt.colorbar(im,orientation='vertical')
	plt.title(run)
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Frame')
	plt.axis([wavelength[0],wavelength[-1],1,N])
	plt.show()
	
#-----------------------------------------------------------------------------------------
			
# Function to plot spectra
def spec_plot(filepath,**kwargs):
	
	# Define keyword arguments
	ROI=kwargs.get('ROI',1)
	frames=kwargs.get('frames','all')
	dc=kwargs.get('disp_cosmics',False)
	bgfile=kwargs.get('bgfile',False)
	df=kwargs.get('diff_frames',False)
	overlay=kwargs.get('overlay',False)
	darkframe=kwargs.get('dark_frame',False)
	bgscale=kwargs.get('bg_scale',1)
	store=kwargs.get('store',False)
	savepath=kwargs.get('savepath',False)
	f_rem=kwargs.get('remove_frames',False) # [[sig frames to remove],[bg frames to remove]]
	yr=kwargs.get('ylim','default')
	
	# Import numpy and pyplot
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut

	if bgfile==False:
		# Remove cosmic rays from data with spec_cosmic_remove
		wavelength,counts,N,M,run=spec_sub_dark(filepath)
	
	else:
		# Import data and remove cosmics from signal and subtract background with spec_sub
		if f_rem==False:
			f_rem=[False,False]	
		wavelength,counts,N,M,run,b_run=spec_sub(filepath,bgfile,ROI,dc,darkframe,bgscale,f_rem[0],f_rem[1])
		
	if frames=='all':
		frames=range(1,N+1)
				
	# Plot the frames desired
	if overlay==True:
		for frame in frames:
			plt.figure('spec_plot')
			plt.clf()
			plt.plot(wavelength,counts[(frame-1),:])
			plt.title(run+': frame '+repr(frame))
			plt.xlabel('Wavelength (nm)')
			plt.ylabel('Counts')
			plt.axis([wavelength[0],wavelength[-1],0,1.1*np.amax(counts[(frame-1),:])])
			if yr!='default':
				plt.ylim(yr[0],yr[1])
			if savepath!=False:
				ut.create_dir(savepath)
				plt.savefig(savepath+str(date)+'_run'+str(run[7:])+'_f'+str(frame)+'.png')
			plt.show()
	if overlay==False:
		fig=plt.figure(run)
		fig.clf()
		counts_max=[]
		for frame in frames:
			plt.clf()
			plt.plot(wavelength,counts[(frame-1),:],label='Frame '+repr(frame))
			plt.title(run)
			plt.xlabel('Wavelength (nm)')
			plt.ylabel('Counts')
			plt.axis([np.amin(wavelength),np.amax(wavelength),0,1.1*np.amax(counts[frame-1,:])])
			plt.legend()
			plt.show()
			if savepath != False:
				ut.create_dir(savepath)
				plt.savefig(savepath+str(date)+'_run'+str(run[7:])+'_f'+str(frame)+'.png')
	
	if store==True:
		return wavelength,counts

#-----------------------------------------------------------------------------------------

# Function to store spectra
def spec_store(filepath,**kwargs):
	
	# Define keyword arguments
	ROI=kwargs.get('ROI',1)
	date_2=kwargs.get('date',date)
	
	# Import numpy and pyplot
	import numpy as np
	import matplotlib.pyplot as plt
	import os

	# Remove cosmic rays from data with spec_cosmic_remove
	wavelength,counts,N,M,run=spec_cosmic_remove(filepath,ROI,[])
		
	return wavelength,counts

#-----------------------------------------------------------------------------------------
# Excitation spectra for discrete steps (laser not in spectrum)
def disc_ex(runlist,bounds,lvnum,step_info,**kwargs):
	
	savepath=kwargs.get('savepath',False)
	p0=kwargs.get('p0','int')
	store=kwargs.get('store',False)
	
	import numpy as np
	import matplotlib.pyplot as plt
	import gen_reader as gr
	import utility as ut
	
	def solo_gauss(x,b,a,c):
		return a*np.exp(-(4*np.log(2)*np.square(x-b))/(np.square(c)))
	def solo_lorentz(x,b,a,c):
		return a*(1./(1+np.square((x-b)/(c/2))))
		
	lvarr=gr.reader('/home/chris/anaconda/data/'+str(date)+'/lv/'+str(lvnum)+'_laser_temp.txt',header=False,delimeter=';')
	col=['b','r','g','c']
	
	# Create excitation wavelength array
	wlength_ex_up=np.arange(step_info[0],step_info[0]+(step_info[1]*(len(runlist)/2+1)),step_info[1])
	wlength_ex_dn=np.flipud(wlength_ex_up[:-1])
	wlength_ex=np.concatenate((wlength_ex_up,wlength_ex_dn))
	wlength_ex=157.3*wlength_ex+459.46
	
	# Powerscale and integrate the spectra
	integral,param_all=[],[]
	for i in runlist:
		wlength,cts=spec_plot(i,ROI=1,store=True,overlay=True)
		bnds=ut.bound_finder(wlength,bounds)
		cts_sc=cts/(1E3*lvarr[runlist.index(i),1])
		if p0=='int':
			ut.vline(wlength[bnds[0]])
			ut.vline(wlength[bnds[1]])
			integral.append(np.sum(cts_sc[:,bnds[0]:bnds[1]],axis=1))
		else:
			params,fits,types=spec_fit_fast(wlength[bnds[0]:bnds[1]],cts_sc[bnds[0]:bnds[1]],p0,'none','sub')
			params=params[0,:,:]
			param_all.append(params)
			for j in range(len(p0)):
				plt.plot(wlength,solo_gauss(wlength,params[j,0],params[j,2],params[j,1]),col[j],label=str(round(params[j,0],2))+'nm')
			plt.plot(wlength,fits[0,:],'k--',label='Full Fit')
			plt.legend()
		if savepath!=False:
			ut.create_dir(savepath)
			#plt.savefig(savepath+'excitation_spectrum.png')
	if p0!='int':
		param_all=np.array(param_all)
	
	# Make a Plot
	plt.figure('disc_ex_spec')
	plt.clf()
	if p0=='int':
		plt.plot(wlength_ex,np.array(integral),'bo')
	if p0!='int':
		for i in range(np.shape(param_all)[1]):
			plt.plot(wlength_ex,param_all[:,i,2],col[i],marker='o',ls=None,label=str(round(param_all[0,i,0],2))+'nm')	
	plt.xlabel('Excitation Wavelength (nm)')
	plt.ylabel('Integrated Counts/mW')
	plt.title(str(date)+' Excitation Spectrum')
	plt.legend()
	if savepath!=False:
		plt.savefig(savepath+str(date)+'excitation_spectrum.png')
	
	if store==True:
		if p0=='int':
			return wlength_ex,np.array(integral)

#-----------------------------------------------------------------------------------------


def excite_bg(signum,lvnumsig,bgnum,param_0,savefile,**kwargs):
	
	import the_one_script as tos
	import gen_reader as gr
	import matplotlib.pyplot as plt
	import numpy as np
	import scipy.special as sp
	import math
	import os
	import shutil
	import asciitable as asc
	from scipy.optimize import minimize
	
	# Keyword Arguments
	shon=kwargs.get('shon',False)
	fitlim=kwargs.get('fitlim','full')
	scalematch=kwargs.get('scale_to','sapphire')
	laser=kwargs.get('laser','green')
	skips=kwargs.get('skips',[])
	
	if laser=='blue':
		lwave=480
		l_ind=0
		u_ind=180
	if laser=='green':
		lwave=555
		l_ind=450
		u_ind=575
	
	# Delete directory if it exists, and make a new one
	if os.path.exists(savefile)==True:
		shutil.rmtree(savefile)
	os.makedirs(savefile)
	os.makedirs(savefile+'/subfits')
	
	# Load the data
	wsig,datasig=tos.spec_plot(signum,store=True,savefile=savefile+'/raw')
	wbg,databg=tos.spec_plot(bgnum,frames=[1],store=True)
	wdk,datadk=tos.spec_store(darkrun)
	siglv=gr.reader('/Users/christopherchambers/Py/data/'+str(date)+'/lv/'+str(lvnumsig)+'_laser_temp.txt',header=False,delimeter=';')
	sigld,sigldamp,duds=laserdata(wsig[l_ind:u_ind],datasig[:,l_ind:u_ind]-125,savefile,signum,skips=skips)
	bgld,bgldamp,duds2=laserdata(wbg[l_ind:u_ind],databg[:,l_ind:u_ind]-125,savefile,bgnum)
	
	# Configure skips to index notation and save
	skip_ind=duds
	for i in range(len(skips)):
		skip_ind.append(skips[i]-1)
	
	skframes=np.zeros((max(len(skips),len(duds)),2))
	skframes[:len(duds),0]=np.array(duds)+1
	skframes[:len(skips),1]=np.array(skips)
	save(skframes,['laser','manual'],savefile+'/skiped_frames.txt')	
	
	# Make frame list 
	frames=np.array(range(1,np.shape(datasig)[0]+1))
	frames=np.delete(frames,skip_ind,axis=0)
	siglv=np.delete(siglv,skip_ind,axis=1)
	
	# Get bg with closest wavelength
	bgpairs=[]
	for i in range(len(sigld)):
		diff=100
		for j in range(len(bgld)):
			diff2=np.absolute(sigld[i]-bgld[j])
			if diff2>=diff:
				continue
			else:
				diff=diff2
				ind=j
		bgpairs.append((i,ind))
	
	# Subtract darkframe from signal and background
	darkav=np.mean(datadk)
	datasig=datasig-darkav
	databg=databg-darkav
	
	#Scale the matched bg by region above the sapphire fluorescence
	subs=[]
	scales=[]
	for (i,j) in bgpairs:
		def findscale(a):
			if scalematch=='sapphire':
				return np.sum(np.absolute(datasig[i,1300:]-a*databg[j,1300:]))
			if scalematch=='raman':
				return np.sum(np.absolute(datasig[i,315:355]-a*databg[j,315:355]))
		res=minimize(findscale,[3])
		subs.append(datasig[i,:]-res.x*databg[j,:])
		scales.append(res.x)
	subs=np.array(subs)
	subs=np.delete(subs,skip_ind,axis=0)
	sigld=np.delete(sigld,skip_ind,axis=0)
	
	# Restrict data to fit limits
	if fitlim!='full':
		# Find indices for fit limits
		lim_index=[]
		for k in range(2):
			j=-1
			for i in wsig:
				j=j+1
				if i<fitlim[k]:
					last_wave=i
				else:
					a=abs(fitlim[k]-last_wave)
					b=abs(i-fitlim[k])
					if a>b:
						lim_index.append(j)
						break
					else:
						lim_index.append(j-1)
						break	
		wsig=wsig[lim_index[0]:lim_index[1]]
		subs=subs[:,lim_index[0]:lim_index[1]]
	
	# Fit gaussian/lorentzian to the peaks with spec_fit
	def solo_gauss(x,b,a,c):
		return a*np.exp(-(4*np.log(2)*np.square(x-b))/(np.square(c)))
	def solo_lorentz(x,b,a,c):
		return a*(1./(1+np.square((x-b)/(c/2))))
	def solo_asym(x,b,a,c,d):
		return a*(1.-sp.erf((x-b)/c))*(1.+sp.erf((x-b)/d))
	param_all,fits,types=spec_fit_fast(wsig,subs,param_0,'sub')

	# Plot the fits
	data=subs
	n=np.shape(param_all)[1]
	for frame in range(1,np.shape(data)[0]+1):
		f=frame-1
		fit_label=repr(n)+' Component Fit'
		plt.figure('run'+str(signum)+' Fit')
		plt.clf()
		plt.plot(wsig,data[f,:],'b',label='Data')
		plt.plot(wsig,fits[f,:],'k--',linewidth=3,label=fit_label)
		# Plot the individual Gaussians
		colors=['green','red','purple','goldenrod','aquamarine','saddlebrown','magenta','orange','blue','black']
		fit_text,fit_text2,fit_text3='','',''
		for j in range(n):
			if n-j<=1:
				ls=':'
			else:
				ls='-'
			if types[j]=='g':
				type='Gaussian'
				plt.plot(wsig,solo_gauss(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],label=type+repr(j+1))
			if types[j]=='l':
				type='Lorentzian'
				plt.plot(wsig,solo_lorentz(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],label=type+repr(j+1))
			if types[j]=='a':
				type='Asymmetric'
				plt.plot(wsig,solo_asym(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1],param_all[f,j,3]),color=colors[j],label=type+repr(j+1))
			if j==0:
				if types[j]=='a':
					fit_text += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				if types[j] in ['g','l']:
					fit_text += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2)) 
			if 1<=j<=2:
				if types[j]=='a':
					fit_text += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				if types[j] in ['g','l']:
					fit_text += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
			if j==3:
				if types[j]=='a':
					fit_text2 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				if types[j] in ['g','l']:
					fit_text2 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2)) 
			if 4<=j<=5:
				if types[j]=='a':
					fit_text2 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				if types[j] in ['g','l']:
					fit_text2 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
			if j==6:
				if types[j]=='a':
					fit_text3 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				if types[j] in ['g','l']:
					fit_text3 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2)) 
			if 7<=j:
				if types[j]=='a':
					fit_text3 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				if types[j] in ['g','l']:
					fit_text3 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
		plt.title('run'+str(signum)+': frame '+repr(frames[f]))
		plt.axis([wsig[0],wsig[-1],0,1.5*np.amax(data[f,:])])
		plt.xlabel('Wavelength (nm)')
		plt.ylabel('Counts')
		plt.legend(fontsize=8)
		# Add a text box with parameters from fit
		bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
		plt.text(.02,.97,fit_text,fontsize=7,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
		if n>3:
			plt.text(.2,.97,fit_text2,fontsize=7,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
		if n>6:
			plt.text(.38,.97,fit_text3,fontsize=7,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
		plt.show()
		# Save plot to file
		plt.savefig(savefile+'/subfits/run'+str(signum)+'_'+repr(frames[f])+'_fit.png')
		
	# Save amp parameters for each gaussian
	names=[]
	for j in range(n):
		names.append(str(param_all[0,j,0])+'nm')
	save(param_all[:,:,2],names,savefile+'/run'+str(signum)+'_params.txt')
	
	# Plot peaks with wavelength
	#Find peaks of asymgauses
	if shon==True:
		n=n-2
	peaks_norm=np.zeros((np.shape(data)[0],n))
	names2=[]
	for i in range(np.shape(data)[0]):
		for j in range(n):
			if types[j]=='a':
				peaks_norm[i,j]=np.amax(solo_asym(wsig,param_all[i,j,0],param_all[i,j,2],param_all[i,j,1],param_all[i,j,3]))/siglv[1,i]
			else:
				peaks_norm[i,j]=param_all[i,j,2]/siglv[1,i]
			if i==0:
				if types[j]=='a': 
					names2.append(types[j]+':('+str(param_0[j][0])+','+str(param_0[j][1])+','+str(param_0[j][2])+')')
				else:	
					names2.append(types[j]+':('+str(param_0[j][0])+','+str(param_0[j][1])+')')
	intdata=np.insert(peaks_norm,0,sigld,axis=1)
	save(intdata,['ExWavelength']+names2,savefile+'/run'+str(signum)+'_ex_spec_peaks.txt')
	plt.figure('run'+str(signum)+'ex_spec_peak')
	plt.clf()
	for j in range(n):
		plt.plot(sigld,peaks_norm[:,j],'-o',color=colors[j],label=str(round(param_all[0,j,0],2))+'nm')
	plt.axis([np.amin(sigld)-1,np.amax(sigld)+1,0,1.1*np.amax(peaks_norm)])
	plt.title('Excitation Spectrum for run'+str(signum))
	plt.xlabel('Excitation Wavelength(nm)')
	plt.ylabel('Peak Counts/mW')
	plt.legend(loc=2)
	plt.savefig(savefile+'/run'+str(signum)+'ex_spec_peaks.png')
	
	# Make plot of peak for shon lines
	if shon==True:
		peaks_shon=np.zeros((np.shape(data)[0],2))
		names_s=[]
		for i in range(np.shape(peaks_shon)[0]):
			for j in range(2):
				peaks_shon[i,j]=param_all[i,j-5,2]/siglv[1,i]
				if i==0:
					names_s.append(str(round(param_0[j-5][0],1))+'nm')
		intdata2=np.insert(peaks_shon,0,sigld,axis=1)
		save(intdata2,['ExWavelength']+names_s,savefile+'/run'+str(signum)+'_ex_spec_shon.txt')
		plt.figure('run'+str(signum)+'ex_spec_shon')
		plt.clf()
		for j in range(2):
			plt.plot(sigld,peaks_shon[:,j],color=colors[j],label=names_s[j])
		plt.axis([np.amin(sigld)-1,np.amax(sigld)+1,0,1.1*np.amax(peaks_shon)])
		plt.title('Excitation Spectrum for "Shon Lines" in run'+str(signum))
		plt.xlabel('Excitation Wavelength(nm)')
		plt.ylabel('Peak Counts/mW')
		plt.legend()
		plt.savefig(savefile+'/run'+str(signum)+'ex_spec_shon.png')
	
	print('\a')
	print('\a')
	print('\a')
	
	return wsig,subs,param_all,intdata

#-----------------------------------------------------------------------------------------
def excite_raw(signum,lvnumsig,skips,param_0,bgparam_0,savefile,**kwargs):
	
	import the_one_script as tos
	import gen_reader as gr
	import matplotlib.pyplot as plt
	import numpy as np
	import scipy.special as sp
	import math
	import os
	import asciitable as asc
	from scipy.optimize import minimize
	
	# Keyword Arguments
	shon=kwargs.get('shon',False)
	fitlim=kwargs.get('fitlim','full')
	laser=kwargs.get('laser','green')
	
	if laser=='blue':
		lwave=480
		l_ind=0
		u_ind=180	
	if laser=='green':
		lwave=555
		l_ind=450
		u_ind=575
	if laser=='6g':
		lwave=567
		l_ind=200
		u_ind=440	
	
	# Load the data
	wsig,datasig=tos.spec_plot(signum,store=True,savefile=savefile+'/raw')
	wdk,datadk=tos.spec_store(darkrun)
	dkcts=np.mean(datadk)
	siglv=gr.reader('/Users/christopherchambers/Py/data/'+str(date)+'/lv/'+str(lvnumsig)+'_laser_temp.txt',header=False,delimeter=';')
	sigld,sigldamp,duds=laserdata(wsig[l_ind:u_ind],datasig[:,l_ind:u_ind]-dkcts,lwave,savefile,signum,skips=skips)
	
	print('Laser Wavelengths Done')
	
	# Configure skips to index notation
	skip_ind=duds
	for i in range(len(skips)):
		skip_ind.append(skips[i]-1)
	
	# Make frame list and remove skipped frames 
	frames=np.array(range(1,np.shape(datasig)[0]+1))
	frames=np.delete(frames,skip_ind,axis=0)
	siglv=np.delete(siglv,skip_ind,axis=0)
	sigld=np.delete(sigld,skip_ind,axis=0)
	
	# Remove skipped frames and subtract dark frame
	datasig=np.delete(datasig-dkcts,skip_ind,axis=0)
	
	# Create directory if it does not exist	
	if os.path.exists(savefile)==False:
		os.makedirs(savefile)
	if os.path.exists(savefile+'/rawfits')==False:
		os.makedirs(savefile+'/rawfits')
	if os.path.exists(savefile+'/subfits')==False:
		os.makedirs(savefile+'/subfits')
	
	# Restrict data to fit limits
	if fitlim!='full':
		# Find indices for fit limits
		lim_index=[]
		for k in range(2):
			j=-1
			for i in wsig:
				j=j+1
				if i<fitlim[k]:
					last_wave=i
				else:
					a=abs(fitlim[k]-last_wave)
					b=abs(i-fitlim[k])
					if a>b:
						lim_index.append(j)
						break
					else:
						lim_index.append(j-1)
						break	
		wsig=wsig[lim_index[0]:lim_index[1]]
		datasig=datasig[:,lim_index[0]:lim_index[1]]
	
	# Fit gaussian/lorentzian to the peaks with spec_fit
	def solo_gauss(x,b,a,c):
		return a*np.exp(-(4*np.log(2)*np.square(x-b))/(np.square(c)))
	def solo_lorentz(x,b,a,c):
		return a*(1./(1+np.square((x-b)/(c/2))))
	def solo_asym(x,b,a,c,d):
		return a*(1.-sp.erf((x-b)/c))*(1.+sp.erf((x-b)/d))
	param_all,fits,bgfit,types,m=spec_fit_fast(wsig,datasig,param_0,bgparam_0,'raw')
	subs=datasig-bgfit
	
	print('Fits Done')

	# Plot the fits
	data=datasig
	n=np.shape(param_all)[1]-m
	for frame in range(1,np.shape(data)[0]+1):
		f=frame-1
		fit_label=repr(n)+' Component Fit'
		plt.figure('run'+str(signum)+' Fit')
		plt.clf()
		plt.plot(wsig,data[f,:],'b',label='Data')
		plt.plot(wsig,fits[f,:],'k--',linewidth=3,label=fit_label)
		# Plot the individual Gaussians
		colors=['green','red','purple','goldenrod','aquamarine','saddlebrown','magenta','orange','blue','black']
		fit_text,fit_text2,fit_text3='','',''
		for j in range(n):
			if n-j<=1:
				ls=':'
			else:
				ls='-'
			if types[j]=='g':
				type='Gaussian'
				plt.plot(wsig,solo_gauss(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],label=type+repr(j+1))
			if types[j]=='l':
				type='Lorentzian'
				plt.plot(wsig,solo_lorentz(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],label=type+repr(j+1))
			if types[j]=='a':
				type='Asymmetric'
				plt.plot(wsig,solo_asym(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1],param_all[f,j,3]),color=colors[j],label=type+repr(j+1))
			if j==0:
				if types[j]=='a':
					fit_text += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				if types[j] in ['g','l']:
					fit_text += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2)) 
			if 1<=j<=2:
				if types[j]=='a':
					fit_text += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				if types[j] in ['g','l']:
					fit_text += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
			if j==3:
				if types[j]=='a':
					fit_text2 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				if types[j] in ['g','l']:
					fit_text2 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2)) 
			if 4<=j<=5:
				if types[j]=='a':
					fit_text2 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				if types[j] in ['g','l']:
					fit_text2 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
			if j==6:
				if types[j]=='a':
					fit_text3 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				if types[j] in ['g','l']:
					fit_text3 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2)) 
			if 7<=j:
				if types[j]=='a':
					fit_text3 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				if types[j] in ['g','l']:
					fit_text3 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
		plt.plot(wsig,bgfit[f,:],'k',label='Background')
		plt.title('run'+str(signum)+': frame '+repr(frames[f]))
		plt.axis([wsig[0],wsig[-1],0,1.5*np.amax(data[f,:])])
		plt.xlabel('Wavelength (nm)')
		plt.ylabel('Counts')
		plt.legend(fontsize=8)
		# Add a text box with parameters from fit
		bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
		plt.text(.02,.97,fit_text,fontsize=7,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
		if n>3:
			plt.text(.2,.97,fit_text2,fontsize=7,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
		if n>6:
			plt.text(.38,.97,fit_text3,fontsize=7,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
		plt.show()
		# Save plot to file
		plt.savefig(savefile+'/rawfits/run'+str(signum)+'_'+repr(frames[f])+'_fit_raw.png')
		
		# Plot bg-subtracted fits
		plt.figure('run'+str(signum)+' Fit Sub')
		plt.clf()
		plt.plot(wsig,data[f,:]-bgfit[f,:],'b',label='Data')
		plt.plot(wsig,fits[f,:]-bgfit[f,:],'k--',linewidth=3,label='Fit')
		#bgsubdata=np.zeros((n+2,np.shape(wsig)[0]),dtype='float')
		#bgsubdata[0,:]=wsig
		#bgsubdata[1,:]=fits[f,:]-bgfit[f,:]
		names3=['wavelength','full_fit']
		for j in range(n):
			if types[j]=='g':
				plt.plot(wsig,solo_gauss(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],label=str(round(param_all[f,j,0],1))+'nm')
				#bgsubdata[j+2,:]=solo_gauss(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1])
				#names3.append(str(round(param_all[f,j,0],1))+'nm')
			if types[j]=='l':
				plt.plot(wsig,solo_lorentz(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],label=str(round(param_all[f,j,0],1))+'nm')
				#bgsubdata[j+2,:]=solo_lorentz(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1])
				#names3.append(str(round(param_all[f,j,0],1))+'nm')
			if types[j]=='a':
				plt.plot(wsig,solo_asym(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1],param_all[f,j,3]),color=colors[j],label=str(round(param_all[f,j,0],1))+'nm')
				#bgsubdata[j+2,:]=solo_asym(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1],param_all[f,j,3])
				#names3.append(str(round(param_all[f,j,0],1))+'nm')
		plt.axis([wsig[0],wsig[-1],0,1.5*np.amax(data[f,:]-bgfit[f,:])])
		plt.xlabel('Wavelength (nm)')
		plt.ylabel('Counts')
		plt.legend(fontsize=12)
		plt.savefig(savefile+'/subfits/run'+str(signum)+'_'+repr(frames[f])+'_fit_sub.png')
		#save(bgsubdata,names3,savefile+'/subfits/run'+str(signum)+'_'+repr(frames[f])+'_fit_sub.txt')
		
	# Plot peaks with wavelength
	# remove the ?-able frames from sigld and g_ints and save
	if shon==True:
		n=n-2
	peaks_norm=np.zeros((np.shape(data)[0],n))
	names2=[]
	for i in range(np.shape(data)[0]):
		for j in range(n):
			peaks_norm[i,j]=param_all[i,j,2]/siglv[i,1]
			if i==0:
				if types[j]=='a': 
					names2.append(types[j]+':('+str(param_0[j][0])+','+str(param_0[j][1])+','+str(param_0[j][2])+')')
				else:	
					names2.append(types[j]+':('+str(param_0[j][0])+','+str(param_0[j][1])+')')
	intdata=np.insert(peaks_norm,0,sigld,axis=1)
	tos.save(np.transpose(intdata),['ExWavelength']+names2,savefile+'/run'+str(signum)+'_ex_spec_peaks.txt')
	plt.figure('run'+str(signum)+'ex_spec_peak')
	plt.clf()
	for j in range(n):
		plt.plot(sigld,peaks_norm[:,j],'-o',color=colors[j],label=str(round(param_all[0,j,0],2))+'nm')
	plt.axis([np.amin(sigld)-1,np.amax(sigld)+1,0,1.1*np.amax(peaks_norm)])
	plt.title('Excitation Spectrum for run'+str(signum))
	plt.xlabel('Excitation Wavelength(nm)')
	plt.ylabel('Peak Counts/mW')
	plt.legend().draw_frame(False)
	plt.savefig(savefile+'/run'+str(signum)+'ex_spec_peaks.png')
	
	# Make plot of peak for shon lines
	if shon==True:
		peaks_shon=np.zeros((np.shape(data)[0],2))
		names_s=[]
		for i in range(np.shape(peaks_shon)[0]):
			for j in range(2):
				peaks_shon[i,j]=param_all[i,j-5,2]/siglv[i,1]
				if i==0:
					names_s.append(str(round(param_0[j-5][0],1))+'nm')
		intdata=np.insert(peaks_shon,0,sigld,axis=1)
		save(np.transpose(intdata),['ExWavelength']+names_s,savefile+'/run'+str(signum)+'_ex_spec_shon.txt')
		plt.figure('run'+str(signum)+'ex_spec_shon')
		plt.clf()
		for j in range(2):
			plt.plot(sigld,peaks_shon[:,j],color=colors[j],label=names_s[j])
		plt.axis([np.amin(sigld)-1,np.amax(sigld)+1,0,1.1*np.amax(peaks_shon)])
		plt.title('Excitation Spectrum for "Shon Lines" in run'+str(signum))
		plt.xlabel('Excitation Wavelength(nm)')
		plt.ylabel('Peak Counts/mW')
		plt.legend().draw_frame(False)
		plt.savefig(savefile+'/run'+str(signum)+'ex_spec_shon.png')
	
	return wsig,datasig,subs,param_all,intdata
#-----------------------------------------------------------------------------------------
# Function to plot peak height vs time in raw spectra

def bleach_fit(signum,dt,param_0,fit_type,savefile,**kwargs):
	
	import the_one_script as tos
	import gen_reader as gr
	import matplotlib.pyplot as plt
	import scipy.special as sp
	import numpy as np
	import math
	import os
	import asciitable as asc
	from scipy.optimize import minimize
	
	# Keyword Arguments
	shon=kwargs.get('shon',False)
	fitlim=kwargs.get('fitlim','full')
	exclude=kwargs.get('exclude',[])
	readout=kwargs.get('readout',.09818)
	plot_type=kwargs.get('plot_type','Integrated')
	bgint=kwargs.get('integral',False)
	
	# Load the data
	wsig,datasig=tos.spec_store(signum,savefile=savefile+'/raw')
	
	# Make time array (in seconds)
	dt=dt+readout
	time=np.arange(dt,dt*(np.shape(datasig)[0])+.1,dt)
	
	# Make frame list 
	frames=np.array(range(1,np.shape(datasig)[0]+1))
	
	# Create directory if it does not exist	
	if os.path.exists(savefile)==False:
		os.makedirs(savefile)
	if os.path.exists(savefile+'/rawfits')==False:
		os.makedirs(savefile+'/rawfits')
	if fit_type in ['blue','green']:
		if os.path.exists(savefile+'/subfits')==False:
			os.makedirs(savefile+'/subfits')
	
	# Restrict data to fit limits
	if fitlim!='full':
		# Find indices for fit limits
		lim_index=[]
		for k in range(2):
			j=-1
			for i in wsig:
				j=j+1
				if i<fitlim[k]:
					last_wave=i
				else:
					a=abs(fitlim[k]-last_wave)
					b=abs(i-fitlim[k])
					if a>b:
						lim_index.append(j)
						break
					else:
						lim_index.append(j-1)
						break	
		wsig=wsig[lim_index[0]:lim_index[1]]
		datasig=datasig[:,lim_index[0]:lim_index[1]]
	# Remove skipped frames and subtract dark frame
	datasig=datasig-125
	
	# Fit gaussian/lorentzian to the peaks with spec_fit
	def solo_gauss(x,b,a,c):
		return a*np.exp(-(4*np.log(2)*np.square(x-b))/(np.square(c)))
	def solo_lorentz(x,b,a,c):
		return a*(1./(1+np.square((x-b)/(c/2))))
	def solo_asym(x,b,a,c,d):
		return a*(1.-sp.erf((x-b)/c))*(1.+sp.erf((x-b)/d))
	if fit_type=='sub':
		param_all,fits,types=spec_fit_fast(wsig,datasig,param_0,fit_type)
	if fit_type in ['blue','green']:
		param_all,fits,bgfit,types=spec_fit_fast(wsig,datasig,param_0,fit_type)
		subs=datasig-bgfit
	
	print('Fits Done')

	# Plot the fits
	data=datasig
	n=len(param_0)
	ints=np.zeros((np.shape(data)[0],len(param_0)))
	if fit_type=='blue':
		n=n-3
	if fit_type=='green':
		n=n-2
	for frame in range(1,np.shape(data)[0]+1):
		f=frame-1
		fit_label=repr(n)+' Component Fit'
		plt.figure('run'+str(signum)+' Fit')
		plt.clf()
		plt.plot(wsig,data[f,:],'b',label='Data')
		plt.plot(wsig,fits[f,:],'k--',linewidth=3,dashes=(2,2),label=fit_label)
		# Plot the individual Gaussians
		colors=['green','red','purple','goldenrod','aquamarine','saddlebrown','magenta','orange','blue','black']
		fit_text,fit_text2,fit_text3='','',''
		for j in range(n):
			if types[j]=='g':
				type='Gaussian'
				plt.plot(wsig,solo_gauss(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],label=type+repr(j+1))
				ints[f,j]=np.sum(solo_gauss(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]))
			if types[j]=='l':
				type='Lorentzian'
				plt.plot(wsig,solo_lorentz(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],label=type+repr(j+1))
				ints[f,j]=np.sum(solo_lorentz(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]))
			if types[j]=='a':
				type='Asymmetric'
				plt.plot(wsig,solo_asym(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1],param_all[f,j,3]),color=colors[j],label=type+repr(j+1))
				ints[f,j]=np.sum(solo_asym(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1],param_all[f,j,3]))
			if j==0:
				if types[j]=='a':
					fit_text += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				else:
					fit_text += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2)) 
			if 1<=j<=2:
				if types[j]=='a':
					fit_text += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				else:
					fit_text += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
			if j==3:
				if types[j]=='a':
					fit_text2 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				else:
					fit_text2 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2)) 
			if 4<=j<=5:
				if types[j]=='a':
					fit_text2 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				else:
					fit_text2 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
			if j==6:
				if types[j]=='a':
					fit_text3 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				else:
					fit_text3 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2)) 
			if 7<=j:
				if types[j]=='a':
					fit_text3 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nWidth1: '+repr(round(param_all[f,j,1],2))+'\nWidth2: '+repr(round(param_all[f,j,3],2))
				else:
					fit_text3 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
		if fit_type=='blue':
			plt.plot(wsig,bgfit[f,:],'k',label='Background')
		if fit_type=='green':
			plt.plot(wsig,bgfit[f,:],'k',label='Background')
		plt.title('run'+str(signum)+': frame '+repr(frames[f]))
		plt.axis([wsig[0],wsig[-1],0,1.5*np.amax(data[f,:])])
		plt.xlabel('Wavelength (nm)')
		plt.ylabel('Counts')
		plt.legend(fontsize=8)
		# Add a text box with parameters from fit
		bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
		plt.text(.02,.97,fit_text,fontsize=7,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
		if n>3:
			plt.text(.2,.97,fit_text2,fontsize=7,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
		if n>6:
			plt.text(.38,.97,fit_text3,fontsize=7,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
		plt.show()
		# Save plot to file
		plt.savefig(savefile+'/rawfits/run'+str(signum)+'_'+repr(frames[f])+'_fit_raw.png')
		
		# Plot bg-subtracted fits
		if fit_type in ['blue','green']:
			plt.figure('run'+str(signum)+' Fit Sub')
			plt.clf()
			plt.plot(wsig,data[f,:]-bgfit[f,:],'b',label='Data')
			plt.plot(wsig,fits[f,:]-bgfit[f,:],'k--',linewidth=3,label='Fit')
			bgsubdata=np.zeros((n+2,np.shape(wsig)[0]),dtype='float')
			bgsubdata[0,:]=wsig
			bgsubdata[1,:]=fits[f,:]-bgfit[f,:]
			names3=['wavelength','full_fit']
			for j in range(n):
				if types[j]=='g':
					plt.plot(wsig,solo_gauss(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],label=str(round(param_all[f,j,0],1))+'nm')
					bgsubdata[j+2,:]=solo_gauss(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1])
					names3.append(str(round(param_all[f,j,0],1))+'nm')
				if types[j]=='l':
					plt.plot(wsig,solo_lorentz(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],label=str(round(param_all[f,j,0],1))+'nm')
					bgsubdata[j+2,:]=solo_lorentz(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1])
					names3.append(str(round(param_all[f,j,0],1))+'nm')
				if types[j]=='a':
					plt.plot(wsig,solo_asym(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1],param_all[f,j,3]),color=colors[j],label=str(round(param_all[f,j,0],1))+'nm')
					bgsubdata[j+2,:]=solo_asym(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1],param_all[f,j,3])
					names3.append(str(round(param_all[f,j,0],1))+'nm')
			plt.axis([wsig[0],690,0,1.5*np.amax(data[f,:]-bgfit[f,:])])
			plt.xlabel('Wavelength (nm)')
			plt.ylabel('Counts')
			plt.legend(fontsize=12).draw_frame(False)
			plt.savefig(savefile+'/subfits/run'+str(signum)+'_'+repr(frames[f])+'_fit_sub.png')
			#save(bgsubdata,names3,savefile+'/subfits/run'+str(signum)+'_'+repr(frames[f])+'_fit_sub.txt')
	
	# Scale 1st frame
	f1scale=np.mean(param_all[1:,-1,2])/param_all[0,-1,2]
	param_all[0,:,2]=f1scale*param_all[0,:,2]
	# Plot peak height or total counts vs time
	plt.figure('run'+str(signum)+' Bleach')
	plt.clf()
	for j in range(n):
		if j+1 in exclude:
			continue
		if plot_type=='Peak':
			plt.plot(time,param_all[:,j,2]/np.amax(param_all[:,j,2]),'-o',color=colors[j],label=str(round(param_all[0,j,0],1))+'nm')
		if plot_type=='Integrated':
			plt.plot(time,ints[:,j]/np.amax(ints[:,j]),'-o',color=colors[j],label=str(round(param_all[0,j,0],1))+'nm')
	plt.axis([0,1.1*time[-1],0,1.1])
	plt.xlabel('Time (s)')
	plt.ylabel('Normalized '+plot_type+' Counts')
	plt.title('Bleaching Curves for run'+str(signum))
	plt.legend(loc=4)
	plt.savefig(savefile+'/run'+str(signum)+'_bleach_norm.png')
	# Plot un-normalized data
	plt.figure('run'+str(signum)+' Bleach2')
	plt.clf()
	for j in range(n):
		if j+1 in exclude:
			continue
		if plot_type=='Peak':
			plt.plot(time,param_all[:,j,2],'-o',color=colors[j],label=str(round(param_all[0,j,0],1))+'nm')
			plt.axis([0,1.1*time[-1],1,1.1*np.amax(param_all[:,:,2])])
		if plot_type=='Integrated':
			plt.plot(time,ints[:,j],'-o',color=colors[j],label=str(round(param_all[0,j,0],1))+'nm')
			plt.axis([0,1.1*time[-1],.1*np.amin(ints[-1,:n]),1.1*np.amax(ints)])
	plt.xlabel('Time (s)')
	plt.ylabel(plot_type+' Counts')
	plt.title('Bleaching Curves for run'+str(signum))
	plt.legend(loc=4)
	plt.gca().set_yscale('log')
	plt.savefig(savefile+'/run'+str(signum)+'_bleach_raw.png')
	
	# Save peak data
	savedata=np.insert(param_all[:,:,2],0,time,axis=1)
	names2=['time(s)']
	for j in range(n):
		if types[j]=='a': 
			names2.append(types[j]+':('+str(param_0[j][0])+','+str(param_0[j][1])+','+str(param_0[j][2])+')')
		else:	
			names2.append(types[j]+':('+str(param_0[j][0])+','+str(param_0[j][1])+')')
	save(savedata,names2,savefile+'/bleach_data_peak.txt')
	
	# Save int data
	savedata2=np.insert(ints,0,time,axis=1)
	save(savedata2,names2,savefile+'/bleach_data_int.txt')
	
	# Find integrated counts in some range
	if bgint!=False:
		# Find indices for fit limits
		int_index=[]
		for k in range(2):
			j=-1
			for i in wsig:
				j=j+1
				if i<bgint[k]:
					last_wave=i
				else:
					a=abs(bgint[k]-last_wave)
					b=abs(i-bgint[k])
					if a>b:
						int_index.append(j)
						break
					else:
						int_index.append(j-1)
						break
		sig_int=np.zeros((np.shape(data)[0],))
		bg_int=np.zeros((np.shape(data)[0],))
		for f in range(np.shape(data)[0]):
			sig_int[f]=np.sum(data[f,int_ind[0]:int_ind[1]]-bgfit[f,int_ind[0]:int_ind[1]])
			bg_int[f]=np.sum(bgfit[f,int_ind[0]:int_ind[1]])
		# Plot the integrals
		plt.plot(time,sig_int,'-o',color='blue',label='Signal')
		plt.plot(time,bg_int,'-o',color='red',label='Background')
		plt.axis([0,1.1*time[-1],.1*np.amin(bg_int[-1]),1.1*np.amax(sig_int)])
		plt.xlabel('Time (s)')
		plt.ylabel('Total Counts in Range')
		plt.title('Integral Curves for run'+str(signum)+' Range: '+str(bgint[0])+'nm - '+str(bgint[1])+'nm')
		plt.legend(loc=4)
		plt.gca().set_yscale('log')
		plt.savefig(savefile+'/int_'+str(bgint[0])+'_'+str(bgint[1])+'.png')
	
		# Save the integrals
		savedata3=np.zeros((2,np.shape(data)[0]))
		savedata3[0,:]=sig_int
		savedata3[1,:]=bg_int
		save(savedata3,['Signal Counts','BG Counts'],savefile+'/int_'+str(bgint[0])+'_'+str(bgint[1])+'.txt')
		
	print('\a')
	print('\a')
	print('\a')
	
	return time,param_all,ints

#-----------------------------------------------------------------------------------------

# Find total counts after some amount of time from integral data

def tot_cts(filepath,column):

	# Import Modules
	import numpy as np
	import matplotlib.pyplot as plt
	import math
	import os
	import the_one_script as tos
	import gen_reader as gr
	from decimal import Decimal

	# import data
	header,data=gr.reader('/Users/christopherchambers/Py/plots/'+str(date)+filepath+'/bleach_data_int.txt')
	
	# summed the counts in the data
	sums=np.zeros((2,np.shape(data)[0]+1))
	sums[0,:]=np.insert(data[:,0],0,0)
	for i in range(1,np.shape(sums)[1]):
		sums[1,i]=np.sum(data[:i,column])
	
	# create fit function
	def fit(x,a,b):
		return (a*b)*(1-np.exp(-x/b))
	
	# fit to the function
	params=tos.fitter(sums,fit,[10000,100],color='red')
	total='%.2E' % Decimal(params[0]*params[1]*(1-np.exp(-1)))
	
	# plot the cutoff level
	cutoff=np.zeros((np.shape(sums)[1],))
	cutoff[:]=params[0]*params[1]*(1-np.exp(-1))
	plt.figure('Fitter')
	plt.plot(sums[0,:],cutoff,'k--',label='sum cutoff')
	plt.legend(loc=4)
	plt.xlabel('Time (t)')
	plt.ylabel=('Total Counts Collected')
	plt.title('Summed Counts in Peak')
	tos.sci()
	fit_text='1/e Time: '+str(round(params[1],1))+'\nTotal Counts: '+str(total)
	bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
	plt.text(.05,.95,fit_text,fontsize=10,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
	plt.savefig('/Users/christopherchambers/Py/plots/'+str(date)+filepath+'/total_cts.png')
	
	return (params[1],params[0]*params[1]*(1-np.exp(-1)))
	
#-----------------------------------------------------------------------------------------

# Script to track peak shifts

def peak_shift(signum,lvnum,spec_type,param_0,frames,savefile,**kwargs):

	# Keyword Arguments
	fitlim=kwargs.get('fitlim','full')
	
	# Import Modules
	import numpy as np
	import matplotlib.pyplot as plt
	import math
	from lmfit import minimize,Parameters,Parameter
	import os
	import the_one_script as tos
	import gen_reader as gr

	# Create directory if it does not exist	
	if os.path.exists(savefile)==False:
		os.makedirs(savefile)
	if os.path.exists(savefile+'/rawfits')==False:
		os.makedirs(savefile+'/rawfits')
	if os.path.exists(savefile+'/subfits')==False:
		os.makedirs(savefile+'/subfits')

	# Load the data
	wsig,datasig=tos.spec_plot(signum,store=True,savefile=savefile+'/raw')
	siglv=gr.reader('/Users/christopherchambers/Py/data/'+str(date)+'/lv/'+str(lvnum)+'_laser_temp.txt',header=False,delimeter=';')
	# Only Fit Relevant Frames for Laser Peak
	sigld,sigldamp,duds=laserdata(wsig[:180],datasig[frames[0]-1:frames[1]-1,:180]-125,480,savefile,signum,shift=frames[0])
		
	print('Laser Wavelengths Done')
	
	# Slice data to get frames/region of interest
	if fitlim!='full':
		# Find indices for fit limits
		lim_index=[]
		for k in range(2):
			j=-1
			for i in wsig:
				j=j+1
				if i<fitlim[k]:
					last_wave=i
				else:
					a=abs(fitlim[k]-last_wave)
					b=abs(i-fitlim[k])
					if a>b:
						lim_index.append(j)
						break
					else:
						lim_index.append(j-1)
						break	
		wsig=wsig[lim_index[0]:lim_index[1]]
		datasig=datasig[:,lim_index[0]:lim_index[1]]
	datasig=datasig[frames[0]-1:frames[1]-1,:]
	# Do fit of all parameters
	
	# Define functions for the fit (n gaussian or lorentzians)
	n=np.shape(param_0)[0]
	def solo_gauss(x,b,a,c):
		return a*np.exp(-(4*np.log(2)*np.square(x-b))/(np.square(c)))
	def solo_lorentz(x,b,a,c):
		return a*(1./(1+np.square((x-b)/(c/2))))
	def full_fit(x,b,a,c,t):
		func=0
		for i in range(n):
			if t[i]=='g':
				func += a[i]*np.exp(-(4*np.log(2)*np.square(x-b[i]))/(np.square(c[i])))
			if t[i]=='l':
				func += a[i]*(1./(1+np.square((x-b[i])/(c[i]/2))))
		return func
	
	# Set initial parameters (only let amplitude change)
	param_all=[]
	ave_res=[]
	data=datasig
	types=[]
	for k in param_0:
		types.append(k[2])
	#General fit case
	for i in range(np.shape(data)[0]):
		param=Parameters()
		# Generate parameters
		for j in range(n):
			param.add('g'+repr(j+1)+'_amp',value=100,min=0)
			param.add('g'+repr(j+1)+'_center',value=param_0[j][0],min=param_0[j][0]-3,max=param_0[j][0]+3)
			param.add('g'+repr(j+1)+'_width',value=param_0[j][1],min=0,max=2*param_0[j][1])
			if param_0[j][3]=='f':
				param['g'+repr(j+1)+'_center'].set(vary=False)
				param['g'+repr(j+1)+'_width'].set(vary=False)
		# Create residual
		def residual(params,x,data,t):
			a,b,c=[],[],[]
			for j in range(n):
				a.append(params['g'+repr(j+1)+'_amp'].value)
				b.append(params['g'+repr(j+1)+'_center'].value)
				c.append(params['g'+repr(j+1)+'_width'].value)
			return data[i,:]-full_fit(x,b,a,c,t)
		
		# Minimize the residual
		result=minimize(residual,param,args=(wsig,data,types))
		param_frame=[]
		for j in range(n):
			param_frame.append([param['g'+repr(j+1)+'_center'].value,param['g'+repr(j+1)+'_width'].value,param['g'+repr(j+1)+'_amp'].value])
		param_all.append(param_frame)
	param_all=np.array(param_all)
	
	# Make array of final fit for each frame
	fits=np.zeros(np.shape(data))
	bgfit=np.zeros(np.shape(data))
	if spec_type=='raw':
		for i in range(np.shape(data)[0]):
			fits[i,:]=full_fit(wsig,param_all[i,:,0],param_all[i,:,2],param_all[i,:,1],types)
			for j in [-1,-2,-3]:
				bgfit[i,:]+=solo_gauss(wsig,param_all[i,j,0],param_all[i,j,2],param_all[i,j,1])
	if spec_type=='sub':
		for i in range(np.shape(data)[0]):
			fits[i,:]=full_fit(wsig,param_all[i,:,0],param_all[i,:,2],param_all[i,:,1],types)
	subs=datasig-bgfit
	
	print('Fits Done')
	
	# Plot the fits
	if spec_type=='raw':
		n=n-3
	
	# Plot the fits
	data=datasig
	n=np.shape(param_all)[1]-3
	for f in range(np.shape(data)[0]):
		frame=f+frames[0]
		fit_label=repr(n)+' Component Fit'
		plt.figure('run'+str(signum)+' Fit')
		plt.clf()
		plt.plot(wsig,data[f,:],'b',label='Data')
		plt.plot(wsig,fits[f,:],'k--',linewidth=3,label=fit_label)
		# Plot the individual Gaussians
		colors=['green','red','purple','goldenrod','aquamarine','saddlebrown','magenta','orange','blue','black']
		fit_text,fit_text2,fit_text3='','',''
		for j in range(n):
			if param_0[j][3]=='s':
				ls=':'
			else:
				ls='-'
			if types[j]=='g':
				type='Gaussian'
				plt.plot(wsig,solo_gauss(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],linestyle=ls,label=type+repr(j+1))
			if types[j]=='l':
				type='Lorentzian'
				plt.plot(wsig,solo_lorentz(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],linestyle=ls,label=type+repr(j+1))
			if j==0:
				fit_text += type+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nAmplitude: '+repr(round(param_all[f,j,2],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2)) 
			if 1<=j<=2:
				fit_text += '\n\n'+type+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nAmplitude: '+repr(round(param_all[f,j,2],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
			if j==3:
				fit_text2 += type+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nAmplitude: '+repr(round(param_all[f,j,2],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
			if 4<=j<=5:
				fit_text2 += '\n\n'+type+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nAmplitude: '+repr(round(param_all[f,j,2],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
			if j==6:
				fit_text3 += type+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nAmplitude: '+repr(round(param_all[f,j,2],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
			if 7<=j<=8:
				fit_text3 += '\n\n'+type+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nAmplitude: '+repr(round(param_all[f,j,2],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
		plt.plot(wsig,bgfit[f,:],'k',label='Background')
		plt.title('run'+str(signum)+': frame '+repr(frame))
		plt.axis([wsig[0],wsig[-1],0,1.5*np.amax(data[f,:])])
		plt.xlabel('Wavelength (nm)')
		plt.ylabel('Counts')
		plt.legend(fontsize=8)
		# Add a text box with parameters from fit
		bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
		plt.text(.02,.97,fit_text,fontsize=7,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
		if n>3:
			plt.text(.2,.97,fit_text2,fontsize=7,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
		if n>6:
			plt.text(.38,.97,fit_text3,fontsize=7,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
		plt.show()
		# Save plot to file
		plt.savefig(savefile+'/rawfits/run'+str(signum)+'_'+repr(frame)+'_fit_raw.png')
		
		# Plot bg-subtracted fits
		plt.figure('run'+str(signum)+' Fit Sub')
		plt.clf()
		plt.plot(wsig,data[f,:]-bgfit[f,:],'b',label='Data')
		plt.plot(wsig,fits[f,:]-bgfit[f,:],'k--',linewidth=3,label='Fit')
		bgsubdata=np.zeros((n+2,np.shape(wsig)[0]),dtype='float')
		bgsubdata[0,:]=wsig
		bgsubdata[1,:]=fits[f,:]-bgfit[f,:]
		names3=['wavelength','full_fit']
		for j in range(n):
			if param_0[j][3]=='s':
				ls=':'
			else:
				ls='-'
			if types[j]=='g':
				plt.plot(wsig,solo_gauss(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],linestyle=ls,label=str(round(param_all[f,j,0],1))+'nm')
				bgsubdata[j+2,:]=solo_gauss(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1])
				names3.append(str(round(param_all[f,j,0],1))+'nm')
			if types[j]=='l':
				plt.plot(wsig,solo_lorentz(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],linestyle=ls,label=str(round(param_all[f,j,0],1))+'nm')
				bgsubdata[j+2,:]=solo_lorentz(wsig,param_all[f,j,0],param_all[f,j,2],param_all[f,j,1])
				names3.append(str(round(param_all[f,j,0],1))+'nm')
		plt.axis([wsig[0],wsig[-1],0,1.5*np.amax(data[f,:]-bgfit[f,:])])
		plt.title('run'+str(signum)+': frame '+repr(frame)+' subtracted')
		plt.xlabel('Wavelength (nm)')
		plt.ylabel('Counts')
		plt.legend(fontsize=12)
		plt.savefig(savefile+'/subfits/run'+str(signum)+'_'+repr(frame)+'_fit_sub.png')
		
	# Plot the peak shifts
	for j in range(n):
		if param_0[j][3]=='s':
			# Plot Center Shift
			plt.figure('run'+str(signum)+' center shift')
			plt.clf()
			plt.plot(sigld,param_all[:,j,0],'-o')
			plt.xlabel('Excitation Wavelength (nm)')
			plt.ylabel('Emission Wavelength (nm)')
			plt.title('Shift of Wavelength for Peak at '+str(round(np.mean(param_all[:,j,0]),1))+'nm')
			plt.savefig(savefile+'/run'+str(signum)+'_cshift_'+str(int(np.mean(param_all[:,j,0])))+'.png')
			# Plot Width Shift
			plt.figure('run'+str(signum)+' width shift')
			plt.clf()
			plt.plot(sigld,param_all[:,j,1],'-o')
			plt.xlabel('Excitation Wavelength (nm)')
			plt.ylabel('Emission FWHM (nm)')
			plt.title('Shift of FWHM for Peak at '+str(round(np.mean(param_all[:,j,0]),1))+'nm')
			plt.savefig(savefile+'/run'+str(signum)+'_wshift_'+str(int(np.mean(param_all[:,j,0])))+'.png')
				
	print('\a')
	print('\a')
	print('\a')	
		
	return wsig,datasig,subs,param_all	
#-----------------------------------------------------------------------------------------
def bg_sub(runnum,param_0,savefile):

	import the_one_script as tos
	import matplotlib.pyplot as plt
	import numpy as np
	import math
	import os
	
	# Create directories to save file (if they don't already exist)
	if os.path.exists(savefile)==False:
		os.makedirs(savefile)
		os.makedirs(savefile+'/bgfits')
		os.makedirs(savefile+'/bgsubs')
		
	# Load data and select regions of interest for fit
	wlength,rawdata=tos.spec_plot(runnum,frames=[1],store=True)
	rawdata=rawdata-124.5
	xdata=np.concatenate((wlength[275:300],wlength[1280:]))
	data=np.concatenate((rawdata[:,275:300],rawdata[:,1280:]),axis=1)
	
	# Define fit functions and do the fit
	n=len(param_0)
	def solo_gauss(x,b,a,c):
		return a*np.exp(-(4*np.log(2)*np.square(x-b))/(np.square(c)))
	def solo_lorentz(x,b,a,c):
		return a*(1./(1+np.square((x-b)/(c/2))))
	def full_fit(x,b,a,c,t):
		func=0
		for i in range(n):
			if t[i]=='g':
				func += a[i]*np.exp(-(4*np.log(2)*np.square(x-b[i]))/(np.square(c[i])))
			if t[i]=='l':
				func += a[i]*(1./(1+np.square((x-b[i])/(c[i]/2))))
		return func
	param_all,fits,types=spec_fit_fast(xdata,data,param_0,'sub')	
	
	# Plot the bg fits with signal, and the subtractions
	for f in range(np.shape(rawdata)[0]):
		frame=f+1
		plt.figure('run'+str(runnum)+' bg fits')
		plt.clf()
		plt.plot(wlength[275:],rawdata[f,275:],'b-',label='Data')
		plt.plot(wlength[275:],full_fit(wlength[275:],param_all[f,:,0],param_all[f,:,2],param_all[f,:,1],types),'k--',linewidth=3,label='Fit')
		colors=['green','red','purple','goldenrod','aquamarine','saddlebrown','magenta','orange','blue','black']
		fit_text=''
		for j in range(n):
			if types[j]=='g':
				type='Gaussian'
				plt.plot(wlength[275:],solo_gauss(wlength[275:],param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],label=type+repr(j+1))
			if types[j]=='l':
				type='Lorentzian'
				plt.plot(wlength[275:],solo_lorentz(wlength[275:],param_all[f,j,0],param_all[f,j,2],param_all[f,j,1]),color=colors[j],label=type+repr(j+1))
			if j==0:
				fit_text += type+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nAmplitude: '+repr(round(param_all[f,j,2],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2)) 
			if j>=1:
				fit_text += '\n\n'+type+repr(j+1)+': \nCenter: '+repr(round(param_all[f,j,0],2))+'\nAmplitude: '+repr(round(param_all[f,j,2],2))+'\nFWHM: '+repr(round(param_all[f,j,1],2))
		plt.axis([wlength[275],wlength[-1],0,1.1*np.amax(rawdata[f,275:])])
		plt.xlabel('Wavelength (nm)')
		plt.ylabel('Counts')
		plt.title('Run'+str(runnum)+' Background Fit')
		plt.legend(fontsize=10)
		# Add a text box with parameters from fit
		bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
		plt.text(.02,.97,fit_text,fontsize=10,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
		plt.savefig(savefile+'/bgfits/run'+str(runnum)+'_'+repr(frame)+'_bgfit.png')
		
		# Plot the Subtraction
		plt.figure('run'+str(runnum)+' bg subs')
		plt.clf()
		plt.plot(wlength[275:],rawdata[f,275:]-full_fit(wlength[275:],param_all[f,:,0],param_all[f,:,2],param_all[f,:,1],types),'b-')
		tos.zero()
		plt.axis([wlength[275],wlength[-1],-50,1.1*np.amax(rawdata[f,275:]-full_fit(wlength[275:],param_all[f,:,0],param_all[f,:,2],param_all[f,:,1],types))])
		plt.xlabel('Wavelength (nm)')
		plt.ylabel('Counts')
		plt.title('Run'+str(runnum)+' Background Subtracted')
		plt.savefig(savefile+'/bgsubs/run'+str(runnum)+'_'+repr(frame)+'_bgsub.png')
		
	print('\a')
	print('\a')
	print('\a')
	
	
#-----------------------------------------------------------------------------------------

def laserdata(xdata,data,c_0,savefile,runnum,**kwargs):
	
	import numpy as np
	import matplotlib.pyplot as plt
	import the_one_script as tos
	from lmfit import Parameter,Parameters,minimize
	import os
	
	# Keyword Arguments
	skips=kwargs.get('skips',[])
	shift=kwargs.get('shift',0)
	
	if os.path.exists(savefile+'/laser')==False:
		os.makedirs(savefile+'/laser')
	
	def solo_gauss(x,b,a,c):
		return a*np.exp(-(4*np.log(2)*np.square(x-b))/(np.square(c)))
		
	frames=range(1,np.shape(data)[0]+1)
	ldat=np.zeros((3,np.shape(data)[0]))
	ldat[0,:]=frames
	duds=[]
	for i in range(np.shape(data)[0]):
		param=Parameters()
		param.add('amp',value=100,min=0)
		param.add('FWHM',value=3,min=0,max=5)
		if i==0:
			param.add('center',value=c_0,min=c_0-5,max=c_0+5)
		if i>=1:
			param.add('center',value=ldat[1,i-1],min=ldat[1,i-1]-5,max=ldat[1,i-1]+5)
		def residual(params,x,dat):
			b=params['center'].value
			a=params['amp'].value
			c=params['FWHM'].value
			return dat-solo_gauss(x,b,a,c)
		result=minimize(residual,param,args=(xdata,data[i,:]))
		ldat[1,i]=param['center'].value
		ldat[2,i]=param['amp'].value
		if ldat[2,i]<=5:
			duds.append(i)
		plt.figure('laser')
		plt.clf()
		plt.plot(xdata,data[i,:],'b-',label='Data')
		plt.plot(xdata,solo_gauss(xdata,param['center'].value,param['amp'].value,param['FWHM'].value),'k--',linewidth=2,label='Fit')
		plt.xlabel('Wavelength (nm)')
		plt.ylabel('Counts')
		plt.title('run'+str(runnum)+'.'+str(frames[i])+' Laser Fit')
		plt.axis([xdata[0],xdata[-1],0,1.1*np.amax(data[i,:])])
		plt.savefig(savefile+'/laser/'+str(runnum)+'_'+str(frames[i])+'laserfit.png')
	save(ldat,['Frame','Wavelength','Amplitude'],savefile+'/laser/'+str(runnum)+'_laserdata.txt')	
		
	return ldat[1,:],ldat[2,:],duds

#-----------------------------------------------------------------------------------------

# Function to fit background subtracted signal to a gaussian
def spec_fit_fast(xdata,spec,param_0,bgparam_0,spec_type):
	
	# Import Modules
	import numpy as np
	import scipy.special as sp
	import matplotlib.pyplot as plt
	import math
	from lmfit import minimize,Parameters,Parameter
	import os
	
	# Do the Fit
	
	# Define functions for the fit (n gaussians)
	n=np.shape(param_0)[0]
	if bgparam_0=='none':
		m=0
	else:
		m=np.shape(bgparam_0)[0]
	N=n+m
	def solo_gauss(x,b,a,c):
		return a*np.exp(-(4*np.log(2)*np.square(x-b))/(np.square(c)))
	def solo_lorentz(x,b,a,c):
		return a*(1./(1+np.square((x-b)/(c/2))))
	def solo_asym(x,b,a,c,d):
		return a*(1.-sp.erf((x-b)/c))*(1.+sp.erf((x-b)/d))
	def full_fit(x,b,a,c,d,t):
		func=0
		for i in range(N):
			if t[i]=='g':
				func += a[i]*np.exp(-(4*np.log(2)*np.square(x-b[i]))/(np.square(c[i])))
			if t[i]=='l':
				func += a[i]*(1./(1+np.square((x-b[i])/(c[i]/2))))
			if t[i]=='a':
				func += a[i]*(1-sp.erf((x-b[i])/c[i]))*(1+sp.erf((x-b[i])/d[i]))
		return func
	
	# Set initial parameters (only let amplitude change)
	param_all=[]
	ave_res=[]
	data=spec
	centers,FWHM,extra,types=[],[],[],[]
	for k in param_0:
		centers.append(k[0])
		FWHM.append(k[1])
		extra.append(k[2])
		types.append(k[3])
	if bgparam_0!='none':
		for l in bgparam_0:
			centers.append(l[0])
			FWHM.append(l[1])
			extra.append(l[2])
			types.append(l[3])
	#General fit case
	for i in range(np.shape(spec)[0]):
		param=Parameters()
		# Generate parameters
		for j in range(N):
			param.add('g'+repr(j+1)+'_amp',value=100,min=0)
		# Create residual
		def residual(params,x,data,b,c,d,t):
			a=[]
			for j in range(N):
				a.append(params['g'+repr(j+1)+'_amp'].value)	
			return data[i,:]-full_fit(x,b,a,c,d,t)
		
		# Minimize the residual
		result=minimize(residual,param,args=(xdata,data,centers,FWHM,extra,types))
		param_frame=[]
		for j in range(N):
			param_frame.append([centers[j],FWHM[j],param['g'+repr(j+1)+'_amp'].value,extra[j]])
		param_all.append(param_frame)
	param_all=np.array(param_all)
	
	# Make array of final fit for each frame
	fits=np.zeros(np.shape(data))
	bgfit=np.zeros(np.shape(data))
	if spec_type=='raw':
		for i in range(np.shape(data)[0]):
			fits[i,:]=full_fit(xdata,param_all[i,:,0],param_all[i,:,2],param_all[i,:,1],param_all[i,:,3],types)
			for j in range(1,m+1):
				if types[-j]=='g':
					bgfit[i,:]+=solo_gauss(xdata,param_all[i,-j,0],param_all[i,-j,2],param_all[i,-j,1])
				if types[-j]=='l':
					bgfit[i,:]+=solo_lorentz(xdata,param_all[i,-j,0],param_all[i,-j,2],param_all[i,-j,1])
				if types[-j]=='a':
					bgfit[i,:]+=solo_asym(xdata,param_all[i,-j,0],param_all[i,-j,2],param_all[i,-j,1],param_all[i,-j,3])
	if spec_type=='sub':
		for i in range(np.shape(data)[0]):
			fits[i,:]=full_fit(xdata,param_all[i,:,0],param_all[i,:,2],param_all[i,:,1],param_all[i,:,3],types)
	
	if spec_type=='raw':		
		return param_all,fits,bgfit,types,m
	if spec_type=='sub':		
		return param_all,fits,types
	
#-----------------------------------------------------------------------------------------
# Script to make Blue Excitation Spectrum

def bg_check(signum,bgnum,savefile):

	import the_one_script as tos
	import gen_reader as gr
	import matplotlib.pyplot as plt
	import numpy as np
	import math
	import os
	from scipy.optimize import minimize
	
	# Load the data
	wsig,datasig=tos.spec_plot(signum,frames=[1],store=True)
	wbg,databg=tos.spec_plot(bgnum,frames=[1],store=True)
	#siglv=gr.reader('/Users/christopherchambers/Py/data/20140916/lv/'+str(lvnumsig)+'_laser_temp.txt',header=False,delimeter=';')
	#bglv=gr.reader('/Users/christopherchambers/Py/data/20140916/lv/'+str(lvnumbg)+'_laser_temp.txt',header=False,delimeter=';')
	names,sigld=gr.reader('/Users/christopherchambers/Py/plots/20140916/run'+str(signum)+'laserdata.txt')
	names,bgld=gr.reader('/Users/christopherchambers/Py/plots/20140916/run'+str(bgnum)+'laserdata.txt')
	
	# Subtract dark frame from sig and bg
	datasig=datasig-125
	databg=databg-125
	
	# Get bg with closest wavelength
	bgpairs=[]
	for i in range(np.shape(sigld)[1]):
		diff=100
		for j in range(np.shape(bgld)[1]):
			diff2=np.absolute(sigld[1,i]-bgld[1,j])
			if diff2>=diff:
				continue
			else:
				diff=diff2
				ind=j
		bgpairs.append((i,ind))

	#Scale the matched bg by raman cutoff shelf
	scales=[]
	for (i,j) in bgpairs:
		def findscale(a):
			#return np.sum(np.absolute(datasig[i,1300:]-a*databg[j,1300:]))
			return np.sum(np.absolute(datasig[i,315:355]-a*databg[j,315:355]))
		res=minimize(findscale,[3])
		scales.append(res.x)
	
	# Save pairs and scales
	if os.path.exists(savefile)==False:
		os.makedirs(savefile)
		
	pairdata=np.array(bgpairs)+1
	tos.save(np.transpose(pairdata),['signal','background'],savefile+'/run'+str(signum)+'bgpairs.txt')
	scaledata=np.zeros((2,len(scales)))
	scaledata[0,:]=range(1,len(scales)+1)
	scaledata[1,:]=scales
	tos.save(scaledata,['frame','bgscale'],savefile+'/run'+str(signum)+'bgscale.txt')
	
	# Plot scaled bg's
	plt.figure('bg_test')
	for i in range(len(bgpairs)):
		plt.clf()
		plt.plot(wsig,datasig[bgpairs[i][0],:],color='Blue',label='Signal')
		plt.plot(wbg,scales[i]*databg[bgpairs[i][1],:],color='Red',label='Bg scaled by'+str(scales[i]))
		plt.xlabel('Wavelength(nm)')
		plt.ylabel('Counts')
		plt.title('BG Scale Check')
		plt.axis([500,np.amax(wsig),0,1.1*np.amax(datasig[bgpairs[i][0],200:])])
		plt.savefig(savefile+'/run'+str(signum)+'_f'+str(i+1)+'_bg_check.png')
	
	print('\a')
	print('\a')
	print('\a')
	
#-----------------------------------------------------------------------------------------
# Fit a single spectrum
def fit_test(xdata,spec,frame,fit_type,param_0,bgparam_0,**kwargs):

	# Keyword Arguments
	fitlim=kwargs.get('fitlim','full')
	c_range=kwargs.get('c_range',10)
	fwhm_range=kwargs.get('fwhm_range',10)
	
	# Import Modules
	import numpy as np
	import scipy.special as sp
	import matplotlib.pyplot as plt
	import math
	from lmfit import minimize,Parameters,Parameter
	import os
	
	# Define functions for the fit (n gaussian/lorentzians)
	n=np.shape(param_0)[0]
	if bgparam_0=='none':
		m=0
	else:
		m=np.shape(bgparam_0)[0]
	N=n+m
	def solo_gauss(x,b,a,c):
		return a*np.exp(-(4*np.log(2)*np.square(x-b)/(np.square(c))))
	def solo_lorentz(x,b,a,c):
		return a*(1./(1+np.square((x-b)/(c/2))))
	def solo_asym(x,b,a,c,d):
		return a*(1-sp.erf((x-b)/c))*(1+sp.erf((x-b)/d))
	def gauss_fit(x,b,a,c,d,t):
		func=0
		for i in range(N):
			if t[i]=='g':
				func += a[i]*np.exp(-(4*np.log(2)*np.square(x-b[i]))/(np.square(c[i])))
			if t[i]=='l':
				func += a[i]*(1./(1+np.square((x-b[i])/(c[i]/2))))
			if t[i]=='a':
				func += a[i]*(1-sp.erf((x-b[i])/c[i]))*(1+sp.erf((x-b[i])/d[i]))
		return func
	
	# Restrict data to fit limits
	if fitlim!='full':
		# Find indices for fit limits
		lim_index=[]
		for k in range(2):
			j=-1
			for i in xdata:
				j=j+1
				if i<fitlim[k]:
					last_wave=i
				else:
					a=abs(fitlim[k]-last_wave)
					b=abs(i-fitlim[k])
					if a>b:
						lim_index.append(j)
						break
					else:
						lim_index.append(j-1)
						break	
		xdata=xdata[lim_index[0]:lim_index[1]]
		spec=spec[:,lim_index[0]:lim_index[1]]
	
	# Set initial parameters
	param_all=[]
	ave_res=[]
	types=[]
	data=spec
	if bgparam_0!='none':
		for i in bgparam_0:
			param_0.append(i)
	for k in param_0:
		types.append(k[3])
	#General fit case
	param=Parameters()
	# Generate parameters
	for j in range(N):
		dwidth=fwhm_range
		dcenter=c_range
		param.add('g'+repr(j+1)+'_amp',100,min=0)	
		param.add('g'+repr(j+1)+'_center',value=param_0[j][0],min=param_0[j][0]-dcenter/2,max=param_0[j][0]+dcenter/2)
		param.add('g'+repr(j+1)+'_width1',value=param_0[j][1],min=param_0[j][1]-dwidth/2,max=param_0[j][1]+dwidth/2)
		param.add('g'+repr(j+1)+'_width2',value=param_0[j][2],min=param_0[j][2]-dwidth/2,max=param_0[j][2]+dwidth/2)
		if param_0[j][4]=='f':
			param['g'+repr(j+1)+'_center'].set(vary=False)
			param['g'+repr(j+1)+'_width1'].set(vary=False)
			param['g'+repr(j+1)+'_width2'].set(vary=False)		
	# Create residual
	def residual(params,x,data,t):
		a,b,c,d=[],[],[],[]
		for j in range(N):
			a.append(params['g'+repr(j+1)+'_amp'].value)
			b.append(params['g'+repr(j+1)+'_center'].value)
			c.append(params['g'+repr(j+1)+'_width1'].value)
			d.append(params['g'+repr(j+1)+'_width2'].value)
		return data[frame-1,:]-gauss_fit(x,b,a,c,d,t)
		
	# Minimize the residual
	result=minimize(residual,param,args=(xdata,data,types))
	param_frame=[]
	for j in range(N):
		param_frame.append([param['g'+repr(j+1)+'_center'].value,param['g'+repr(j+1)+'_width1'].value,param['g'+repr(j+1)+'_amp'].value,param['g'+repr(j+1)+'_width2'].value])
	param_all=np.array(param_frame)
	
	
	# Make array of final fit
	fit=gauss_fit(xdata,param_all[:,0],param_all[:,2],param_all[:,1],param_all[:,3],types)
	if fit_type!='bg':
		bgfit=0
		for j in range(1,m+1):
			if types[-j]=='g':
				bgfit+=solo_gauss(xdata,param_all[-j,0],param_all[-j,2],param_all[-j,1])
			if types[-j]=='l':
				bgfit+=solo_lorentz(xdata,param_all[-j,0],param_all[-j,2],param_all[-j,1])
			if types[-j]=='a':
				bgfit+=solo_asym(xdata,param_all[-j,0],param_all[-j,2],param_all[-j,1],param_all[-j,3])
	# Plot the fits
	fig=plt.figure('Test Fit')
	f=frame-1
	fit_label=repr(n)+' Gaussian Fit'
	fig.clf()
	ax=fig.add_subplot(111)
	ax.plot(xdata,data[f,:],'b',label='Data')
	ax.plot(xdata,fit,'k--',linewidth=3,dashes=(2,2),label=fit_label)
	# Plot the individual Gaussians
	colors=['green','red','purple','goldenrod','aquamarine','saddlebrown','magenta','orange','navy','limegreen']
	fit_text,fit_text2='',''
	for j in range(n):
		if types[j]=='g':
			type='Gaussian'
			plt.plot(xdata,solo_gauss(xdata,param_all[j,0],param_all[j,2],param_all[j,1]),color=colors[j],label=type+repr(j+1))
		if types[j]=='l':
			type='Lorentzian'
			plt.plot(xdata,solo_lorentz(xdata,param_all[j,0],param_all[j,2],param_all[j,1]),color=colors[j],label=type+repr(j+1))
		if types[j]=='a':
			type='Asymmetric'
			plt.plot(xdata,solo_asym(xdata,param_all[j,0],param_all[j,2],param_all[j,1],param_all[j,3]),color=colors[j],label=type+repr(j+1))
		if j==0:
			if types[j]=='a':
				fit_text += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[j,0],2))+'\nWidth1: '+repr(round(param_all[j,1],2))+'\nWidth2: '+repr(round(param_all[j,3],2))
			else:
				fit_text += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[j,0],2))+'\nFWHM: '+repr(round(param_all[j,1],2)) 
		if 1<=j<=3:
			if types[j]=='a':
				fit_text += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[j,0],2))+'\nWidth1: '+repr(round(param_all[j,1],2))+'\nWidth2: '+repr(round(param_all[j,3],2))
			else:
				fit_text += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[j,0],2))+'\nFWHM: '+repr(round(param_all[j,1],2))
		if j==4:
			if types[j]=='a':
				fit_text2 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[j,0],2))+'\nWidth1: '+repr(round(param_all[j,1],2))+'\nWidth2: '+repr(round(param_all[j,3],2))
			else:
				fit_text2 += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[j,0],2))+'\nFWHM: '+repr(round(param_all[j,1],2)) 
		if 5<=j:
			if types[j]=='a':
				fit_text2 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[j,0],2))+'\nWidth1: '+repr(round(param_all[j,1],2))+'\nWidth2: '+repr(round(param_all[j,3],2))
			else:
				fit_text2 += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[j,0],2))+'\nFWHM: '+repr(round(param_all[j,1],2))
	if bgparam_0!='none':
		plt.plot(xdata,bgfit,'k',label='Background Fit')
	plt.title('Test')
	plt.axis([np.amin(xdata),np.amax(xdata),0,1.5*np.amax(data[f,:])])
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Counts')
	ax.legend(fontsize=8)
	# Add a text box with parameters from fit
	bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
	ax.text(.05,.95,fit_text,fontsize=8,bbox=bbox_props,va='top',ha='left',transform=ax.transAxes)
	if n>=4:
		ax.text(.20,.95,fit_text2,fontsize=8,bbox=bbox_props,va='top',ha='left',transform=ax.transAxes)
	plt.show()
	
	print('\a')
	print('\a')	
	print('\a')		
	
#-----------------------------------------------------------------------------------------
# Fit a single spectrum
def fit_test_ends(xdatas,spec,frame,fit_type,param_0):

	# Import Modules
	import numpy as np
	import matplotlib.pyplot as plt
	import math
	from lmfit import minimize,Parameters,Parameter
	import os

	# Define functions for the fit (n gaussian/lorentzians)
	n=np.shape(param_0)[0]
	def solo_gauss(x,b,a,c):
		return a*np.exp(-(4*np.log(2)*np.square(x-b)/(np.square(c))))
	def solo_lorentz(x,b,a,c):
		return a*(1./(1+np.square((x-b)/(c/2))))
	def gauss_fit(x,b,a,c,t):
		func=0
		for i in range(n):
			if t[i]=='g':
				func += a[i]*np.exp(-(4*np.log(2)*np.square(x-b[i]))/(np.square(c[i])))
			if t[i]=='l':
				func += a[i]*(1./(1+np.square((x-b[i])/(c[i]/2))))
		return func
	
	# Set initial parameters
	param_all=[]
	ave_res=[]
	xdata=np.concatenate((xdatas[275:300],xdatas[1280:]))
	data=np.concatenate((spec[:,275:300],spec[:,1280:]),axis=1)
	types=[]
	for k in param_0:
		types.append(k[2])
	#General fit case
	param=Parameters()
	# Generate parameters
	for j in range(n):
		dwidth=20
		dcenter=20
		param.add('g'+repr(j+1)+'_amp',100,min=0)	
		param.add('g'+repr(j+1)+'_center',value=param_0[j][0],min=param_0[j][0]-dcenter/2,max=param_0[j][0]+dcenter/2)
		param.add('g'+repr(j+1)+'_width',value=param_0[j][1],min=param_0[j][1]-dwidth/2,max=param_0[j][1]+dwidth/2)
		if param_0[j][3]=='f':
			param['g'+repr(j+1)+'_center'].set(vary=False)
			param['g'+repr(j+1)+'_width'].set(vary=False)		
	# Create residual
	def residual(params,x,data,t):
		a,b,c=[],[],[]
		for j in range(n):
			a.append(params['g'+repr(j+1)+'_amp'].value)
			b.append(params['g'+repr(j+1)+'_center'].value)
			c.append(params['g'+repr(j+1)+'_width'].value)
		return data[frame-1,:]-gauss_fit(x,b,a,c,t)
		
	# Minimize the residual
	result=minimize(residual,param,args=(xdata,data,types))
	param_frame=[]
	for j in range(n):
		param_frame.append([param['g'+repr(j+1)+'_center'].value,param['g'+repr(j+1)+'_width'].value,param['g'+repr(j+1)+'_amp'].value])
	param_all=np.array(param_frame)
	
	# Make array of final fit
	fit=gauss_fit(xdatas,param_all[:,0],param_all[:,2],param_all[:,1],types)
	if fit_type=='raw':
		bgfit=0
		for j in [-1,-2,-3]:
			bgfit+=solo_gauss(xdatas,param_all[j,0],param_all[j,2],param_all[j,1])
	# Plot the fits
	fig=plt.figure('Test Fit')
	f=frame-1
	fit_label=repr(n)+' Gaussian Fit'
	fig.clf()
	ax=fig.add_subplot(111)
	ax.plot(xdatas,spec[f,:],'b',label='Data')
	ax.plot(xdatas,fit,'k--',linewidth=3,label=fit_label)
	# Plot the individual Gaussians
	colors=['green','red','purple','goldenrod','aquamarine','saddlebrown','magenta','orange']
	fit_text=''
	if fit_type=='raw':
		m=n-2
	else:
		m=n	
	for j in range(m):
		if types[j]=='g':
			type='Gaussian'
			plt.plot(xdatas,solo_gauss(xdatas,param_all[j,0],param_all[j,2],param_all[j,1]),color=colors[j],label=type+repr(j+1))
		if types[j]=='l':
			type='Lorentzian'
			plt.plot(xdatas,solo_lorentz(xdatas,param_all[j,0],param_all[j,2],param_all[j,1]),color=colors[j],label=type+repr(j+1))
		if j==0:
			if types[j]=='e':
				fit_text += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[j,0],2))+'\nWidth1: '+repr(round(param_all[j,1],2))+'\nWidth2: '+repr(round(param_all[j,3],2))
			else:
				fit_text += type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[j,0],2))+'\nFWHM: '+repr(round(param_all[j,1],2)) 
		if j>=1:
			if types[j]=='e':
				fit_text += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[j,0],2))+'\nWidth1: '+repr(round(param_all[j,1],2))+'\nWidth2: '+repr(round(param_all[j,3],2))
			else:
				fit_text += '\n\n'+type+' '+repr(j+1)+': \nCenter: '+repr(round(param_all[j,0],2))+'\nFWHM: '+repr(round(param_all[j,1],2))
	if fit_type=='raw':
		plt.plot(xdata,bgfit,'k',linewidth=2,label='Background Fit')
	plt.title('Test')
	plt.axis([np.amin(xdatas),np.amax(xdatas),0,1.5*np.amax(data[f,:])])
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Counts')
	ax.legend(fontsize=8)
	# Add a text box with parameters from fit
	bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
	ax.text(.05,.97,fit_text,fontsize=8,bbox=bbox_props,va='top',ha='left',transform=ax.transAxes)
	plt.show()
	
	
	
	print('\a')
	print('\a')	
	print('\a')			
#-----------------------------------------------------------------------------------------
# Script to save to ascii
def save(data,names,outfile):

	# Import modules
	import numpy as np
	import asciitable
	import os
	
	# Make a file to write the data to
	open(outfile,'a')
	
	# Find the datatype for each set of data
	data_type=[]
	type_str=['float','int']
	for i in data:
		k=0
		for j in [float,int]:
			if isinstance(i[0],j)==True:
				data_type.append(type_str[k])
				break
			k += 1

	# Initialize array
	data_array=np.zeros(len(data[0]),dtype=zip(names,data_type))
	data_array[:]=map(None,*data)
	
	# Delete outfile if it exists
	for i in range(1,len(outfile)+1):
		if outfile[-i]=='/':
			outdir=outfile[:-i+1]
			outname=outfile[-i+1:]
			break
		else:
			continue
	
	for file in os.listdir(outdir):
		if file==outname:
			os.listdir(outdir).remove(outname)
			break
		else:
			continue
		
	# Write data with asciitable
	asciitable.write(data_array,outfile)
