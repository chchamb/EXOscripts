# Scripts to read, process and plot scope data
date=20170630

#-----------------------------------------------------------------------------------------
# Script to read scope channels

def sc_read(run,**kwargs):
	
	folder=kwargs.get('folder',1)
	chan=kwargs.get('chan','all')
	
	import numpy as np
	import glob
	
	if chan=='comp':
		chan=[3,4]
	
	if chan=='comp+':
		chan=[1,3,4]
		
	if chan=='all':
		chan=np.arange(1,len(glob.glob('/home/chris/anaconda/data/'+str(date)+'/scope/ALL'+str(run).zfill(4)+'/*.CSV'))+1)
		if len(chan)==0:
			print('run not in file')
		
	# Open the file
	volts=[]
	for i in chan:	
		if folder==1:	
			datafile=open('/home/chris/anaconda/data/'+str(date)+'/scope/ALL'+str(run).zfill(4)+'/F'+str(run).zfill(4)+'CH'+str(i)+'.CSV','r')
		if folder==2:	
			datafile=open('/home/chris/anaconda/data/'+str(date)+'/scope_2/ALL'+str(run).zfill(4)+'/F'+str(run).zfill(4)+'CH'+str(i)+'.CSV','r')

		# Initalize data arrays
		time,volt=[],[]
		for line in datafile:
				line=line.strip()
				values=line.split(',')
				if values==[]:
					break
				time.append(float(values[3]))
				volt.append(float(values[4]))
		# Close the datafile, what were you, raised in MATLAB?
		datafile.close()
	
		time=np.array(time)
		volts.append(volt)
	time=time-time[0]
	volts=np.array(volts)

	return time,volts
	
#-----------------------------------------------------------------------------------------
# Script to smooth data (moving average)
def smoother(xdata,ydata,bin,**kwargs):

	savefile=kwargs.get('savefile',False)
	
	import numpy as np
	import matplotlib.pyplot as plt

	# Do a moving average
	ysmooth=np.zeros(np.shape(ydata))
	for i in range(np.shape(ydata)[0]):
		if i<bin:
			continue
		if i>=np.shape(ydata)[0]-(bin+1):
			continue
		else:
			ysmooth[i]=np.mean(ydata[i-bin:i+bin+1])
	# Slice off unsmoothed parts
	xdata2=xdata[bin:-(bin+1)]
	ysmooth=ysmooth[bin:-(bin+1)]
		
	# Plot original and smoothed
	plt.figure('smoother')
	plt.clf()
	plt.plot(xdata,ydata,'b-',label='Raw')
	plt.plot(xdata2,ysmooth,'k--',dashes=(2,2),label='Smoothed: bin size= '+str(2*bin))
	plt.legend(fontsize=12)
	plt.xlim(xdata[0],xdata[-1])
	plt.title('Smoothed Data')
	if savefile!=False:
		plt.savefig(savefile)
	
	return xdata2,ysmooth

#-----------------------------------------------------------------------------------------
# Numerical differentiator

def num_diff(xdata,ydata,bin,**kwargs):
	
	plot=kwargs.get('plot',True)

	import numpy as np
	import matplotlib.pyplot as plt

	# Do a moving average
	dy=np.zeros(np.shape(ydata))
	for i in range(np.shape(ydata)[0]):
		if i<bin:
			dy[i]=(ydata[i+bin]-ydata[0])/(xdata[i+bin]-xdata[0])
		if i>=np.shape(ydata)[0]-(bin+1):
			dy[i]=(ydata[-1]-ydata[i-bin])/(xdata[-1]-xdata[i-bin])
		else:
			dy[i]=(ydata[i+bin+1]-ydata[i-bin])/(xdata[i+bin+1]-xdata[i-bin])
			
	# Plot function and derivative in separate windows
	if plot==True:
		plt.figure('function')
		plt.clf()
		plt.plot(xdata,ydata,'b-',label='Function')
		plt.legend(fontsize=12)
		plt.xlim(xdata[0],xdata[-1])
		plt.figure('derivative')
		plt.clf()
		plt.plot(xdata,dy,'b-',label='Derivative')
		plt.legend(fontsize=12)
		plt.xlim(xdata[0],xdata[-1])
	
	return dy	
		
#-----------------------------------------------------------------------------------------
# Numerical differentiator using dt

def num_diff_dt(dt,ydata,bin,**kwargs):
	
	plot=kwargs.get('plot',True)

	import numpy as np
	import matplotlib.pyplot as plt

	# Do a moving average
	dy=np.zeros(np.shape(ydata))
	for i in range(np.shape(ydata)[0]):
		if i<bin:
			dy[i]=(ydata[i+bin]-ydata[0])/((i+bin)*dt)
		if i>=np.shape(ydata)[0]-(bin+1):
			dy[i]=(ydata[-1]-ydata[i-bin])/((np.shape(ydata)[0]-1-i+bin)*dt)
		else:
			dy[i]=(ydata[i+bin+1]-ydata[i-bin])/(2*bin*dt)
			
	# Plot function and derivative in separate windows
	if plot==True:
		plt.figure('function')
		plt.clf()
		plt.plot(xdata,ydata,'b-',label='Function')
		plt.legend(fontsize=12)
		plt.xlim(xdata[0],xdata[-1])
		plt.figure('derivative')
		plt.clf()
		plt.plot(xdata,dy,'b-',label='Derivative')
		plt.legend(fontsize=12)
		plt.xlim(xdata[0],xdata[-1])
	
	return dy	
		
#-----------------------------------------------------------------------------------------
# Script to transform preamp signal to real signal

def preamp_recon(run,**kwargs):
		
	s_bin=kwargs.get('s_bin',2)
	d_bin=kwargs.get('d_bin',1)
	int_w=kwargs.get('int_w',2)
	RM=kwargs.get('RM',9.61E6)       #.2mv/fC: 8.9E6 
	tau=kwargs.get('tau',10.3E-6)	 #.2mV/fC: 43.6E-6
	savepath=kwargs.get('savepath',False)
	bndfind=kwargs.get('bounds','auto')
	folder=kwargs.get('folder',1)
	inout=kwargs.get('in_out','in')
	cutoff=kwargs.get('cutoff',7.)
	chan=kwargs.get('chan',3)
	chantype=kwargs.get('chantype','comp')

	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from lmfit import minimize,Parameters,Parameter
	import scipy.optimize as opt
	from decimal import Decimal
	import math
	import os
	import utility as ut
	
	
	# Get data from scope
	time,v_preamp_raw=sc_read(run[0],folder=folder,chan=chantype)
	time,v_null_raw=sc_read(run[1],folder=folder,chan='comp')
	if chantype=='all':
		v_preamp_raw=v_preamp_raw[chan-1,:]
		v_null_raw=v_null_raw[chan-1,:]
	if chantype=='comp':
		v_preamp_raw=v_preamp_raw[chan-3,:]
		v_null_raw=v_null_raw[chan-3,:]
		
	# Calculate dt
	dt=ut.avediff(time)
	
	# Smooth the data
	v_null=ut.smoother(v_null_raw,s_bin)
	if savepath==False:
		v_preamp=ut.smoother(v_preamp_raw,s_bin)
	if savepath!=False:
		v_preamp=ut.smoother(v_preamp_raw,s_bin)#,savefile=savefile+'run'+str(run[0])+'/run'+str(run[0])+'_smoothed.png')
	
	# Average all points before trigger
	i=0
	v_neg,v_neg0=[],[]
	while time[i]<(0-80*dt):
		v_neg.append(v_preamp[i])
		v_neg0.append(v_null[i])
		i+=1	
	
	# Zero the signal before the trigger
	v_preamp=v_preamp-np.mean(v_neg)
	v_null=v_null-np.mean(v_neg0)
	
	# Subtract null (no ions, cup 1 in) run
	v_preamp=v_preamp-v_null
	
	# Do the reconstruction using parameters found in calibration
	i_recon=-(1/RM)*(v_preamp+tau*num_diff(time,v_preamp,d_bin,plot=False))
	
	# Slice noise off of the front
	t_2us=ut.bound_finder(time,[cutoff*1E-6])
	time=time[t_2us[0]:]
	i_recon=i_recon[t_2us[0]:]
	if inout=='out':
		i_recon=-1*i_recon
	
	# Integrate peak to get charge
	# Fit peak to find bounds or use manual ones
	if bndfind=='auto':
		peak_t=time[np.where(i_recon==np.amax(i_recon))[0][0]]
		def gauss(x,a,b,c):
			return a*np.exp(-(np.square(x-b))/(2*np.square(c)))
		param,cov=opt.curve_fit(gauss,time,i_recon,[np.amax(i_recon),peak_t,1E-6])
		bnds=ut.bound_finder(time,[param[1]-int_w*np.absolute(param[2]),param[1]+int_w*np.absolute(param[2])])
	if bndfind!='auto':
		bnds=ut.bound_finder(time,[bndfind[0],bndfind[1]])
	charge=np.sum(i_recon[bnds[0]:bnds[1]])*dt

	# Convert to convenient units 
	i_pow,i_fix=ut.prefix(np.amax(i_recon))
	t_pow,t_fix=ut.prefix(np.median(time))
	time_c=np.power(10,t_pow)*time
	i_recon_c=np.power(10,i_pow)*i_recon
	if bndfind=='auto':
		amp=param[0]*np.power(10,i_pow)
		center=param[1]*np.power(10,t_pow)
		width=param[2]*np.power(10,t_pow)

	# Save data as txt
	savedata=np.zeros((np.shape(time)[0],2))
	savedata[:,0]=time
	savedata[:,1]=i_recon
	names=['Time(s)','I_recon(A) ['+ut.conv(RM)+'ohm,'+ut.conv(tau)+'s,'+ut.conv(charge)+'C]']

	# Plot the reconstructed current
	plt.figure('reconstruction')
	plt.clf()
	plt.plot(time_c,i_recon_c,'b-')
	plt.axvline(x=time_c[bnds[0]],color='k',ls='--',dashes=(2,2))
	plt.axvline(x=time_c[bnds[1]],color='k',ls='--',dashes=(2,2))
	plt.xlim(time_c[0],time_c[-1])
	plt.xlabel('Time ('+t_fix+'s)')
	plt.ylabel('Current ('+i_fix+'A)')
	plt.title('Reconstructed Signal run'+str(run))
	bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
	fit_text='Tau: '+ut.conv(tau)+'s\nRM: '+ut.conv(RM)+'ohm\nQ: '+ut.conv(charge)+'C'
	int_text='Integral Bounds: ('+ut.conv(time[bnds[0]])+'s , '+ut.conv(time[bnds[1]])+'s)'
	plt.text(.7,.95,fit_text,fontsize=12,bbox=bbox_props,va='top',ha='left',transform=plt.gca().transAxes)
	if bndfind=='auto':
		plt.text(center+(2*int_w*width),.5*amp,int_text,fontsize=12,bbox=bbox_props,va='top',ha='left')
	if savepath!=False:
		plt.savefig(savepath+'run'+str(run[0])+'_'+str(chan)+'_recon.png')
	
	return time,i_recon,charge	


#-----------------------------------------------------------------------------------------
# Script to fit parameters for calibration

def preamp_cal(run,real_ch,preamp_ch,param_0,R,savefile,**kwargs):

	s_bin=kwargs.get('s_bin',12)
	d_bin=kwargs.get('d_bin',1)

	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from lmfit import minimize,Parameters,Parameter
	import scipy.optimize as opt
	from decimal import Decimal
	import math
	import the_one_script as tos
	import os
	
	# Make directory for data
	if os.path.exists(savefile)==False:
		os.makedirs(savefile)
	
	# Get data from scope
	time_raw,v_real_raw=sc_read(run,real_ch)
	time_raw,v_preamp_raw=sc_read(run,preamp_ch)
	dt=time_raw[1]-time_raw[0]
	
	# Smooth the data
	time,v_real=smoother(time_raw,v_real_raw,s_bin,savefile=savefile+'/run'+str(run)+'_real_smoothed.png')
	time,v_preamp=smoother(time_raw,v_preamp_raw,s_bin,savefile=savefile+'/run'+str(run)+'_preamp_smoothed.png')
	
	# Average all points before trigger
	i=0
	v_p_neg,v_r_neg=[],[]
	while time[i]<(0-2*dt):
		v_p_neg.append(v_preamp[i])
		v_r_neg.append(v_real[i])
		i+=1	
	
	# Zero the signal before the trigger
	v_real=v_real-np.mean(v_r_neg)
	v_preamp=v_preamp-np.mean(v_p_neg)
	
	# Define function for preamp sig
	def preamp_trans(tv,a,b,R):
		return -(1/a)*(tv[1]+b*num_diff(tv[0],tv[1],d_bin,plot=False))
	# Define function for decay of signal
	def decay(x,a,b,c):
		return a*np.exp(-(x-b)/c)
		
	# Find the pulse signal
	edge=[]
	edge.append(np.where(v_real==np.amax(v_real))[0][0])
	edge.append(np.where(v_real==np.amin(v_real))[0][0])
	width=edge[1]-edge[0]
	
	# Slice data to just get decay curve
	t_dec=time[edge[1]+10:]
	v_dec=v_preamp[edge[1]+10:]
	
	# Find guess for tau
	tau_ind=0
	while v_dec[tau_ind]<v_dec[0]/math.e:
		tau_ind+=1	
	
	# Fit the decay curve 
	p_dec,cov_dec=opt.curve_fit(decay,t_dec,v_dec,[v_dec[0],t_dec[0],t_dec[tau_ind]])
	tau_pow,tau_fix=tos.prefix(p_dec[2])
	tau_c=round(p_dec[2]*np.power(10,tau_pow),1)
	
	# Plot decay fit
	plt.figure('decay')
	plt.clf()
	plt.plot(t_dec,v_dec,'b-',label='Preamp Decay')
	plt.plot(t_dec,decay(t_dec,p_dec[0],p_dec[1],p_dec[2]),'k--',label='Fit: tau= '+str(tau_c)+tau_fix+'s')
	plt.legend(fontsize=12,loc=2)
	plt.title('Time Constant Fit')
	plt.savefig(savefile+'/run'+str(run)+'_preamp_cal_tau.png')
	
	plt.xlim(t_dec[0],t_dec[-1])
	plt.legend(loc=2,fontsize=12)
	
	#Convert v_real to current
	i_real=v_real/R	
	# Fit the function to the input data
	preamp_sig=[time,v_preamp]
	# Find the pulse signal
	edge=[]
	edge.append(np.where(v_real==np.amax(v_real))[0][0])
	edge.append(np.where(v_real==np.amin(v_real))[0][0])
	width=edge[1]-edge[0]
	inds=[edge[0]+width/3,edge[1]-width/3]
	# Slice the data to get middle third of pulse
	time_fit=time[inds[0]:inds[1]]
	i_real_fit=i_real[inds[0]:inds[1]]
	preamp_sig_fit=[time_fit,v_preamp[inds[0]:inds[1]]]
	
	# Set initial parameters
	param=Parameters()
	param.add('RM',value=param_0)
	
	# Create residual
	def residual(params,x,data,tau,R):
		a=params['RM'].value	
		return data-preamp_trans(x,a,tau,R)
	
	# Do the fit
	result=minimize(residual,param,args=(preamp_sig_fit,i_real_fit,p_dec[2],R))
	
	# Find convenient units and convert
	t_pow,t_fix=tos.prefix(np.median(time_fit))
	i_pow,i_fix=tos.prefix(np.amax(i_real_fit))
	time_c=np.power(10,t_pow)*time
	i_real_c=np.power(10,i_pow)*i_real
	i_recon_c=np.power(10,i_pow)*preamp_trans(preamp_sig,param['RM'].value,p_dec[2],R)
	
	# Save data as txt
	savedata=np.zeros((np.shape(time)[0],3))
	savedata[:,0]=time
	savedata[:,1]=i_real
	savedata[:,2]=preamp_trans(preamp_sig,param['RM'].value,p_dec[2],R)
	names=['Time(s)','I_real(A)','I_recon(A)']
	tos.save(np.transpose(savedata),names,savefile+'/run'+str(run)+'_preamp_cal_data.txt')
	
	# Plot results 
	plt.figure('calibration')
	plt.clf()
	plt.plot(time_c,i_real_c,'b-',label='Input Signal')
	plt.plot(time_c,i_recon_c,'k--',dashes=(2,2),label='Reconstructed Signal')
	plt.legend(fontsize=12)
	plt.xlabel('Time('+t_fix+'s)')
	plt.ylabel('Current('+i_fix+'A)')
	plt.title('Reconstucted Signal')
	plt.xlim(time_c[0],time_c[-1])
	bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
	fit_text='Tau: '+str(tau_c)+' '+tau_fix+'s\nRM: '+tos.conv(param['RM'].value)+'ohm'
	plt.text(.7,.05,fit_text,fontsize=12,bbox=bbox_props,va='bottom',ha='left',transform=plt.gca().transAxes)
	plt.savefig(savefile+'/run'+str(run)+'_preamp_cal.png')

	return 	preamp_trans(preamp_sig,param['RM'].value,p_dec[2],R),param['RM'].value,p_dec[2]
	
#-----------------------------------------------------------------------------------------
# Fourier Transform

def fourier(run,chan,**kwargs):

	# Define keywords
	xlim=kwargs.get('xlim','default')
	smooth=kwargs.get('smooth',False)

	import numpy as np
	import scipy.fftpack as fourier
	import matplotlib.pyplot as plt
	import gen_reader as gr
	import utility as ut
	
	# Read the data
	t,func=sc_read(run,chan)
	
	# Smooth the data
	if smooth!=False:
		t,func=smoother(t,func,smooth)
	
	# Determine the spacing of the time domain
	t=t-t[0]
	dt,j=0.,0
	for i in range(1,np.shape(t)[0]):
		dt+=t[i]-t[i-1]
		j+=1
	dt=dt/j
	
	# Construct frequency domain (in hertz)
	frq=np.linspace(0,1./(2*dt),np.shape(t)[0])
	
	# Do the Fourier transform
	fftfunc=np.abs(fourier.rfft(func))
	
	# Plot the Fourier Transform and cutoff
	plt.figure('fourier_'+str(run))
	plt.clf()
	plt.plot(frq,fftfunc,'b-',label='Fourier Transform')
	plt.xlabel('Frequency (Hz)')
	plt.ylabel('Fourier Transform')
	if xlim!='default':
		plt.xlim(xlim)
	plt.ylim(1.1*np.amin(fftfunc[1:]),1.1*np.amax(fftfunc[1:]))
	ut.logx()
	plt.title('Fourier Transform for Run '+str(run))
	
	return frq,fftfunc
	
#-----------------------------------------------------------------------------------------
# Lowpass Filter

def lowpass(t,func,cutoff):

	import numpy as np
	import scipy.fftpack as fourier
	import matplotlib.pyplot as plt
	import the_one_script as tos 
	
	# Determine the spacing of the time domain
	t=t-t[0]
	dt,j=0.,0
	for i in range(1,np.shape(t)[0]):
		dt+=t[i]-t[i-1]
		j+=1
	dt=dt/j
	
	# Construct frequency domain (in hertz)
	frq=np.linspace(0,1./(2*dt),np.shape(t)[0])
	
	# Do the Fourier transform
	fftfunc=fourier.rfft(func)
	
	# Plot the Fourier Transform and cutoff
	plt.figure('fourier')
	plt.clf()
	plt.plot(frq,fftfunc,'b-',label='Fourier Transform')
	plt.axvline(x=cutoff,color='k',ls='--',dashes=(2,2),label='Cutoff')
	plt.xlabel('Frequency (Hz)')
	plt.ylabel('Fourier Transform')
	plt.legend()
	tos.logx()
	
	# Apply the cut
	for i in range(np.shape(fftfunc)[0]):
		if frq[i]>cutoff:
			fftfunc[i]=0
	
	# Do the inverse Fourier transform
	func_cut=fourier.irfft(fftfunc)
	
	# Plot original and cut function
	plt.figure('fourier2')
	plt.clf()
	plt.plot(t,func,'b-',label='Original Function')
	plt.plot(t,func_cut,'r-',label='Filtered Function')
	plt.xlabel('Time (s)')
	plt.ylabel('Function')
	plt.legend()
	
	return func_cut

#------------------------------------------------------------------------------------------
# Integrate induction plate signal to get fC/pulse

def pulse_charge_lin(runs,null_pulse,**kwargs):
	
	s_bin=kwargs.get('s_bin',20)
	d_bin=kwargs.get('d_bin',1)
	int_w=kwargs.get('int_w',3)
	bndfind=kwargs.get('bounds','auto')
	savepath=kwargs.get('savepath',False)
	folder=kwargs.get('folder',1)
	chan=kwargs.get('chan',3)
	
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	
	if chan==3:
		RM=9.61E6
		tau=10.3E-6
		parity=1
	if chan==4:
		RM=8.9E6
		tau=43.6E-6
		parity=-1
		
	# Use the reconstruction program to get real signals
	charge=np.zeros((2,np.shape(runs)[1]))
	for i in range(np.shape(runs)[1]):
		if runs[1,i]==0:
			charge[0,i]=runs[0,i]
			charge[1,i]=0.
		else:
			t,I,q=preamp_recon([runs[1,i],null_pulse],savepath=savepath,s_bin=s_bin,d_bin=d_bin,int_w=int_w,bounds=bndfind,folder=folder,chan=chan,RM=RM,tau=tau)
			charge[0,i]=runs[0,i]
			charge[1,i]=q*1e15*parity
	if savepath!=False:
		ut.save(charge,['Run','fC_per_pulse'],savepath+'pulse_charge_fC.txt') 	
		
	return charge

#-----------------------------------------------------------------------------------------
# Integrate induction plate signal to get fC/pulse

def pulse_charge(runs,null_pulse,**kwargs):
	
	s_bin=kwargs.get('s_bin',2)
	d_bin=kwargs.get('d_bin',1)
	int_w=kwargs.get('int_w',3)
	bndfind=kwargs.get('bounds','auto')
	savepath=kwargs.get('savepath',False)
	folder=kwargs.get('folder',1)
	inout=kwargs.get('in_out','in')
	chan=kwargs.get('source','cup')
	chan2=kwargs.get('chan','comp')
	cutoff=kwargs.get('cutoff',10.)
	
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	
	if chan=='cup':
		RM=9.61E6
		tau=10.3E-6
		ch=3
	if chan=='plates':
		RM=8.9E6
		tau=43.6E-6
		ch=4
		
	# Use the reconstruction program to get real signals
	charge=np.zeros((2,len(runs)))
	for i in range(len(runs)):
		t,I,q=preamp_recon([runs[i],null_pulse],savepath=savepath,cutoff=cutoff,s_bin=s_bin,d_bin=d_bin,int_w=int_w,bounds=bndfind,folder=folder,RM=RM,tau=tau,in_out=inout,chantype=chan2,chan=ch)
		charge[0,i]=runs[i]
		charge[1,i]=q*1e15
	if savepath!=False:
		ut.save(charge,['Run','fC_per_pulse'],savepath+'pulse_charge_fC.txt') 	
		
	return charge[1,:]

#-----------------------------------------------------------------------------------------
# Script to transform preamp signal to real signal

def sc_plot(run,**kwargs):
		
	s_bin=kwargs.get('s_bin',False)
	bgrun=kwargs.get('bgrun',False)
	savepath=kwargs.get('savepath',False)
	overlay=kwargs.get('overlay',False)
	chan=kwargs.get('chan','all')
	store=kwargs.get('store',False)
	
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
			
	# Get data from scope
	time,volt=sc_read(run,chan=chan)
	time=time-time[0]
	if bgrun!=False:
		time_bg,volt_bg=sc_read(bgrun,chan=chan)
		volt=volt-volt_bg
			
	# Smooth, if desired
	if s_bin!=False:
		volt2=np.zeros(np.array(np.shape(volt))-np.array([0,2*s_bin+1]))
		for i in range(len(chan)):
			time,volt2[i,:]=smoother(time,volt[i,:],s_bin)
	
	# Plot it up
	plt.figure('scope_plot')
	if overlay==False:
		plt.clf()
	for i in range(len(chan)):
		plt.plot(time,volt[i-1,:])
	plt.axis([time[0],time[-1],1.1*np.amin(volt),1.1*np.amax(volt)])
	if savepath!=False:
		# Make directory for data
		ut.create_dir(savepath+'scope/run'+str(run)+'/')
		plt.savefig(savepath+'scope/run'+str(run)+'/run_'+str(run)+'_scope_plot.png')
	
	if store==True:
		return time,volt

#------------------------------------------------------------------------------------------
# Script to test threshold trigger

def thresh_trig(xdata,ydata,thresh,**kwargs):
	
	delay=kwargs.get('delay',0)
	t0=kwargs.get('t0',0.)
	sign=kwargs.get('sign','pos')
	savepath=kwargs.get('savepath',False)
	lab=kwargs.get('label','')
	
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	
	# Look for triigers at threshhold
	if sign=='pos':
		trigs,j=[0],0
		for i in range(len(ydata)):
			if xdata[i]<=t0:
				continue
			if j!=0 and xdata[i]-trigs[-1]<=delay:
				continue
			if ydata[i]>=thresh:
				trigs.append(xdata[i])
				j+=1
	if sign=='neg':
		trigs,j=[0],0
		for i in range(len(ydata)):
			if xdata[i]<=t0:
				continue
			if j!=0 and xdata[i]-trigs[-1]<=delay:
				continue
			if ydata[i]<=thresh:
				trigs.append(xdata[i])
				j+=1
				
	# Throw out trig at 0
	trigs=np.array(trigs[1:])
	
	# get trigger separation
	dtrigs=np.zeros((len(trigs)-1,))
	for i in range(len(trigs)-1):
		dtrigs[i]=trigs[i+1]-trigs[i]
	
	#Plot signal and triggers
	plt.figure('trig')
	plt.clf()
	plt.plot(xdata,ydata,'b')
	for i in trigs:
		ut.vline(i,legend=False)
	if sign=='pos': 
		ut.textbox('Threshold: '+str(thresh)+'\nDelay: '+str(delay),[.05,.05])
	if sign=='neg': 
		ut.textbox('Threshold: '+str(thresh)+'\nDelay: '+str(delay),[.05,.95])
	if savepath!=False:
		ut.create_dir(savepath)
		plt.savefig(savepath+'thresh_trig_'+lab+'.png')
	
	return trigs,dtrigs
			
#------------------------------------------------------------------------------------------
def plate_v(run,**kwargs):
	
	folder=kwargs.get('folder',1)
	
	import numpy as np
	import matplotlib.pyplot as plt

	# Read in the pulsing plate waveforms
	pos=sc_read(run,1,folder=folder)
	neg=sc_read(run,2,folder=folder)
	
	posmin=np.amin(pos)
	negmax=np.amax(neg)
	
	return [posmin,negmax,np.mean([posmin,negmax])]
	
#-----------------------------------------------------------------------------------------
def pulse_count(fC,area):
	
	import numpy as np
	  
	ppa=1./(.8*(fC/1.602E-4)*(area/(2*np.pi*np.square(1.4E3))))
	atoms=np.array([.5,1,2,5,10,20])
	out=[[str(i),str(round(ppa*i))] for i in atoms]
	print(np.array(out))
	
def ions(q,pulses,area):
	import numpy as np
	return pulses*(.8*(q/1.602E-19)*(area/(2*np.pi*np.square(1.4E3))))

#------------------------------------------------------------------------------------------
def sc_rename(date):
	
	import subprocess as sp
	import numpy as np
	import glob
	import utility as utS
	
	f1=glob.glob('/home/chris/anaconda/data/'+str(date)+'/scope/*')
	f2=glob.glob('/home/chris/anaconda/data/'+str(date)+'/scope2/*')

	n1,n2=[],[]
	for i in f1:
		n1.append(int(i[-4:]))
	for i in f2:
		n2.append(int(i[-4:]))
	for i in range(len(n2)):
		ut.create_dir('/home/chris/anaconda/data/'+str(date)+'/scope/ALL'+str(n2[i]+1+np.amax(n1)).zfill(4)+'/')
		for j in glob.glob('/home/chris/anaconda/data/'+str(date)+'/scope2/ALL'+str(n2[i]).zfill(4)+'/*'):
			newpath=j.replace('scope2','scope')
			sp.run(['cp',j,newpath.replace(str(n2[i]).zfill(4),str(n2[i]+1+np.amax(n1)).zfill(4))]);
	
#------------------------------------------------------------------------------------------
def p_overlap(run,**kwargs):
	
	sbin=kwargs.get('sbin',4)
	tnewl=kwargs.get('interp_points',50.)
	rnd=kwargs.get('epsilon',4)
	bounds=kwargs.get('bounds','full')
	
	import numpy as np
	from decimal import Decimal
	import scipy.interpolate as interp
	import utility as ut
	import matplotlib.pyplot as plt
	
	t,v=sc_read(run)
	
	#slice out the relevant part
	if bounds!='full':
		bounds=ut.bound_finder(t,bounds)
		t2,v2=t[bounds[0]:bounds[1]],v[:,bounds[0]:bounds[1]]
	if bounds=='full':
		t2,v2=t,v
		
	# create new time and voltage arrays with tnewl times as many points for interpolation
	tnew=np.linspace(t2[0],t2[-1],np.shape(v2)[1]*tnewl)
	vnew=np.zeros((2,len(tnew)))
	v_sm=np.zeros((2,len(t2)))
	for i in range(2):
		v_sm[i,:]=ut.smoother(v2[i,:],sbin)
		f=interp.interp1d(t2,v_sm[i,:])
		vnew[i,:]=f(tnew)
		# round to nearest mV
		for j in range(len(tnew)):
			vnew[i,j]=Decimal(vnew[i,j]).quantize(Decimal('1e-'+str(rnd)))
	
	# look for when the interpolated curves are equal
	lint,rint=None,None
	for i in range(len(tnew)):
		# intercept from left
		if vnew[0,i]==vnew[1,i] and lint==None:
			lint=tnew[i]
		#intercept from right
		if vnew[0,-i]==vnew[1,-i] and rint==None and i>0:
			rint=tnew[-i]
	
	plt.figure('overlap')
	plt.clf()
	plt.plot(t,v[0,:],'b')
	plt.plot(t,v[1,:],'g')
	ut.vline(lint,legend=False)
	ut.vline(rint,legend=False)
	
	return rint-lint
			
#------------------------------------------------------------------------------------------
def pulse_width(run):
	
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	
	t,v=sc_read(run,chan=[1,3,4])
	dt=ut.avediff(t)
	
	lb,rb=[],[]
	for i in range(len(v[0,:])):
		if v[0,i-1]<=2.5 and v[0,i]>=2.5:
			lb=t[i-1]
		if v[0,i-1]>=2.5 and v[0,i]<=2.5:
			rb=t[i]
	
	plt.figure('pulsewidth')
	plt.clf()
	plt.plot(t,v[0,:])
	ut.vline(lb,legend=False)
	ut.vline(rb,legend=False)
	plt.xlim(lb-dt*50,rb+dt*50)
	
		
	return rb-lb

#------------------------------------------------------------------------------------------
def pulse_xsec(deflec,chg,p_0,**kwargs):
	
	cup_r=kwargs.get('cup_r',1.4)
	fix_wy=kwargs.get('wy','fit')
	
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	import funk
	
	x,y=np.linspace(-1.25*cup_r,1.25*cup_r,100),np.linspace(-1.25*cup_r,1.25*cup_r,100)
	deflec=deflec-deflec[np.where(chg==np.amax(chg))[0][0]]
	
	# Create 2d gaussian and sum a region (i.e. the cup area)
	if fix_wy=='fit':
		def cup_sig(cen,wx,wy,amp):
			sig=[]
			my,mx=np.meshgrid(x,y)
			for i in cen:
				pulse=funk.gauss2d_fwhm(x,y,amp,i,0.,wx,wy)
				mask=mx**2+my**2<=cup_r**2
				sig.append(np.sum(mask*pulse))
			return np.array(sig)
	if fix_wy!='fit':
		def cup_sig(cen,wx,amp):
			sig=[]
			my,mx=np.meshgrid(x,y)
			for i in cen:
				pulse=funk.gauss2d_fwhm(x,y,amp,i,0.,wx,fix_wy)
				mask=mx**2+my**2<=cup_r**2
				sig.append(np.sum(mask*pulse))
			return np.array(sig)
	
	# Fit to scan data
	params=ut.fitter1d(cup_sig,np.vstack((deflec,chg)),p_0)
	
	#print('main width: '+str(round(params[0],2)))
	return params
