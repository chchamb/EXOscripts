#----------------------------------------------------------------------------------------
# General use scripts
#----------------------------------------------------------------------------------------

# Script to draw a horizontal line at 0
def zero(**kwargs):

	# Define Keyword Argument
	fig=kwargs.get('figure','current')
	leg=kwargs.get('legend',False)
	bounds=kwargs.get('bounds','default')

	# Import modules
	import numpy as np
	import matplotlib.pyplot as plt
	
	if bounds=='default':
		# Make an array from one end of the axis to the other, and one of zeroes
		xlims=plt.gca().get_xlim()
	else:
		xlims=bounds
	x=np.linspace(xlims[0],xlims[1],100)	
	y=np.zeros(np.shape(x)[0])

	# Plot on the named axis
	if fig!='current':
		plt.figure(fig)
	if leg==True:
		plt.plot(x,y,'k--',linewidth=2,label='Zero')
		plt.legend()
	if leg==False:
		plt.plot(x,y,'k--',linewidth=2)
	
#----------------------------------------------------------------------------------------	
	
# Script to draw a horizontal line at some value
def hline(value,**kwargs):

	# Define Keyword Arguments
	col=kwargs.get('color','Black')
	units=kwargs.get('units','')
	fig=kwargs.get('figure','current')
	lab=kwargs.get('label','default')
	leg=kwargs.get('legend',True)

	# Import modules
	import numpy as np
	import matplotlib.pyplot as plt
	
	# Get figure info
	if fig!='current':
		plt.figure(fig)
		
	# Make an array from 0 to 1E4, and one of zeroes
	x=np.linspace(plt.xlim()[0],plt.xlim()[1],50)	
	y=np.zeros(np.shape(x)[0])
	y[:]=value
	
	label=scinot(value)+units
	if lab!='default':
		label=lab
	plt.plot(x,y,'--',color=col,label=label,linewidth=2)
	if leg==True:
		plt.legend()
	plt.draw()	
	
#----------------------------------------------------------------------------------------

# Script to draw a vertical line at some value
def vline(value,**kwargs):

	# Define Keyword Arguments
	col=kwargs.get('color','Black')
	units=kwargs.get('units','')
	fig=kwargs.get('figure','current')
	lab=kwargs.get('label','default')
	leg=kwargs.get('legend',True)

	# Import modules
	import numpy as np
	import matplotlib.pyplot as plt
	
	# Get figure info
	if fig!='current':
		plt.figure(fig)
	
	label=scinot(value)+units
	if lab!='default':
		label=lab
	plt.axvline(value,color=col,ls='dashed',lw=2,label=label)
	if leg==True:
		plt.legend()
	plt.draw()	
	
#----------------------------------------------------------------------------------------


# Script to analyze labview data
def lv_ave(date,**kwargs):

	# Define keyword arguments
	df=kwargs.get('darkframe',False)
	s=kwargs.get('save',True)
	rep=kwargs.get('repeats',[])
	ratio=kwargs.get('ratio',1)
	f1=kwargs.get('f1_ratio',False)
	
	# Import Modules
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	from decimal import Decimal
	
	# Make a directory to store plots and averages in
	create_dir('/home/chris/anaconda2/plots/'+str(date)+'/lv')
	# Extract all the filenames in the folder .../date/lv
	filenames=os.listdir('/home/chris/anaconda2/data/'+str(date)+'/lv')
	filenames.sort()
	# Add the full path onto the front
	filepaths=[]
	for name in filenames:
		if name=='ave_pows.txt':
			filenames.remove(name)
			continue
		if name[-3:]=='txt':
			filepaths.append('/home/chris/anaconda2/data/'+str(date)+'/lv/'+name)
			 
	# Get lv number from file name
	lv_nums=[]
	for path in filepaths:
		j=1                                                        # j is the incrementor for characters in the run name/number
		for i in range(1,len(path)+1):
			lv=''
			if path[-i]=='/':                                      # Steps backwards in path until it finds a '/'
				while path[-i+j]!='_':
					lv += path[-i+j]                               # Adds on characters after the '/' until it finds a '_'
					j+=1
				break
			else:
				continue
		lv_nums.append(lv)
	lv_nums=np.array(lv_nums)
	
	# Insert new number for repeat files into lv_nums
	for i in rep:
		lv_nums=np.insert(lv_nums,np.where(lv_nums==str(i))[0]+1,int(i*10+1))
	
	# Import data from file(s)
	av_pow=np.zeros(len(filepaths)+len(rep))
	f1ratio=[]
	for i in filepaths:
		datafile=open(i,'r')
		frame=[]
		power=[]
		n=1
		for line in datafile:
			line=line.strip()
			columns=line.split(';')
			frame.append(int(columns[0]))
			power.append(1000*ratio*float(columns[1])) # 1000 converts to mW
		if lv_nums[filepaths.index(i)] in rep:
			n=2
		for j in range(n):
			power_part=power[j*(len(power)/n):(j+1)*(len(power)/n)]
			if df==True:
				av_pow[filepaths.index(i)+j]=np.mean(power_part[1:])
			if df==False:
				av_pow[filepaths.index(i)+j]=np.mean(power_part)
			# Plot the power for each frame
			plt.figure('lv')
			plt.clf()
			plt.plot(frame,power_part,'bo',label='Measured Power')
			plt.xlabel('Frame')
			plt.ylabel('Measured Power (mW)')
			plt.title(str(date)+': lv '+str(lv_nums[filepaths.index(i)+j]))
			avestr='%.2E' % Decimal(av_pow[filepaths.index(i)+j])
			h_line(av_pow[filepaths.index(i)+j],label='Average: '+avestr+'mW')
			plt.axis([0,frame[-1]+1,0,1.1*np.amax(power_part)])
			plt.legend(loc=4)
			plt.savefig('/Users/christopherchambers/Py/plots/'+str(date)+'/lv/lv'+str(lv_nums[filepaths.index(i)+j])+'.png')
			plt.draw()
			if f1==True:
				f1ratio.append(np.mean(power[1:])/power[0])
		
	if s==True:
		# Make a file to write to
		open('/Users/christopherchambers/Py/plots/'+str(date)+'/lv/ave_pows.txt','w')
		import asciitable
		# Make into a structured numpy array
		if ratio==1:
			p_type='Ave Power (ref)'
		else:
			p_type='Ave Power (trans)'
		data=np.zeros(len(av_pow),dtype=[('lv','S4'),(p_type,'float64')])
		data[:]=zip(lv_nums,av_pow)
		asciitable.write(data,'/Users/christopherchambers/Py/plots/'+str(date)+'/lv/ave_pows.txt')
		
	if f1==True:
		return lv_nums,av_pow,f1ratio
	else:
		return lv_nums,av_pow

#----------------------------------------------------------------------------------------

# Script to read temp data from fringe_temp.vi
def temp_read(filepath,**kwargs):

	# Define Keyword Arguments
	col=kwargs.get('color','Black')
	tag=kwargs.get('tag','')
	bounds=kwargs.get('bounds','all')
	o=kwargs.get('overlay',False)
	
	# Import Modules
	import numpy as np
	import matplotlib.pyplot as plt
	
	# Open the file and parse data
	datafile=open(filepath,'r')
	timestamp,temp=[],[]
	for line in datafile:
		line=line.strip()
		columns=line.split()
		timestamp.append(columns[4])
		temp.append(float(columns[5]))
	
	# Convert timestamp to seconds
	hms=[3600,60,1]
	time_raw=np.zeros(len(timestamp))
	for i in range(len(timestamp)):
		time_raw[i]=sum(a*b for a,b in zip(hms,map(float,timestamp[i].split(':'))))
	# Set first time to t=0
	time_0=time_raw[0]
	time= time_raw-time_0
	dt=np.zeros(np.shape(time)[0]-1)
	for i in range(1,np.shape(time)[0]):
		dt[i-1]=time[i]-time[i-1]
	dt=np.mean(dt)
	 
	# Find the index of the time closest to bounds
	bound_index=np.zeros(np.shape(bounds),dtype='int')
	if bounds=='all':
		bound_index=[0,len(time)]	
	else:	
		for k in range(2):
			j=-1
			for i in time[1:]:
				j=j+1
				if i<bounds[k]:
					last_time=i
				else:
					a=abs(bounds[k]-last_time)
					b=abs(i-bounds[k])
					if a>b:
						bound_index[k]=j
						break
					else:
						bound_index[k]=j-1
						break
	N=bound_index[1]-bound_index[0]
	
	# Plot time vs. temp on current axis
	if o==False:
		plt.plot(time,temp,color=col,label=tag)
	if o==True:
		plt.plot(np.arange(0,N*dt,dt),temp[bound_index[0]:bound_index[1]],color=col,label=tag)
	plt.xlabel('Time(s)')
	plt.ylabel('Temperature (K)')
	plt.show()
	
#----------------------------------------------------------------------------------------
# Script to plot fringe data from lv files

def fringe_plot(filenum,cols,**kwargs):

	# Define keyword Arguments
	TR=kwargs.get('plot','both')
	store=kwargs.get('store',False)
	# Import Modules
	import numpy as np
	import matplotlib.pyplot as plt
	
	# Get filepath and open file
	filepath='/Users/christopherchambers/Py/data/'+str(date)+'/lv/fringes_'+str(filenum)+'.txt'
	datafile=open(filepath,'r')
	
	# Read the Data from File
	timestamp,ref,trans=[],[],[]
	for line in datafile:
		line=line.strip()
		columns=line.split()
		timestamp.append(columns[1])
		ref.append(float(columns[2]))
		trans.append(float(columns[5]))
	datafile.close()
	ref=np.array(ref)
	trans=np.array(trans)
	
	# Convert timestamp to time elapsed
	hms=[3600,60,1]
	time_raw=np.zeros(len(timestamp))
	for i in range(len(timestamp)):
		time_raw[i]=sum(a*b for a,b in zip(hms,map(float,timestamp[i].split(':'))))
	# Set first time to t=0
	time_0=time_raw[0]
	time= time_raw-time_0
	
	# Plot on current axis
	if TR=='both': 
		plt.plot(time,ref,color=cols[0],label='fringes_'+repr(filenum)+' (R)')
		plt.plot(time,trans,color=cols[1],label='fringes_'+repr(filenum)+' (T)')
		plt.axis([0,time[-1],0,1.1*np.amax(np.concatenate((ref,trans)))])
	if TR=='R':
		plt.plot(time,ref,color=cols[0],label='fringes_'+repr(filenum)+' (R)')
		plt.axis([0,time[-1],0,1.1*np.amax(ref)])
	if TR=='T':
		plt.plot(time,trans,color=cols[0],label='fringes_'+repr(filenum)+' (T)')
		plt.axis([0,time[-1],0,1.1*np.amax(trans)])
	plt.legend()
	plt.show()
	
	if store==True:
		data=np.zeros((2,len(time)))
		data[0,:],data[1,:]=ref,trans
		return time,data
	
#----------------------------------------------------------------------------------------

# Script to fit a function to data
# [data] should be in 2 rows, row 0: x values, row 1: y values
def fitter1d(func,data,param_0,**kwargs):

	# Define Keyword Arguments
	col=kwargs.get('color','Black')
	fitlim=kwargs.get('fitlim','full')
	fig=kwargs.get('fig','default')
	disp_param=kwargs.get('d_param',False)
	savepath=kwargs.get('savepath',False)
	label=kwargs.get('label','')
	numpts=kwargs.get('numpts',5)
	
	# Import Modules
	import numpy as np
	import matplotlib.pyplot as plt
	from scipy.optimize import curve_fit
	import utility as ut
	
	# Slice data to get frames/region of interest
	if fitlim!='full':
		# Find index for upper fit limit
		j=-1
		for i in data[0,:]:
			j=j+1
			if i<fitlim:
				last_wave=i
			else:
				a=abs(fitlim-last_wave)
				b=abs(i-fitlim)
				if a>b:
					lim_index=j
					break
				else:
					lim_index=j-1
					break	
		data=data[:,0:lim_index]
	
	# Do the fit
	param,sigma=curve_fit(func,data[0,:],data[1,:],p0=param_0)
	
	# Make a more finely spaced x data for fit plot
	xfit=np.linspace(np.amin(data[0,:]),np.amax(data[0,:]),numpts*np.shape(data)[1])
	# Plot fit and data
	if fig=='default':
		plt.figure('Fitter')
	else:
		plt.figure(fig)
	plt.clf()
	plt.plot(data[0,:],data[1,:],'bo',label='Data')
	# Dynamically create a plot command with as many parameters as needed
	plt_str='plt.plot(xfit,func(xfit,'
	p_txt=''
	for i in range(len(param)):
		if i==0:
			plt_str += 'param[%d]' %i
			if disp_param!=False:
				p_txt+=disp_param[i]+': '+scinot(param[i])
		else:
			plt_str += ',param[%d]' %i
			if disp_param!=False:
				p_txt+='\n'+disp_param[i]+': '+scinot(param[i])
	plt_str += '),"--",color=col,linewidth=2,label="Model")'
	exec(plt_str)
	if disp_param!=False:
		textbox(p_txt,[.75,.8])
	plt.legend()
	if savepath!=False:
		ut.create_dir(savepath)
		plt.savefig(savepath+'fit_'+label+'.png')
	
	return param

#-----------------------------------------------------------------------------------------
	
def gaussfit2d(data,x,y,p0,**kwargs):
	
	plot=kwargs.get('plot',True)
	
	import numpy as np
	from lmfit import Parameters,Parameter,minimize,Model
	import funk
	import matplotlib.pyplot as plt
	
	def gauss2d_flat(x,y,amp,xcenter,ycenter,xwidth,ywidth):  # This ends up as f(y,x)
		my,mx=np.meshgrid(x,y)
		xwidth,ywidth=xwidth/2.,ywidth/2.	
		return np.ravel(amp*np.exp(-2*(np.square(mx-ycenter)/(np.square(ywidth))))*np.exp(-2*(np.square(my-xcenter))/(np.square(xwidth))))
	
	# Define x,y grid and p_0 names
	p_0s=['amp','xcenter','ycenter','xwidth','ywidth']
	func_model=Model(gauss2d_flat,independent_vars=['x','y'],param_names=p_0s)
	func_model.set_param_hint('amp',value=p0[0],min=0)
	func_model.set_param_hint('xcenter',value=p0[1])
	func_model.set_param_hint('ycenter',value=p0[2])
	func_model.set_param_hint('xwidth',value=p0[3],min=0)
	func_model.set_param_hint('ywidth',value=p0[4],min=0)
	result=func_model.fit(np.ravel(data),x=x,y=y,verbose=False)
	
	if plot==True:
		plt.figure('best_fit')
		plt.clf()
		plt.imshow(np.flipud(np.reshape(result.eval(),np.shape(data))),interpolation='nearest',extent=[x[0],x[-1],y[0],y[-1]])
		plt.colorbar(orientation='vertical')
		
	return result.best_values
	
#----------------------------------------------------------------------------------------

def gaussfit2d_multi(data,x,y,p0,**kwargs):
	
	plot=kwargs.get('plot',True)
	
	import numpy as np
	from lmfit import Parameters,Parameter,minimize,Model
	import funk
	import matplotlib.pyplot as plt
	
	def gauss2d_flat(x,y,amp,xcenter,ycenter,xwidth,ywidth):  # This ends up as f(y,x)
		my,mx=np.meshgrid(x,y)
		xwidth,ywidth=xwidth/2.,ywidth/2.	
		return np.ravel(amp*np.exp(-2*(np.square(mx-ycenter)/(np.square(ywidth))))*np.exp(-2*(np.square(my-xcenter))/(np.square(xwidth))))
	def zero(x,y,a):
		return np.ravel(a*np.zeros((np.shape(x)[0],np.shape(y)[0])))
	
	n=len(p0[0])
	fullmod=Model(zero,independent_vars=['x','y'],param_names='a')
	fullmod.set_param_hint('a',value=0,max=.001,min=-.001)
	for i in range(n):
		# Define x,y grid and p_0 names
		p_0s=['amp','xcenter','ycenter','xwidth','ywidth']
		func_model=Model(gauss2d_flat,independent_vars=['x','y'],param_names=p_0s,prefix='g'+str(i+1)+'_')
		fullmod+=func_model
		fullmod.set_param_hint('g'+str(i+1)+'_amp',value=p0[0][i],min=0)
		fullmod.set_param_hint('g'+str(i+1)+'_xcenter',value=p0[1][i])
		fullmod.set_param_hint('g'+str(i+1)+'_ycenter',value=p0[2][i])
		fullmod.set_param_hint('g'+str(i+1)+'_xwidth',value=p0[3][i],min=0)
		fullmod.set_param_hint('g'+str(i+1)+'_ywidth',value=p0[4][i],min=0)
		
	result=fullmod.fit(np.ravel(data),x=x,y=y)
	if plot==True:
		plt.figure('best_fit')
		plt.clf()
		plt.imshow(np.flipud(np.reshape(result.eval(),np.shape(data))),interpolation='nearest',extent=[x[0],x[-1],y[0],y[-1]])
		plt.colorbar(orientation='vertical')
	
	return result.best_values,np.flipud(np.reshape(result.eval(),np.shape(data)))
	
#----------------------------------------------------------------------------------------
# Script to save to ascii
def save(data,names,outfile):

	# Import modules
	import numpy as np
	import asciitable
	import os
	
	# Make a file to write the data to
	open(outfile,'a')

	# Create dictionary of 'name':data
	data_dict={}
	if len(names)==1:
		data_dict[names[0]]=data
	if len(names)>=2:
		for i in range(len(names)):
				if len(np.shape(data))>1:
					data_dict[names[i]]=data[i,:]
				if len(np.shape(data))==1:
					data_dict[names[i]]=[data[i]]
	# Write data with asciitable
	asciitable.write(data_dict,outfile,names=names)
		
#----------------------------------------------------------------------------------------

# Script to save to ascii
def spec_save(data,outfile):

	# Import modules
	import numpy as np
	import asciitable
	import os
	
	# Make a file to write the data to
	open(outfile,'a')
	
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
	# Rearrange so that each frame is column, instead of row
	datalist=[]
	names=[]
	for i in range(np.shape(data)[0]):
		names.append('Frame '+repr(i+1))
	for j in range(np.shape(data)[1]):
		datalist.append(list(data[:,j]))
		
	# Write data with asciitable
	asciitable.write(datalist,outfile,names=names)

#-----------------------------------------------------------------------------------------	

# Script to remove line from a plot	
def rem(figname,typ,num):

	# Import Modules
	import matplotlib.pyplot as plt
	
	# Remove and redraw the plot
	if figname!='c':
		fig=plt.figure(figname)
	ax=plt.gca()
	if typ=='line':
		ax.lines.pop(num-1)
	if typ=='patch':
		ax.patches.pop(num-1)
	if typ=='text':
		ax.texts.pop(num-1)
	plt.draw()

#-----------------------------------------------------------------------------------------
# Script to change figure to log scale

def log():

	import matplotlib.pyplot as plt
	
	plt.gca().set_yscale('log')
	plt.draw()
	
#-----------------------------------------------------------------------------------------

# Script to change figure to log scale

def loglog():

	import matplotlib.pyplot as plt
	
	plt.gca().set_yscale('log')
	plt.gca().set_xscale('log')
	plt.draw()

#-----------------------------------------------------------------------------------------
# Script to change figure to log scale in x

def logx():

	import matplotlib.pyplot as plt
	
	plt.gca().set_xscale('log')
	plt.draw()
		
#-----------------------------------------------------------------------------------------

# Script to change figure to linear scale

def lin():

	import matplotlib.pyplot as plt
	
	plt.gca().set_yscale('linear')
	plt.draw()
	
#-----------------------------------------------------------------------------------------
	
# Script to change scale type to scientific
def sci(**kwargs):
	
	ax=kwargs.get('axis','y')
	
	import matplotlib.pyplot as plt
	
	plt.ticklabel_format(style='sci',axis=ax,scilimits=(0,0))
	plt.draw()	
	
	
#-----------------------------------------------------------------------------------------

# Script to convert to good unit prefix

def prefix(num):
	
	import numpy as np
	if num==0:
		pow=0
	else:
		pow=np.floor(np.log10(num))
	
	if -15>pow>=-18:
		fix='a'
		pow2=18
	if -12>pow>=-15:
		fix='f'
		pow2=15
	if -9.>pow>=-12:
		fix='p'
		pow2=12
	if -6>pow>=-9:
		fix='n'
		pow2=9
	if -3>pow>=-6:
		fix='u'
		pow2=6
	if 0>pow>=-3:
		fix='m'
		pow2=3
	if 3>pow>=0:
		fix=''
		pow2=0
	if 6>pow>=3:
		fix='k'
		pow2=-3
	if 9>pow>=6:
		fix='M'
		pow2=-6
	if 12>pow>=9:
		fix='G'
		pow2=-9
	if pow>=12:
		fix='T'
		pow2=-12
		
	return pow2,fix 		
		
#-----------------------------------------------------------------------------------------

def conv(num):

	import numpy as np
	num2=np.absolute(num)
	pow,fix=prefix(num2)
	if pow>=0:
		num_c=num*np.power(10,pow)
	if pow<0:
		num_c=num/np.power(10,np.absolute(pow))
	
	return str(round(num_c,2))+' '+fix
	
#-----------------------------------------------------------------------------------------

def scinot(num):
	
	from decimal import Decimal
	
	return '%.2E' % (num)

#------------------------------------------------------------------------------------------
# Find indices for bounds

def bound_finder(arr,bnds):

	bnd_index=[]
	for k in range(len(bnds)):
		j=-1
		for i in arr:
			j=j+1
			if i<bnds[k]:
				last_wave=i
			else:
				a=abs(bnds[k]-last_wave)
				b=abs(i-bnds[k])
				if a>b:
					bnd_index.append(j)
					break
				else:
					bnd_index.append(j-1)
					break
	return bnd_index	
	
#-----------------------------------------------------------------------------------------

# Find indices for bounds

def bound_finder_neg(arr,bnds):

	bnd_index=[]
	for k in range(len(bnds)):
		j=-1
		for i in arr:
			j=j+1
			if i>bnds[k]:
				last_wave=i
			else:
				a=abs(bnds[k]-last_wave)
				b=abs(i-bnds[k])
				if a>b:
					bnd_index.append(j)
					break
				else:
					bnd_index.append(j-1)
					break
	return bnd_index	
	
#-----------------------------------------------------------------------------------------
# Make directory(ies) if it(they) doesn't exist

def create_dir(savefile,**kwargs):
	
	subdirs=kwargs.get('subdirs','none')
	import os
	
	if os.path.exists(savefile)==False:
		os.makedirs(savefile)
	
	if subdirs!='none':
		
		for i in subdirs:
			if os.path.exists(savefile+'/'+i)==False:
				os.makedirs(savefile+'/'+i)
#-----------------------------------------------------------------------------------------
# Remove directory(ies) if it(they) exist

def rm_dir(savefile,**kwargs):
	
	subdirs=kwargs.get('subdirs','none')
	import os
	import glob
	
	if os.path.exists(savefile)==True:
		for i in glob.glob(savefile+'*'):
			os.remove(i)
	
	if subdirs!='none':
		
		for i in subdirs:
			if os.path.exists(savefile+'/'+i)==False:
				os.makedirs(savefile+'/'+i)
#-----------------------------------------------------------------------------------------
# Shift 2d array by fraction of a pixel (re-binned) note: x is axis 1 due to how things are read

def frac_shift(data,xshift,yshift,sign):
	
	import numpy as np
	
	if sign=='np':
		xshift=-1*xshift
	if sign=='pn':
		yshift=-1*yshift
	if sign=='nn':
		xshift=-1*xshift
		yshift=-1*yshift
	
	# Shift each row in y
	if yshift==0:
		data_new=data
	if yshift!=0:
		data_new=np.zeros(np.shape(data))
		for j in range(np.shape(data)[1]):
			if yshift>0:
				data_new[0,j]=data[0,j]*(yshift)
				for i in np.arange(1,np.shape(data)[0]):
					data_new[i,j]=data[i,j]*yshift+data[i-1,j]*(1-yshift)
			if yshift<0:
				yshifta=np.absolute(yshift)
				data_new[-1,j]=data[-1,j]*(1-yshifta)
				for i in np.arange(np.shape(data)[0]-1):
					data_new[i,j]=data[i+1,j]*yshifta+data[i,j]*(1-yshifta)
	# Shift each column in x
	if xshift==0:
		data_new2=data_new
	if xshift!=0:
		data_new2=np.zeros(np.shape(data_new))
		for j in range(np.shape(data_new)[0]):
			if xshift>0:
				data_new2[j,0]=data_new[j,0]*(xshift)
				for i in np.arange(1,np.shape(data)[1]):
					data_new2[j,i]=data_new[j,i]*xshift+data_new[j,i-1]*(1-xshift)
			if xshift<0:
				xshifta=np.absolute(xshift)
				data_new2[j,-1]=data_new[j,-1]*(1-xshifta)
				for i in np.arange(np.shape(data)[1]-1):
					data_new2[j,i]=data_new[j,i+1]*xshifta+data_new[j,i]*(1-xshifta)
	
	return data_new2

#-----------------------------------------------------------------------------------------
# Shift 2d array by integer pixels (no re-binning) note: x is axis 1 due to how things are read

def int_shift(data,xshift,yshift,**kwargs):
	
	val=kwargs.get('val',0)
	
	import numpy as np
	
	data2=data
	xzeros=np.zeros((np.shape(data)[0],1))
	yzeros=np.zeros((np.shape(data)[1]))
	xzeros[:],yzeros[:]=val,val
	
	# Shift in x
	if xshift>0:
		for i in range(int(xshift)):
			data=np.hstack((xzeros,data))
			data=np.delete(data,-1,axis=1)
	if xshift<0:
		for i in range(int(np.absolute(xshift))):
			data=np.hstack((data,xzeros))
			data=np.delete(data,0,axis=1)
	# Shift in y
	if yshift>0:
		for i in range(int(yshift)):
			data=np.vstack((yzeros,data))
			data=np.delete(data,-1,axis=0)
	if yshift<0:
		for i in range(int(np.absolute(yshift))):
			data=np.vstack((data,yzeros))
			data=np.delete(data,0,axis=0)
	
	return data
#------------------------------------------------------------------------------------------				

# Script to combine both shifts

def shifter(data,xshift,yshift,sign):

	import numpy as np
	
	# Break up shift into int and frac
	xshift_i=np.floor(xshift)
	xshift_f=xshift-xshift_i
	yshift_i=np.floor(yshift)
	yshift_f=yshift-yshift_i
	
	# Do the shifts (int first since it does not re-bin)
	data=int_shift(data,xshift_i,yshift_i)
	data=frac_shift(data,xshift_f,yshift_f,sign)
	
	return data

#------------------------------------------------------------------------------------------
	
# Calculate beam width at some z (z units variable, w_0 in micron, wlenght in nm)
def beam_w(z,w_0,wlength,**kwargs):
	
	import numpy as np
	
	zunits=kwargs.get('units','inches')
	n=kwargs.get('n',1.)
	
	# Convert z to meters
	if zunits=='inches':
		z=z*(25.4/1000)
	if zunits=='mm':
		z=z*1E-3
	if zunits=='micron':
		z=z*1E-3
	wlength=wlength*1E-9
	w_0=w_0*1E-6
		
	z_0=(np.pi*np.square(w_0)*n)/(wlength)
	return conv(w_0*np.sqrt(1+np.square(z/z_0)))

#------------------------------------------------------------------------------------------

# Draw a rectangle
def rect(bounds,**kwargs):
	
	import matplotlib.pyplot as plt
	import matplotlib.patches as pa
	
	figure=kwargs.get('fig','current')
	bndtype=kwargs.get('boundtype','edge')
	
	# Make bounds
	if bndtype=='edge':
		[xl,xh,yl,yh]=bounds
	if bndtype=='center':
		xl=bounds[0]-.5*bounds[2]
		xh=bounds[0]+.5*bounds[2]
		yl=bounds[1]-.5*bounds[3]
		yh=bounds[1]+.5*bounds[3]
	# Axes definition from plotter puts integer in middle of bin, so have to shift back by .5
	xl,xh,yl,yh=xl-.5,xh-.5,yl-.5,yh-.5
	
	if figure!='current':
		fig=plt.figure(figure)
	
	area=plt.gca().add_patch(pa.Rectangle((xl,yl),xh-xl,yh-yl,ls='dashed',fill=False,ec='k',lw=2))
	
	plt.draw()

#------------------------------------------------------------------------------------------

def textbox(txt,loc,**kwargs):
	
	f=kwargs.get('fig','current')
	fontsize=kwargs.get('fontsize',10)
	color=kwargs.get('color','k')
	
	import matplotlib.pyplot as plt
	
	if f!='current':
		plt.figure(f)
	ax=plt.gca()
	bbox_props = dict(boxstyle='square', fc='.9', ec='k', alpha=1.0)
	plt.text(loc[0],loc[1],txt,fontsize=fontsize,color=color,bbox=bbox_props,va='top',ha='left',transform=ax.transAxes)		

#-----------------------------------------------------------------------------------------
# Fourier Transform

def fourier(time,data,**kwargs):

	# Define keywords
	xlim=kwargs.get('xlim','default')
	smooth=kwargs.get('smooth',False)
	store=kwargs.get('store',True)

	import numpy as np
	import scipy.fftpack as fourier
	import matplotlib.pyplot as plt
	import gen_reader as gr
	import utility as ut
	
	# Smooth the data
	if smooth!=False:
		data=smoother(data,smooth)
	
	# Determine the spacing of the time domain
	time=time-time[0]
	dt,j=0.,0
	for i in range(1,np.shape(time)[0]):
		dt+=time[i]-time[i-1]
		j+=1
	dt=dt/j
	
	# Construct frequency domain (in hertz)
	frq=np.linspace(0,1./(2*dt),np.shape(time)[0])
	#frq[1:]=1./frq[1:]
	
	# Do the Fourier transform
	fftfunc=np.abs(fourier.rfft(data))
	
	# Plot the Fourier Transform and cutoff
	plt.figure('fourier')
	plt.clf()
	plt.plot(frq[1:],fftfunc[1:],'b-',label='Fourier Transform')
	plt.xlabel('Frequency (Hz)')
	plt.ylabel('Fourier Transform')
	if xlim!='default':
		plt.xlim(xlim)
	plt.ylim(1.1*np.amin(fftfunc[1:]),1.1*np.amax(fftfunc[1:]))
	ut.logx()
	plt.title('Fourier Transform')
	
	
	if store==True:
		return frq,fftfunc

#-----------------------------------------------------------------------------------------

# Inverse Fourier Transform

def ifourier(frq,data,**kwargs):

	# Define keywords
	xlim=kwargs.get('xlim','default')
	smooth=kwargs.get('smooth',False)
	store=kwargs.get('store',False)

	import numpy as np
	import scipy.fftpack as fourier
	import matplotlib.pyplot as plt
	import gen_reader as gr
	import utility as ut
	
	# Smooth the data
	if smooth!=False:
		time,data=smoother(time,data,smooth)
	
	# Determine the spacing of the frequency domain
	df,j=0.,0
	for i in range(1,np.shape(frq)[0]):
		df+=frq[i]-frq[i-1]
		j+=1
	df=df/j
	
	# Construct time domain in seconds
	time=np.linspace(0,2*df,np.shape(frq)[0])
	
	# Do the Fourier transform
	ifftfunc=np.abs(fourier.rfft(data))
	
	# Plot the Fourier Transform and cutoff
	plt.figure('ifourier')
	plt.clf()
	plt.plot(time[1:],ifftfunc[1:],'b-',label='Inverse Fourier Transform')
	plt.xlabel('Time (s)')
	plt.ylabel('Inverse Fourier Transform')
	if xlim!='default':
		plt.xlim(xlim)
	plt.ylim(1.1*np.amin(ifftfunc[1:]),1.1*np.amax(ifftfunc[1:]))
	plt.title('Inverse Fourier Transform')
	
	if store==True:
		return time,ifftfunc,

#-----------------------------------------------------------------------------------------

# Script to smooth data (moving average)
def smoother(ydata,binw,**kwargs):

	savefile=kwargs.get('savefile',False)
	
	import numpy as np
	import matplotlib.pyplot as plt

	# Do a moving average
	ysmooth=np.zeros(np.shape(ydata))
	for i in range(len(ydata)):
		if i<binw:
			ysmooth[i]=np.mean(ydata[:i+binw+1])
		if i>=np.shape(ydata)[0]-(binw+1):
			ysmooth[i]=np.mean(ydata[i-binw:])
		if binw<=i<np.shape(ydata)[0]-(binw+1):
			ysmooth[i]=np.mean(ydata[i-binw:i+binw+1])
		
	# Plot original and smoothed
	plt.figure('smoother')
	plt.clf()
	plt.plot(ydata,'b-',label='Raw')
	plt.plot(ysmooth,'k--',dashes=(2,2),lw=2,label='Smoothed: bin size= '+str(2*binw))
	plt.legend(fontsize=12)
	plt.title('Smoothed Data')
	if savefile!=False:
		plt.savefig(savefile)
	
	return ysmooth
 #-----------------------------------------------------------------------------------------
def cbar():
	
	import matplotlib.pyplot as plt
	
	cb=plt.colorbar(orientation='vertical')
	cb.formatter.set_powerlimits((1,1))
	cb.update_ticks
	
	plt.draw()

#------------------------------------------------------------------------------------------
# Plot contour(s)
def cont(array,contours,**kwargs):
	
	fig=kwargs.get('figure','current')
	custlab=kwargs.get('custom_label',False)
	fntsz=kwargs.get('fontsize',12)
	xdata=kwargs.get('xdata','default')
	ydata=kwargs.get('ydata','default')
	cols=kwargs.get('colors',['black'])
	
	import numpy as np
	import matplotlib.pyplot as plt
	
	if fig!='current':
		plt.figure(fig)
	
	contours_real=np.array(contours)*np.amax(array)
	if xdata=='default':
		cs=plt.contour(array,levels=contours_real,colors=cols,linewidths=[2],linestyles=['dashed'])
	else:
		cs=plt.contour(xdata,ydata,array,levels=contours_real,colors=cols,linewidths=[2],linestyles=['dashed'])
		
	fmt={}
	if custlab==False:
		contours_str=str(contours)
	else:
		contours_str=custlab
	for i,j in zip(cs.levels,contours_str):
		fmt[i]=j
	plt.clabel(cs,inline=1,fmt=fmt,fontsize=fntsz,manual=True)
	plt.draw()
	
#------------------------------------------------------------------------------------------
# Make arrays for data entry
def data_lists(date):
	
	import numpy as np
	import gen_reader as gr
	
	header,full_list=gr.reader('/home/chris/anaconda/data/'+str(date)+'/run_data.csv',dtype=int,delimeter=',')
	
	runlist=full_list[:,0]
	lvlist=full_list[:,1]
	scplist=np.vstack((runlist,full_list[:,2]))
	pulselist=full_list[:,3]
	
	return runlist,lvlist,scplist,pulselist
#-----------------------------------------------------------------------------------------
# Lowpass Filter

def lowpass(t,func,cutoff):

	import numpy as np
	import scipy.fftpack as fourier
	import matplotlib.pyplot as plt
	
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
	logx()
	
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
# Highpass Filter

def highpass(t,func,cutoff):

	import numpy as np
	import scipy.fftpack as fourier
	import matplotlib.pyplot as plt
	
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
	logx()
	
	# Apply the cut
	for i in range(np.shape(fftfunc)[0]):
		if frq[i]<cutoff:
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
# timestamp to time array
def timestamp(array,**kwargs):
	
	t0=kwargs.get('t0',False)
	
	import numpy as np
	
	time=[]
	for i in array:
		hms=i.split(':')
		time.append(float(hms[0])*3600+float(hms[1])*60+float(hms[2]))
	if t0==True:
		return np.reshape(np.array(time),(len(time),1))-time[0]
	if t0==False:
		return np.array(time)
 
 #-----------------------------------------------------------------------------------------
 
def bool_read(array,delay):
	 
	import numpy as np
	 
	time=[]
	j=0
	for i in range(np.shape(array)[1]):
		if i<=j+delay:
			continue
		if array[1,i]==1.:
			time.append(array[0,i])
			j=i	
	gap=[]
	for i in range(len(time)-1):
		gap.append(time[i+1]-time[i])
	return time,gap		 
 
#------------------------------------------------------------------------------------------

def e_fold_area(data,x,y):
	
	import numpy as np
	
	above=0
	for i  in np.ravel(data):
		if i>=np.amax(data)*np.exp(-1):
			above+=1
	dx,dy=x[1]-x[0],y[1]-y[0]
	area=above*dx*dy
	
	return area

#------------------------------------------------------------------------------------------
def multi_plot(x,y1,y2,**kwargs):
	
	styles1=kwargs.get('styles1','default')
	styles2=kwargs.get('styles2','default')
	leg=kwargs.get('legend',False)
	ylims=kwargs.get('ylims',False)
	ylab=kwargs.get('ylabels',False)
	
	import matplotlib.pyplot as plt
	
	cols=['b','r','g','c','m','k']
	if styles1=='default':
		styles1=[]
		for i in range(len(y1)):
			styles1.append(cols[i]+'-')
	if styles2=='default':
		styles2=[]
		for i in range(len(y2)):
			styles2.append(cols[i]+'-')
	
	fig=plt.gcf()
	plt.clf()
	ax1=fig.add_subplot(111)
	for i in range(len(y1)):
		ax1.plot(x,y1[i],styles1[i])
		if leg!=False:
			plt.legend(leg[i],loc=2,numpoints=1,fontsize=12)
		if ylab!=False:
			plt.ylabel(ylab[0])
		if ylims!=False:
			plt.ylim(ylims[0],ylims[1])
	ax2=ax1.twinx()
	for i in range(len(y2)):
		ax2.plot(x,y2[i],styles2[i])
		if leg!=False:
			plt.legend(leg[i],loc=1,numpoints=1,fontsize=12)
		if ylab!=False:
			plt.ylabel(ylab[1])
		if ylims!=False:
			plt.ylim(ylims[2],ylims[3])

#------------------------------------------------------------------------------------------
def integrate(array,bounds,**kwargs):
	
	bndtype=kwargs.get('boundtype','center')
	
	import numpy as np
	
	# Integrate over the bounds
	if bndtype=='edge':
		xbounds=np.array(bounds[:2])
		ybounds=np.array(bounds[-2:])
	if bndtype=='center':
		xbounds=np.array([bounds[0]-(bounds[2]/2.),bounds[0]+(bounds[2]/2.)])
		ybounds=np.array([bounds[1]-(bounds[3]/2.),bounds[1]+(bounds[3]/2.)])
	xbounds_i=np.floor(xbounds)
	ybounds_i=np.floor(ybounds)
	xbounds_f=xbounds-np.floor(xbounds)
	ybounds_f=ybounds-np.floor(ybounds)
	fullpix=np.sum(array[ybounds_i[0]+1:ybounds_i[1],xbounds_i[0]+1:xbounds_i[1]])
	sidepixl=(1-xbounds_f[0])*np.sum(array[ybounds_i[0]+1:ybounds_i[1],xbounds_i[0]])
	sidepixr=xbounds_f[1]*np.sum(array[ybounds_i[0]+1:ybounds_i[1],xbounds_i[1]])
	sidepixb=(1-ybounds_f[0])*np.sum(array[ybounds_i[0],xbounds_i[0]+1:xbounds_i[1]])
	sidepixt=ybounds_f[1]*np.sum(array[ybounds_i[1],xbounds_i[0]+1:xbounds_i[1]])
	cornerbl=(1-xbounds_f[0])*(1-ybounds_f[0])*array[ybounds_i[0],xbounds_i[0]]
	cornertl=(1-xbounds_f[0])*ybounds_f[1]*array[ybounds_i[1],xbounds_i[0]]
	cornertr=xbounds_f[1]*ybounds_f[1]*array[ybounds_i[1],xbounds_i[1]]
	cornerbr=xbounds_f[1]*(1-ybounds_f[0])*array[ybounds_i[0],xbounds_i[1]]
	
	integral=fullpix+sidepixl+sidepixr+sidepixt+sidepixb+cornerbl+cornertl+cornerbr+cornertr
	
	return integral
	
#------------------------------------------------------------------------------------------

def t_over_r(date):
	
	import gen_reader as gr
	import numpy as np
	
	T=gr.wf_reader('/home/chris/anaconda/data/'+str(date)+'/lv/T.txt',header=False)
	R=gr.wf_reader('/home/chris/anaconda/data/'+str(date)+'/lv/R.txt',header=False)
	
	return np.mean(T[:,1])/np.mean(R[:,1])	

#------------------------------------------------------------------------------------------
def avediff(array):
	
	import numpy as np
	
	diff=[]
	for i in range(len(array)-1):
		diff.append(array[i+1]-array[i])
	
	return np.mean(diff)

#------------------------------------------------------------------------------------------
def max_finder(arr,xarr,bnds):
	
	import numpy as np
	
	if bnds=='all':
		bnds_i=[0,-1]
	else:
		bnds_i=bound_finder(xarr,bnds)

	arr=arr[bnds_i[0]:bnds_i[1]]
	max_i=np.where(arr==np.amax(arr))[0][0]+bnds_i[0]
	return np.amax(arr),xarr[max_i],max_i

#------------------------------------------------------------------------------------------
def rebin(array,binsize):
	
	import numpy as np
	
	return array.reshape(np.shape(array)[1]//binsize[1],binsize[1],np.shape(array)[0]//binsize[0],binsize[0]).sum(3).sum(1)

#------------------------------------------------------------------------------------------
def cup_sig(cen,wx,wy,amp,cup_r):
	
	import numpy as np
	import funk
	
	x,y=np.linspace(-5,5,100),np.linspace(-5,5,100)
	sig=[]
	my,mx=np.meshgrid(x,y)
	for i in cen:
		pulse=funk.gauss2d_fwhm(x,y,amp,i,0.,wx,wy)
		mask=mx**2+my**2<=cup_r**2
		sig.append(np.sum(mask*pulse))
	return np.array(sig)
 
 #------------------------------------------------------------------------------------------
def rnd(num,dec):
	 
	numstr=str(num)
	for i in range(len(numstr)):
		if numstr[i]=='.':
			rnd_num=numstr[:i+1+dec]
			break
	return rnd_num

#-----------------------------------------------------------------------------------------
# Numerical differentiator

def num_diff(xdata,ydata,binsz,**kwargs):
	
	plot=kwargs.get('plot',True)

	import numpy as np
	import matplotlib.pyplot as plt

	# Do a moving average
	dy=np.zeros(np.shape(ydata))
	for i in range(np.shape(ydata)[0]):
		if i<binsz:
			dy[i]=(ydata[i+binsz]-ydata[0])/(xdata[i+binsz]-xdata[0])
		if i>=np.shape(ydata)[0]-(binsz+1):
			dy[i]=(ydata[-1]-ydata[i-binsz])/(xdata[-1]-xdata[i-binsz])
		else:
			dy[i]=(ydata[i+binsz+1]-ydata[i-binsz])/(xdata[i+binsz+1]-xdata[i-binsz])
			
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

#------------------------------------------------------------------------------------------

def time(ts):
	
	hms=ts.split(':')
	return float(hms[0])*3600+float(hms[1])*60+float(hms[2])

#------------------------------------------------------------------------------------------

def twin_plot(data1,data2,**kwargs):
	
	labels=kwargs.get('axlabels',['','',''])
	lims=kwargs.get('axlims','auto')
	cols=kwargs.get('colors',['blue','red'])
	
	import matplotlib.pyplot as plt
	
	fig=plt.figure('twin')
	ax=fig.add_subplot(111)
	ax.plot(data1[0],data1[1],color=cols[0])
	ax2 = ax.twinx()
	ax2.plot(data2[0],data2[1],color=cols[1])
	ax.set_xlabel(labels[0])
	ax.set_ylabel(labels[1])
	ax2.set_ylabel(labels[2])
	if lims!='auto':
		ax.set_xlim(lims[0])
		ax2.set_ylim(lims[2])
		ax.set_ylim(lims[1])
	plt.show()

#-------------------------------------------------------------------------------------------

def save_i(data,savepath):
	
	import numpy as np
	
	array=np.zeros((2,len(data)))
	array[0,:]=np.real(data)
	array[1,:]=np.imag(data)
	
	save(array,['Re','Im'],savepath)
