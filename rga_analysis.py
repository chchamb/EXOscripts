#=========================================================================================
# Scripts for RGA analysis
#=========================================================================================
date_g=20170426

# Script to parse and plot RGA data

def rga_read(filepath,**kwargs):
	
	date2=kwargs.get('date','same')
	
	if date2 != 'same':
		date=date2
	if date2 == 'same':
		date=date_g
	
	# Import Modules
	import numpy as np
		
	# Import and parse data (remove the top part of header)
	datafile=open('/home/chris/anaconda/data/'+str(date)+'/rga/pvt'+str(filepath)+'.txt','r')
	species=[]
	p_data=[]
	chan=0
	str_line=0
	for line in datafile:
		line=line.strip()
		columns=line.split()
		# This grabs the name of the channel in the header (Nitrogen...) and counts the number of channels
		if len(columns)==6:
			species.append(columns[2])
			chan +=1
			continue
		# This puts the data (except the first line) into a list of lists
		if len(columns)==(chan+1):
			if str_line>0:
				for i in range(len(columns)):
					a=columns[i]
					columns[i]=float(a[:-1])
				p_data.append(columns)
			str_line+=1
	# Put the data into an array
	p_data=np.array(p_data)
	
	return species,p_data

#-----------------------------------------------------------------------------------------

# Script to plot RGA data

def rga_plot(filepath,**kwargs):

	# Define keyword arguments
	chan2plot=kwargs.get('chan2plot','all')
	bounds=kwargs.get('bounds','full')
	yscale=kwargs.get('yscale','linear')
	tag=kwargs.get('tag','')
	o=kwargs.get('overlay',False)
	col=kwargs.get('colors','default')
	style=kwargs.get('style','-')
	date2=kwargs.get('date','same')
	
	if date2 != 'same':
		date=date2
	if date2 == 'same':
		date=date_g
	
	# Import modules
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	
	# Get the data from the file using rga_read
	species,p_data=rga_read(filepath,date=date2)
	species,p_data=rga_read(filepath,date=date2)
	
	# Separate time and pressure data
	time=p_data[:,0]
	p_data=p_data[:,1:]
	dt=ut.avediff(time)
	
	# Deal with plotting bounds
	if bounds=='full':
		bounds_s=[np.amin(time),np.amax(time)]
		bound_i=[0,len(time)-1]
		
	# Convert time to seconds
	if bounds!='full':
		bounds_s=ut.timestamp(bounds)
		# Find the index of the time closest to bounds
		bound_index=[]
		for k in range(2):
			j=-1
			for i in time[1:]:
				j=j+1
				if i<bounds_s[k]:
					last_time=i
				else:
					a=abs(bounds_s[k]-last_time)
					b=abs(i-bounds_s[k])
					if a>b:
						bound_index.append(j)
						break
					else:
						bound_index.append(j-1)
						break
		bound_i=ut.bound_finder(time,bounds_s)
				
	# Deal with which channels to intt
	if chan2plot=='all':
		chan2plot=species
	if chan2plot=='AllXenon':
		chan2plot=['Xenon129','Xenon131','Xenon132','Xenon134']
	
	# Plot the lines
	plt.figure('rga_plot')
	plt.clf()
	j=0
	figtext=''
	for i in chan2plot:
		if i=='Nitrogen':
			col='FireBrick'
		if i=='Oxygen':
			col='Turquoise'
		if i=='Hydrogen':
			col='Goldenrod'
		if i=='CarbonDioxide':
			col='DarkGreen'
		if i=='Water':
			col='MidnightBlue'
		if i=='Xenon129':
			col='Violet'
		if i=='Xenon131':
			col='DarkViolet'
		if i=='Xenon132':
			col='DarkMagenta'
		if i=='Xenon134':
			col='DarkSlateBlue'
		if i=='Argon':
			col='DeepPink'
		if i=='mass83':
			col='DodgerBlue'
		if i=='Barium138':
			col='DodgerBlue'		
		plt.plot(time[bound_i[0]:bound_i[1]],p_data[bound_i[0]:bound_i[1],species.index(i)],color=col,label=i)
	plt.xlim(time[bound_i[0]],time[bound_i[1]])
	if yscale=='log':
		ut.log()
	if yscale=='linear':
		ut.sci()
	plt.xlabel('Time (s)')
	plt.ylabel('Partial Pressure (Torr)')
	plt.title(str(date)+' pvt'+str(filepath))
	plt.legend(fontsize=12)
	plt.show()
	

#----------------------------------------------------------------------------------------
# Script to average RGA data

def rga_ave(filepath,**kwargs):
	
	# Define keyword arguments
	channels=kwargs.get('channels','all')
	bounds=kwargs.get('bounds',[])
	yscale=kwargs.get('yscale','linear')
	tag=kwargs.get('tag','')
	col=kwargs.get('colors','default')
	style=kwargs.get('style','-')
	
	# Import modules
	import numpy as np
	import matplotlib.pyplot as plt
	
	# Get run number from data file name
	run=' '
	j=1                                                        # j is the incrementor for characters in the run name/number
	for i in range(1,len(filepath)+1):
		if filepath[-i]=='/':                                      # Steps backwards in path until it finds a '/'
			while filepath[-i+j]!='_':
				run += filepath[-i+j]   						   # Adds on characters after the '/' until it finds a '_'
				j+=1
			break
		else:
			continue
	
	# Get the data from the file using rga_read
	species,p_data,chan=rga_read(filepath)
	
	# Separate time and pressure data
	time=p_data[:,0]
	p_data=p_data[:,1:]
	dt=time[1]-time[0]
	
	# Convert time to seconds
	bounds=np.array(bounds)
	bounds_s=np.zeros(np.shape(bounds))
	for i in range(np.shape(bounds)[0]):
		for j in range(np.shape(bounds)[1]):
			ftr=[3600,60,1]
			bounds_s[i,j]=sum([a*b for a,b in zip(ftr, map(int,bounds[i,j].split(':')))])
	
	# Find the index of the time closest to bounds
	bound_index=np.zeros(np.shape(bounds))
	for l in range(np.shape(bounds)[0]):
		for k in range(2):
			j=-1
			for i in time[1:]:
				j=j+1
				if i<bounds_s[l,k]:
					last_time=i
				else:
					a=abs(bounds_s[l,k]-last_time)
					b=abs(i-bounds_s[l,k])
					if a>b:
						bound_index[l,k]=j
						break
					else:
						bound_index[l,k]=j-1
						break
						
	# Make array with just data from within bounds (channel,bounds,values)
	num_points=bound_index[0,1]-bound_index[0,0]
	p_bounds=np.zeros((len(channels),np.shape(bound_index)[0],num_points))
	for i in range(np.shape(p_bounds)[0]):
		for j in range(np.shape(p_bounds)[1]):
			p_bounds[i,j,:]=p_data[bound_index[j,0]:bound_index[j,0]+num_points,species.index(channels[i])]
	# Average across the bounds
	p_aves=np.mean(p_bounds,axis=1)
	
	# Plot the averages on current axis
	for i in range(np.shape(p_aves)[0]):
		plt.plot(np.arange(0,dt*np.shape(p_aves)[1],dt),p_aves[i,:],style,color=col[i],label=channels[i]+': '+tag)
		plt.xlabel('Time (s)')
		plt.ylabel('Average Partial Pressure (Torr)')
		ax=plt.gca()
		ax.set_yscale(yscale)
		if yscale=='linear':
			plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
		plt.legend()
		plt.draw()
	
	return p_aves
	
#----------------------------------------------------------------------------------------

# Script to plot the ratio of peaks due to impurities during leaks

def rga_peak(filepath,imp,peaks):
	
	# Import modules
	import numpy as np
	
	# Get run number from data file name
	run=' '
	j=1
	for i in range(1,len(filepath)+1):
		if filepath[-i]=='/':
			while filepath[-i+j]!='_':
				run=run+filepath[-i+j]
				j=j+1
			break
		else:
			continue
	
	# Get the data from the file using rga_read
	species,p_data,chan=rga_read(filepath)
	
	# Separate time and pressure data
	time=p_data[:,0]
	p_data=p_data[:,1:]
	
	# Convert time to seconds
	peaks_s=[]
	for i in peaks:
		ftr=[3600,60,1]
		peaks_s.append(sum([a*b for a,b in zip(ftr, map(int,i.split(':')))]))
	
	# Find the index of the time closest to xenon peak
	peak_index=[]
	for k in range(len(peaks_s)):
		j=-1
		for i in time[1:]:
			j=j+1
			if i<peaks_s[k]:
				last_time=i
			else:
				a=abs(peaks_s[k]-last_time)
				b=abs(i-peaks_s[k])
				if a>b:
					peak_index.append(j)
					break
				else:
					peak_index.append(j-1)
					break
	
	# Get the maximum in region around peak subtracting the baseline
	xe_max=[]
	imp_max=[]
	for i in range(len(peaks_s)):
		xe_sub_max=[]
		for k in ['Xenon129','Xenon131','Xenon132']:
			xe_sub_max.append(np.amax(p_data[peak_index[i]-10:peak_index[i]+20,species.index(k)]-np.mean(p_data[peak_index[i]-25:peak_index[i]-5,species.index(k)])))
		xe_max.append(np.sum(xe_sub_max))
		imp_sub_max=[]
		for j in imp:
			imp_sub_max.append(np.amax(p_data[peak_index[i]-10:peak_index[i]+20,species.index(j)]-np.mean(p_data[peak_index[i]-25:peak_index[i]-5,species.index(j)])))
		imp_max.append(imp_sub_max)
	
	# Make into numpy arrays
	xe_max=np.array(xe_max)
	imp_max=np.array(imp_max)
	
	# Make array with ratios
	ratio=[]
	for i in range(len(peaks)):
		ratio_leak=[]
		for j in range(len(imp)):
			ratio_leak.append(imp_max[i,j]/xe_max[i])
		ratio.append(ratio_leak)
	ratio=np.array(ratio)
	
	return imp,ratio,imp_max,xe_max

#----------------------------------------------------------------------------------------								
								
# Script to plot the ratio of peaks due to impurities during leaks

def rga_int(filepath,imp,bounds):
	
	# Import modules
	import numpy as np
	
	# Get run number from data file name
	run=' '
	j=1
	for i in range(1,len(filepath)+1):
		if filepath[-i]=='/':
			while filepath[-i+j]!='_':
				run=run+filepath[-i+j]
				j=j+1
			break
		else:
			continue
	
	# Get the data from the file using rga_read
	species,p_data,chan=rga_read(filepath)
	
	# Separate time and pressure data
	time=p_data[:,0]
	p_data=p_data[:,1:]
	bounds=np.array(bounds)
	
	# Convert time to seconds
	bounds_s=np.zeros(np.shape(bounds))
	for i in range(np.shape(bounds)[0]):
		for j in range(2):
			ftr=[3600,60,1]
			bounds_s[i,j]=sum([a*b for a,b in zip(ftr, map(int,bounds[i,j].split(':')))])

	# Find the index of the time closest to bounds
	bound_index=np.zeros(np.shape(bounds))
	for l in range(np.shape(bounds)[0]): 	
		for k in range(2):
			j=-1
			for i in time[1:]:
				j=j+1
				if i<bounds_s[l,k]:
					last_time=i
				else:
					a=abs(bounds_s[l,k]-last_time)
					b=abs(i-bounds_s[l,k])
					if a>b:
						bound_index[l,k]=j
						break
					else:
						bound_index[l,k]=j-1
						break
	
	# Get the integral in region around peak subtracting the baseline
	xe_int=np.zeros(np.shape(bounds)[0])
	imp_int=np.zeros((len(imp),np.shape(bounds)[0]))
	for i in range(np.shape(bounds_s)[0]):
		l=0
		for k in ['Xenon129','Xenon131','Xenon132']:
			xe_int[i] += np.sum(p_data[bound_index[i,0]:bound_index[i,1],species.index(k)]-np.mean(p_data[bound_index[i,0]-25:bound_index[i,0]-5,species.index(k)]))
		for j in imp:
			imp_int[l,i]=np.sum(p_data[bound_index[i,0]:bound_index[i,1],species.index(j)]-np.mean(p_data[bound_index[i,0]-25:bound_index[i,0]-5,species.index(j)]))
			l +=1
			
	# Make array with ratios
	ratio=np.zeros(np.shape(imp_int))
	for i in range(np.shape(bounds)[0]):
		for j in range(len(imp)):
			ratio[j,i]=imp_int[j,i]/xe_int[i]
		
	return imp,ratio,imp_int,xe_int

#------------------------------------------------------------------------------------------
# Script to integrate RGA data

def rga_int(filepath,bounds,**kwargs):

	# Define keyword arguments
	chan2int=kwargs.get('chan2int','all')
	yscale=kwargs.get('yscale','log')
	col=kwargs.get('colors','default')
	output=kwargs.get('store',False)
	pltsub=kwargs.get('plt_sub',False)
	bndtype=kwargs.get('boundtype','edge')
	
	# Import modules
	import numpy as np
	import matplotlib.pyplot as plt
	import utility as ut
	
	# Get the data from the file using rga_read
	species,p_data=rga_read(filepath)
	
	# Separate time and pressure data
	time=p_data[:,0]
	p_data=p_data[:,1:]
	dt=ut.avediff(time)
	
	# Convert time to seconds
	if type(bounds[0])==str:
		bounds=ut.timestamp(bounds)
	
	if bndtype=='front':
		bounds=[bounds[0],bounds[0]+bounds[1]]
	
	# Find the index of the time closest to bounds
	bound_i=ut.bound_finder(time,bounds)
				
	# Deal with which channels to integrate
	if chan2int=='all':
		chan2int=species
	if chan2int=='AllXenon':
		chan2int=['Xenon129','Xenon131','Xenon132','Xenon134']
	
	# Integrate selected channels over bounds
	ints,bg=[],[]
	for i in chan2int:
		bg.append(np.mean(p_data[bound_i[0]-20:bound_i[0],species.index(i)]))
		ints.append(np.sum(p_data[bound_i[0]:bound_i[1],species.index(i)])-bg[-1])
		
	# Plot the lines
	plt.figure('rga_int')
	plt.clf()
	j=0
	figtext=''
	for i in chan2int:
		if i=='Nitrogen':
			col='FireBrick'
		if i=='Oxygen':
			col='Turquoise'
		if i=='Hydrogen':
			col='Goldenrod'
		if i=='CarbonDioxide':
			col='DarkGreen'
		if i=='Water':
			col='MidnightBlue'
		if i=='Xenon129':
			col='Violet'
		if i=='Xenon131':
			col='DarkViolet'
		if i=='Xenon132':
			col='DarkMagenta'
		if i=='Xenon134':
			col='DarkSlateBlue'
		if i=='Argon':
			col='DeepPink'
		if i=='mass83':
			col='DodgerBlue'	
		if pltsub==True:
			plt.plot(time[bound_i[0]-10:bound_i[1]+10],p_data[bound_i[0]-10:bound_i[1]+10,species.index(i)]-bg[j],color=col,label=i)	
		if pltsub==False:
			plt.plot(time[bound_i[0]-10:bound_i[1]+10],p_data[bound_i[0]-10:bound_i[1]+10,species.index(i)],color=col,label=i)
		if j==0:
			figtext+=i+': '+ut.conv(ints[j])+'Torr'
		if j>0:
			figtext+='\n'+i+': '+ut.conv(ints[j])+'Torr'
		j+=1
	plt.xlim(time[bound_i[0]-5],time[bound_i[1]+5])
	plt.axvline(bounds[0],color='black',ls='--')
	plt.axvline(bounds[1],color='black',ls='--')
	if yscale=='log':
		ut.log()
	if yscale=='linear':
		ut.sci()
	plt.xlabel('Time (s)')
	plt.ylabel('Partial Pressure (Torr)')
	plt.title(str(date)+' RGA integral')
	ut.textbox(figtext,[.05,.95])
	plt.legend()
	plt.show()
	
	if output==True:
		return ints
	
