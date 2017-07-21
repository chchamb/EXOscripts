# Script to read random crap

def reader(filename,**kwargs):
	
	header=kwargs.get('header',True)
	d_type=kwargs.get('dtype','float')
	delim=kwargs.get('delimeter',' ')
	
	import numpy as np
	
	datafile=open(filename,'r')
	
	columns=[]
	for line in datafile:
		line=line.strip()
		columns.append(line.split(delim))
	datafile.close()
	
	if header==True:
		matrix=np.array(columns[1:],dtype=d_type)
		return columns[0],matrix
	if header==False:
		matrix=np.array(columns,dtype=d_type)
		return matrix
		
#-----------------------------------------------------------------------------------------

def reader_list(filename,**kwargs):
	
	header=kwargs.get('header',False)
	d_type=kwargs.get('dtype','float')
	delim=kwargs.get('delimeter',' ')
	
	import numpy as np
	
	datafile=open(filename,'r')
	
	columns=[]
	for line in datafile:
		line=line.strip()
		columns.append(line.split(delim))
	datafile.close()
	
	if header==True:
		return columns[0],columns
	if header==False:
		return columns	
		
#-----------------------------------------------------------------------------------------	
def bound_read(filename):
	import numpy as np
	
	datafile=open(filename,'r')
	
	lbounds=[]
	ubounds=[]
	for line in datafile:
		line=line.strip()
		columns=line.split()
		lbounds.append(columns[0])
		ubounds.append(columns[1])
	
	data=map(list,zip(lbounds[1:],ubounds[1:]))
	
	return data

def matrix_reader(filename,**kwargs):

	l=kwargs.get('skip_line',0)
	
	import numpy as np
	
	datafile=open(filename,'r')
	
	rows=[]
	i=0
	for line in datafile:
		line=line.strip()
		if i>(l-1):
			rows.append(line.split())
		i += 1
	datafile.close()
	matrix=np.transpose(np.array(rows))
	
	return matrix
	
#-----------------------------------------------------------------------------------------
# Script to read waveforms

def wf_reader(filename,**kwargs):
	
	header=kwargs.get('header',False)
	delim=kwargs.get('delimeter',None)
	
	import numpy as np
	import utility as ut
	
	datafile=open(filename,'r')
	
	time,columns=[],[]
	for line in datafile:
		line=line.strip()
		if len(line.split(delim))>=2:
			time.append(line.split(delim)[1])
			columns.append(line.split(delim)[2:])
	datafile.close()
	
	if header==True:
		matrix=np.array(columns[1:],dtype='float')
		matrix=np.hstack((ut.timestamp(time[1:],t0=True),matrix))
		return columns[0],matrix
	if header==False:
		matrix=np.array(columns,dtype='float')
		matrix=np.hstack((ut.timestamp(time,t0=True),matrix))
		return matrix
		
#------------------------------------------------------------------------------------------
