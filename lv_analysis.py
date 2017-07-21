# Script to read and process labview data
date=20170720
#------------------------------------------------------------------------------------------
def lv_ave(runtable,**kwargs):
	
	savefile=kwargs.get('savefile',False)
	
	import matplotlib.pyplot as plt
	import numpy as np
	import asciitable as asc
	import gen_reader as gr

	datafile=open(runtable,'r',)
	
	runnum,lvnum=[],[]
	for line in datafile:
		line=line.strip()
		column=line.split(',')
		runnum.append(float(column[0]))
		lvnum.append(column[1])
	datafile.close()
	runnum=np.array(runnum)
	
	# Do the average on the lv files
	avpow=[]
	for i in lvnum:
		lv_matrix=gr.reader('/home/chris/anaconda/data/'+str(date)+'/lv/'+str(i)+'_laser_temp.txt',header=False,delimeter=';')
		avpow.append(float(np.mean(lv_matrix[:,1])))
	avpow=np.array(avpow)
	runtable=np.vstack((runnum,avpow))
	if savefile!=False:
		plt.figure('lv_ave')
		plt.clf()
		plt.plot(runnum,1000*avpow,'bo')
		plt.xlabel('Run')
		plt.ylabel('Laser Power (mW)')
		plt.savefig('/home/chris/anaconda/plots/'+str(date)+'/'+savefile+'.png')
	return runtable
	
#------------------------------------------------------------------------------------------

def get_avpow(lvnum,**kwargs):
	
	T_R=kwargs.get('T_R',1.)
	op_flat=kwargs.get('op_flat',True)
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	
	if op_flat==True:
		flat_T=.89
	if op_flat==False:
		flat_T=1.
	
	if len(str(lvnum))==1:
		lvnum='000'+str(lvnum)
	if len(str(lvnum))==2:
		lvnum='00'+str(lvnum)
	if len(str(lvnum))==3:
		lvnum='0'+str(lvnum)
	lv_matrix=gr.reader('/home/chris/anaconda/data/'+str(date)+'/lv/'+str(lvnum)+'_laser_temp.txt',header=False,delimeter=';')
	avpow=np.mean(lv_matrix[:,1])*T_R*flat_T #1.57 is power meter compensation
	
	return avpow
	
#------------------------------------------------------------------------------------------

def lv_energy(lvnum,**kwargs):
	
	T_R=kwargs.get('T_R',1.)
	op_flat=kwargs.get('op_flat',True)
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	
	if op_flat==True:
		flat_T=.89
	if op_flat==False:
		flat_T=1.
	
	if len(str(lvnum))==1:
		lvnum='000'+str(lvnum)
	if len(str(lvnum))==2:
		lvnum='00'+str(lvnum)
	if len(str(lvnum))==3:
		lvnum='0'+str(lvnum)
	lv_matrix=gr.reader('/home/chris/anaconda/data/'+str(date)+'/lv/'+str(lvnum)+'_laser_temp.txt',header=False,delimeter=';')
	
	return lv_matrix[:,1]*T_R*flat_T
#------------------------------------------------------------------------------------------

def get_avtime(lvnum):
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	
	if len(str(lvnum))==1:
		lvnum='000'+str(lvnum)
	if len(str(lvnum))==2:
		lvnum='00'+str(lvnum)
	if len(str(lvnum))==3:
		lvnum='0'+str(lvnum)
	lv_matrix=gr.reader('/home/chris/anaconda/data/'+str(date)+'/lv/'+str(lvnum)+'_laser_temp.txt',header=False,delimeter=';')
	avpow=np.mean(lv_matrix[:,-1])
	
	return avpow
#------------------------------------------------------------------------------------------

def lv_time(lvnum):
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	
	if len(str(lvnum))==1:
		lvnum='000'+str(lvnum)
	if len(str(lvnum))==2:
		lvnum='00'+str(lvnum)
	if len(str(lvnum))==3:
		lvnum='0'+str(lvnum)
	lv_matrix=gr.reader('/home/chris/anaconda/data/'+str(date)+'/lv/'+str(lvnum)+'_laser_temp.txt',header=False,delimeter=';')
	avpow=lv_matrix[:,-1]
	
	return avpow
	
#------------------------------------------------------------------------------------------

def lv_energy_bad(lvnum,**kwargs):
	
	T_R=kwargs.get('T_R',1.)
	op_flat=kwargs.get('op_flat',True)
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	import glob
	
	lvnum=str(lvnum).zfill(4)
	f=len(glob.glob('/home/chris/anaconda/data/'+str(date)+'/lv/'+lvnum+'_*_laser_temp.txt'))
	
	if op_flat==True:
		flat_T=.89
	if op_flat==False:
		flat_T=1.
	
	powsum=np.zeros(f)
	for i in range(f):
		frnum=str(i+1).zfill(4)
		lv_matrix=gr.wf_reader('/home/chris/anaconda/data/'+str(date)+'/lv/'+lvnum+'_'+frnum+'_laser_temp.txt',header=False)
		powsum[i]=np.sum(lv_matrix[:,2])*T_R*flat_T
	
	dt=[]
	for i in range(np.shape(lv_matrix)[0]-1):
		dt.append(lv_matrix[i+1,0]-lv_matrix[i,0])
	dt=np.mean(dt)

	energy=dt*powsum
	
	return energy
	
#------------------------------------------------------------------------------------------

def lv_avpow(lvnum,**kwargs):
	
	T_R=kwargs.get('T_R',1.)
	op_flat=kwargs.get('op_flat',True)
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	import glob
	
	lvnum=str(lvnum).zfill(4)
	f=len(glob.glob('/home/chris/anaconda/data/'+str(date)+'/lv/'+lvnum+'_*_laser_temp.txt'))
	
	if op_flat==True:
		flat_T=.89
	if op_flat==False:
		flat_T=1.
	
	avpow=np.zeros(f)
	for i in range(f):
		frnum=str(i+1).zfill(4)
		lv_matrix=gr.wf_reader('/home/chris/anaconda/data/'+str(date)+'/lv/'+lvnum+'_'+frnum+'_laser_temp.txt',header=False)
		avpow[i]=np.mean(lv_matrix[:,2])*T_R*flat_T
	
	return avpow
#------------------------------------------------------------------------------------------

def trig_check(trignum,**kwargs):
	
	ret_full=kwargs.get('return_full',False)
	
	import numpy as np
	import matplotlib.pyplot as plt
	import gen_reader as gr
	
	trigdat=gr.wf_reader('/home/chris/anaconda/data/'+str(date)+'/lv/trig_'+str(trignum)+'.txt')
	
	# Calculate dt and find the accelerometer level during trigger window
	dt,acc_sig=[],[]
	for i in range(np.shape(trigdat)[0]):
		if i>0:
			dt.append(trigdat[i,0]-trigdat[i-1,0])
		if trigdat[i,2]==1:
			acc_sig_1.append(trigdat[i,1])
	
	acc_sig_ave=np.mean(acc_sig)
	dt_ave=np.mean(dt)
	
	plt.figure('trig')
	plt.clf()
	plt.plot(acc_sig,'b')
	plt.title(str(date)+' run'+str(trignum)+' Trigger Signal while Laser Shutter Open')
	
	if ret_full==True:
		return dt,acc_sig
	return dt_ave,acc_sig_ave
	
#------------------------------------------------------------------------------------------

def trig_check2(trignum,**kwargs):
	
	savepath=kwargs.get('savepath',False)
	
	import numpy as np
	import matplotlib.pyplot as plt
	import gen_reader as gr
	import utility as ut
	
	trigdat=gr.wf_reader('/home/chris/anaconda/data/'+str(date)+'/lv/trig_'+str(trignum)+'.txt')
	
	# Calculate dt and find the accelerometer level during trigger window
	acc_sig,acc_sig_1=[],[]
	j=0
	for i in range(np.shape(trigdat)[0]):
		if trigdat[i,2]==0:
			j=0
		if trigdat[i,2]==1:
			if j==0:
				acc_sig.append(np.mean(acc_sig_1))
			acc_sig_1.append(trigdat[i,1])
			j+=1
	
	plt.figure('trig2')
	plt.clf()
	plt.plot(acc_sig[1:],'bo',markersize=3)
	plt.title(str(date)+' run'+str(trignum)+' Average Trigger Signal while Laser Shutter Open')
	plt.xlabel('Trigger Number')
	plt.ylabel('Accelerometer Signal (V)')
	if savepath!=False:
		ut.create_dir(savepath)
		plt.savefig(savepath+'trig_'+str(trignum )+'_acc_sig.png')
	
	return acc_sig[1:]	
