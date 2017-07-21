# Script containing common (but not too common) functions
import numpy as np
#-----------------------------------------------------------------------------------------
def gauss(x,amp,cen,width,wtype):
	if wtype=='sigma':
		width=width
	if wtype=='fwhm':
		width=width/(2.*np.sqrt(2*np.log(2)))
	if wtype=='w':
		width=width/2.
	return amp*np.exp(-np.square(x-cen)/(2*np.square(width)))
#-----------------------------------------------------------------------------------------
def gauss_w(x,amp,cen,width,offset):
	width=width/2.
	return amp*np.exp(-np.square(x-cen)/(2*np.square(width)))+offset
#-----------------------------------------------------------------------------------------
def gauss_fwhm(x,amp,cen,width,offset):
	width=width/(2.*np.sqrt(2*np.log(2)))
	return amp*np.exp(-np.square(x-cen)/(2*np.square(width)))+offset
#-----------------------------------------------------------------------------------------
def m_gauss(x,amp1,amp2,cen1,cen2,width1,width2):
	return gauss_w(x,amp1,cen1,width1,0)+gauss_w(x,amp2,cen2,width2,0)
#-----------------------------------------------------------------------------------------
def lorentz(x,cen,fwhm,offset):
	return (2/(np.pi*fwhm))*(1./(1+np.square((x-cen)/(fwhm/2))))+offset
#-----------------------------------------------------------------------------------------
def asym(x,amp,cen,leftw,rightw,offset):
	import scipy.special as sp
	return amp*(1.-sp.erf((x-cen)/rightw))*(1.+sp.erf((x-cen)/leftw))+offset
#-----------------------------------------------------------------------------------------
def gauss2d(x,y,amp,xcenter,ycenter,width,wtype):  # This ends up as f(y,x)
	my,mx=np.meshgrid(x,y)
	if wtype=='sigma':
		width=width
	if wtype=='fwhm':
		width=width/(2.*np.sqrt(2*np.log(2)))
	if wtype=='w':
		width=width/2.	
	return amp*np.exp(-(np.square(mx-ycenter)+np.square(my-xcenter))/(2*np.square(width)))
	
#-----------------------------------------------------------------------------------------
def gauss2dw(x,y,amp,xcenter,ycenter,xwidth,ywidth):  # This ends up as f(y,x)
	my,mx=np.meshgrid(x,y)
	xwidth,ywidth=xwidth/2.,ywidth/2.	
	return amp*np.exp(-(np.square(mx-ycenter)/(2*np.square(ywidth))))*np.exp(-(np.square(my-xcenter))/(2*np.square(xwidth)))
#-----------------------------------------------------------------------------------------
def gauss2d_fwhm(x,y,amp,xcenter,ycenter,xwidth,ywidth):  # This ends up as f(y,x)
	my,mx=np.meshgrid(x,y)
	xwidth,ywidth=xwidth/(2.*np.sqrt(2*np.log(2))),ywidth/(2.*np.sqrt(2*np.log(2)))	
	return amp*np.exp(-(np.square(mx-ycenter)/(2*np.square(ywidth))))*np.exp(-(np.square(my-xcenter))/(2*np.square(xwidth)))
#-----------------------------------------------------------------------------------------
def decay(x,*args):      # if only 2 args:(amp,tau), if 3 args:(amp,tau,offset)
	if len(args)==3:
		return args[0]*np.exp(-x/args[1])+args[2]
	else:		
		return args[0]*np.exp(-x/args[1])
#-----------------------------------------------------------------------------------------
def beam(z,w_0,f,wlength):
	z_0=(np.pi*np.square(w_0))/(wlength*1E-6)
	return w_0*np.sqrt(1+np.square((z-f)/z_0))
#------------------------------------------------------------------------------------------
def beamxy(z,w_0,n,wlength):
	wlength=float(wlength*1E-3)  # convert to micron
	theta=wlength/(np.pi*w_0)
	alpha=(np.arcsin((1./n)*np.sin((np.pi/4)+theta))-np.arcsin((1./n)*np.sin((np.pi/4)-theta)))/2
	w_0n=np.array([w_0,wlength/(np.pi*n*alpha)])
	z_0=(np.pi*np.square(w_0n)*n)/(wlength)
	w_x=w_0n[0]*np.sqrt(1+np.square((n*z*1.074)/z_0[0]))
	w_y=w_0n[1]*np.sqrt(1+np.square((3.*z)/z_0[1]))
	return 'w_x: '+str(w_x),'w_y: '+str(w_y)
#------------------------------------------------------------------------------------------
def poly(x,*args):
	func=0
	for i in range(len(args)):
		func+=args[i]*np.power(x,len(args)-(i+1))
	return func
#------------------------------------------------------------------------------------------
def line0(x,m):
	return m*x
#------------------------------------------------------------------------------------------
def sin(x,amp,frq,phi,offset):
	return amp*np.sin((x*frq)-phi)+offset
#------------------------------------------------------------------------------------------
def gaussarray(shape,spacing,xwidth,ywidth,**kwargs):
	
	import matplotlib.pyplot as plt
	
	start=kwargs.get('start',False)
	# Create x and y arrays
	x=np.linspace(-(shape[0]-1)*spacing/1.5,(shape[0]-1)*spacing/1.5,100)
	y=np.linspace(-(shape[1]-1)*spacing/1.5,(shape[1]-1)*spacing/1.5,100)

	# Create array of centers
	centers=[]
	for i in range(shape[0]):
		centers2=[]
		for j in range(shape[1]):
			centers2.append([i*spacing,j*spacing])
		centers.append(centers2)
	
	# Shift center of array to (0,0)
	centers=np.array(centers)-[((shape[0]-1)*spacing)/2.,((shape[1]-1)*spacing)/2.]
	
	# Make the gaussian array with the centers
	garray=0
	for i in range(np.shape(centers)[0]):
		for j in range(np.shape(centers)[1]):
			garray+=gauss2dw(x,y,1.,centers[i,j,0],centers[i,j,1],xwidth,ywidth)
	
	if start==True:
		return garray,centers[0,0,:]
	else:
		return garray

#------------------------------------------------------------------------------------------
def gtrain(x,amp,x0,dx,w,N):
	
	func=0
	for i in range(N):
		func+=gauss_w(x,amp,x0+i*dx,w)
	return func

#------------------------------------------------------------------------------------------
def exp_on(x,amp,tau,offset):
	return amp*(1-np.exp(-x/tau))+offset
#------------------------------------------------------------------------------------------
def gauss2d_flat(x,y,amp,xcenter,ycenter,xwidth,ywidth):  # This ends up as f(y,x)
	my,mx=np.meshgrid(x,y)
	xwidth,ywidth=xwidth/2.,ywidth/2.	
	return np.ravel(amp*np.exp(-2*(np.square(mx-ycenter)/(np.square(ywidth))))*np.exp(-2*(np.square(my-xcenter))/(np.square(xwidth))))
#------------------------------------------------------------------------------------------
def root(x,x0,b,m):
	return np.sqrt(m*x-x0)+b

#------------------------------------------------------------------------------------------
def erf(x,a,b):
	import scipy.special as sp
	return a*sp.erf(x*b)
#------------------------------------------------------------------------------------------
def rmse(x):
	return np.sqrt(np.mean(np.square(x-np.mean(x))))
