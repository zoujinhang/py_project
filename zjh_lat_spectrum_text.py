from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import Data_analysis.file as myfile
from Fermi_tool.spectrum import Fit,Spectrum

def band(e,cube):
	k = cube[0]
	a = cube[1]
	b = cube[2]
	Ec = cube[3]
	#print(k,a,b,Ec)
	if a<b:
		return [np.nan]
	val = 10**Ec*(a-b)
	a1 = k*(e[e<val]/100)**a*np.exp(-(e[e<val]/10**Ec))
	a2 = k*np.exp(b-a)*(val/100)**(a-b)*(e[e>=val]/100)**b
	A = np.concatenate((a1,a2))
	return A

def prior(cube,ndim,nparams):
	cube[0] = 10000*cube[0]
	cube[1] = 10*cube[1]-5
	cube[2] = 10*cube[2]-8
	cube[3] = (8.7*cube[3]+1.3)

def band2(e,cube):
	
	k = cube[0]
	a = cube[2]
	b = cube[3]
	c = cube[1]
	E0 = cube[4]
	Ec = cube[5]
	#print(k,a,b,Ec)
	#print(c)
	val = Ec*(a-b)
	e0 = e[np.where(e<E0)[0]]
	e1 = e[np.where((e>=E0)&(e<val))[0]]
	e2 = e[np.where(e>=val)[0]]
	try:
		a0 = k * (E0 / 100) ** (a - c) * np.exp(-E0 / Ec) * (e0 / 100) ** c
		a1 = k * (e1 / 100) ** a * np.exp(-(e1 / Ec))
		a2 = k*np.exp(b-a)*(val/100)**(a-b)*(e2/100)**b
	except:
		print(c,a,b,val,E0,Ec)
	A = np.concatenate((a0,a1,a2))
	return A

def prior2(cube,ndim,nparams):
	cube[0] = 100000*cube[0]#k
	cube[1] = 10*cube[1]#c
	cube[2] = 30*cube[2]-20#a
	if cube[2]<0:
		cube[3] = (cube[2]+30)*cube[3]-30#b
	else:
		cube[3] = 30*cube[3]-30#b
	cube[5] = 10**(20*cube[5]-10)#Ec
	if cube[5]*(cube[2]-cube[3])<100:
		cube[4] = cube[5]*(cube[2]-cube[3]) * cube[4]#E0
	else:
		cube[4] = 100*cube[4]

def bb(e,K,kt):
	A = K*8.0525*e**2/(kt)**4/(np.exp(e/kt)-1)
	return A

def model_bb2(e,cube):
	K1 = cube[0]
	kt1 = cube[1]
	K2 = cube[2]
	kt2 = cube[3]
	if kt1<kt2:
		return [np.nan]
	return bb(e,K1,kt1)+bb(e,K2,kt2)

def prior_bb2(cube,ndim,nparams):
	cube[0] = 10000*cube[0]
	cube[1] = 10**(20*cube[1]-10)
	cube[2] = 10000*cube[2]
	cube[3] = 10**(20*cube[3]-10)
	
rsp_link_b0 = '/home/laojin/my_lat/spectrum/response_matrix/glg_cspec_b0_bn150330828_v00.rsp'
rsp_link_n2 = '/home/laojin/my_lat/spectrum/response_matrix/glg_cspec_n2_bn150330828_v00.rsp'
savedir = '/home/laojin/my_lat/spectrum/A_spectrum/'
if os.path.exists(savedir)==False:
	os.makedirs(savedir)


for i in range(10):
	y_sp,y_sp_err = myfile.readcol(savedir +'ZZ_spectrum_n2_'+str(i)+'.txt')
	y_sp_n2 = np.array(y_sp)
	y_sp_err_n2 = np.array(y_sp_err)
	
	spectrum_n2 = Spectrum(y_sp_n2,y_sp_err_n2,rsp_link_n2,effective_band=[10,800],spectrum_name = 'n2')
	
	
	y_sp,y_sp_err = myfile.readcol(savedir +'ZZ_spectrum_b0_'+str(i)+'.txt')
	y_sp_b0 = np.array(y_sp)
	y_sp_err_b0 = np.array(y_sp_err)
	
	spectrum_b0 = Spectrum(y_sp_b0,y_sp_err_b0,rsp_link_b0,effective_band=[400,20000],spectrum_name = 'b0')
	
	parameters = ['K','a','b','logEc']
	parameters2 = ['K','a0','a','b','E0','Ec']
	parameters_bb2 = ['K1','kt1','K2','kt2']
	
	#fit = Fit([spectrum_n2],band2,prior2,parameters2)
	fit = Fit([spectrum_n2,spectrum_b0],band,prior,parameters)
	#fit = Fit([spectrum_n2,spectrum_b0],model_bb2,prior_bb2,parameters_bb2)
	mul_dir = savedir+'A_n'+str(i)+'_out/'
	if os.path.exists(mul_dir) ==False:
		os.makedirs(mul_dir)
	A = fit.run(outputfiles_basename=mul_dir+'A_n'+str(i)+'_',resume = False, verbose = True)
	
	fig,ax = plt.subplots()
	fit.plot_model(A,n=2,ax = ax,reference = True)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim(10**-2,)
	ax.legend()
	fig.savefig(savedir + 'B_plot_model_'+str(i)+'.png')
	plt.close(fig)
	
	fig,ax = plt.subplots()
	fit.plot_data(A,n=2,ax = ax,reference = True)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim(10**-2,)
	ax.legend()
	fig.savefig(savedir + 'A_plot_data_'+str(i)+'.png')
	plt.close(fig)
	
	fit.plot_corner(A)
	plt.savefig(savedir + 'C_plot_corner_'+str(i)+'.png')
	plt.close()


