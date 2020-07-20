from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import Data_analysis.file as myfile
from Fermi_tool.spectrum import Fit,Spectrum,Prior
from matplotlib.gridspec import GridSpec
import h5py


def band(e,cube):
	k = cube[0]
	a = cube[1]
	b = cube[2]
	Ec = cube[3]
	#print(k,a,b,Ec)
	if a<b:
		return [np.nan]
	val = 10**Ec*(a-b)
	a1 = 10**k*(e[e<val]/100)**a*np.exp(-(e[e<val]/10**Ec))
	a2 = 10**k*np.exp(b-a)*(val/100)**(a-b)*(e[e>=val]/100)**b
	A = np.concatenate((a1,a2))
	return A

def prior(cube,ndim,nparams):
	#print('ndim',ndim)
	#print('nparams',nparams)
	cube[0] = (12*cube[0]-5)
	cube[1] = 8*cube[1]-3
	cube[2] = 7*cube[2]-5
	cube[3] = (11.3*cube[3]-1.3)

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
	B = np.exp(e/kt)-1
	A = K*8.0525*e**2/(kt)**4/(B)
	#print(A)
	return A

def model_bb2_pl(e,cube):
	K1 = cube[0]
	kt1 = cube[1]
	K2 = cube[2]
	kt2 = cube[3]
	K3 = cube[4]
	a = cube[5]
	#print('kt1:',10**kt1,'kt2:',10**kt2,'K1',10**K1,'K2',10**K2)
	if kt1<kt2:
		return [np.nan]
	return bb(e,10**K1,10**kt1)+bb(e,10**K2,10**kt2) + 10**K3*e**(-a)




def prior_bb2(cube,ndim,nparams):
	cube[0] = 10000*cube[0]
	cube[1] = 10**(20*cube[1]-10)
	cube[2] = 10000*cube[2]
	cube[3] = 10**(20*cube[3]-10)
	
#rsp_link_b0 = '/home/laojin/my_lat/spectrum/response_matrix/glg_cspec_b0_bn150330828_v00.rsp'
#rsp_link_n2 = '/home/laojin/my_lat/spectrum/response_matrix/glg_cspec_n2_bn150330828_v00.rsp'
#rsp_link_n1 = '/home/laojin/my_lat/spectrum/response_matrix/glg_cspec_n1_bn150330828_v00.rsp'

rsp_link_b0 = '/media/laojin/Elements/trigdata/2009/bn090809978/glg_cspec_b0_bn090809978_v01.rsp'
rsp_link_n3 = '/media/laojin/Elements/trigdata/2009/bn090809978/glg_cspec_n3_bn090809978_v01.rsp'
rsp_link_n4 = '/media/laojin/Elements/trigdata/2009/bn090809978/glg_cspec_n4_bn090809978_v01.rsp'



datadir = '/home/laojin/my_lat/spectrum/bn090809978/'
#datadir = '/home/laojin/my_lat/spectrum/A_spectrum1/'
#savedir = '/home/laojin/my_lat/spectrum/A_spectrum1/'
savedir = '/home/laojin/my_lat/spectrum/A_spectrum_bn090809978/'
if os.path.exists(savedir)==False:
	os.makedirs(savedir)

band_prior_list = [
	Prior([-5,4]),
	Prior([-3,8]),
	Prior([-5,2]),
	Prior([-1.3,10])
]


model_bb2_pl_prior_list = [
	Prior([-2,10]),
	Prior([1,5]),
	Prior([-2,10]),
	Prior([1,5]),
	Prior([-5,5]),
	Prior([-3,5])
]

#n2 = h5py.File(datadir +'ZZ_spectrum_n2.hdf5','r')
n2 = h5py.File(datadir +'ZZ_spectrum_n4.hdf5','r')
n2_spc = n2['spc'][()]
n2_spc_err = n2['spc_err'][()]
n2.close()

b0 = h5py.File(datadir +'ZZ_spectrum_b0.hdf5','r')
b0_spc = b0['spc'][()]
b0_spc_err = b0['spc_err'][()]
b0.close()

#n1 = h5py.File(datadir +'ZZ_spectrum_n1.hdf5','r')
n1 = h5py.File(datadir +'ZZ_spectrum_n3.hdf5','r')
n1_spc = n1['spc'][()]
n1_spc_err = n1['spc_err'][()]
n1.close()

Ep_list = []
E0_list = []

kt1 = []
kt2 = []
a_list = []
for i in range(n2_spc.shape[0]):
#for i in range(567):
	#y_sp,y_sp_err = myfile.readcol(savedir +'ZZ_spectrum_n2_'+str(i)+'.txt')
	#y_sp_n2 = np.array(y_sp)
	#y_sp_err_n2 = np.array(y_sp_err)
	y_sp_n2 = n2_spc[i]
	y_sp_err_n2 = n2_spc_err[i]
	
	
	#spectrum_n2 = Spectrum(y_sp_n2,y_sp_err_n2,rsp_link_n2,effective_band=[8,800],spectrum_name = 'n2')
	spectrum_n2 = Spectrum(y_sp_n2,y_sp_err_n2,rsp_link_n4,effective_band=[8,800],spectrum_name = 'n4')
	
	#y_sp,y_sp_err = myfile.readcol(savedir +'ZZ_spectrum_b0_'+str(i)+'.txt')
	#y_sp_b0 = np.array(y_sp)
	#y_sp_err_b0 = np.array(y_sp_err)
	y_sp_b0 = b0_spc[i]
	y_sp_err_b0 = b0_spc_err[i]
	
	
	spectrum_b0 = Spectrum(y_sp_b0,y_sp_err_b0,rsp_link_b0,effective_band=[600,30000],spectrum_name = 'b0')
	
	#y_sp,y_sp_err = myfile.readcol(savedir +'ZZ_spectrum_n1_'+str(i)+'.txt')
	#y_sp_n1 = np.array(y_sp)
	#y_sp_err_n1 = np.array(y_sp_err)
	y_sp_n1 = n1_spc[i]
	y_sp_err_n1 = n1_spc_err[i]
	
	#spectrum_n1 = Spectrum(y_sp_n1,y_sp_err_n1,rsp_link_n1,effective_band=[8,800],spectrum_name = 'n1')
	spectrum_n1 = Spectrum(y_sp_n1,y_sp_err_n1,rsp_link_n3,effective_band=[8,800],spectrum_name = 'n3')
	parameters = ['log K','a','b','log Ec']
	parameters2 = ['K','a0','a','b','E0','Ec']
	parameters_bb2 = ['log K1','log kt1','log K2','log kt2']
	parameters_bb2_pl = ['log K1','log kt1','log K2','log kt2','log K3','a']
	#fit = Fit([spectrum_n2],band2,prior2,parameters2)
	fit = Fit([spectrum_n2,spectrum_n1,spectrum_b0],band,band_prior_list,parameters,reference=False)
	#fit = Fit([spectrum_n2,spectrum_n1,spectrum_b0],model_bb2_pl,model_bb2_pl_prior_list,parameters_bb2_pl,reference=False)
	mul_dir = savedir+'A_n'+str(i)+'_out/'
	if os.path.exists(mul_dir) ==False:
		os.makedirs(mul_dir)
	A = fit.run(outputfiles_basename=mul_dir+'A_n'+str(i)+'_',resume = False, verbose = True)
	
	equ_w = A.get_equal_weighted_posterior()
	c = A.get_stats()['modes'][0]
	'''
	mean = np.array(c['mean'])
	sigma = np.array(c['sigma'])
	best = np.array(c['maximum'])
	dx = mean - best
	sigmal = sigma - dx
	sigmah = sigma + dx
	kt1.append([best[1],sigmal[1],sigmah[1]])
	kt2.append([best[3],sigmal[3],sigmah[3]])
	a_list.append([best[5],sigmal[5],sigmah[5]])
	'''
	logkM,am,bm,logEcM = np.array(c['maximum'])
	mean = np.array(c['mean'])[3]
	sigma = np.array(c['sigma'])[3]
	Epm = (2+am)/(am-bm)*10**logEcM
	dx = mean - logEcM
	sigmal = sigma - dx
	sigmah = sigma + dx
	Ec = 10**logEcM
	E0_list.append([Ec,Ec-10**(logEcM-sigmal),10**(logEcM+sigmah)-Ec])
	
	loglike = equ_w[:,-1]
	
	a= equ_w[:,1]
	b = equ_w[:,2]
	log_Ec = equ_w[:,3]
	loglikemin = loglike.min()
	loglikemax = loglike.max()
	d_like = loglikemax-loglikemin
	vv = loglikemax-0.305*d_like
	good_index = np.where(loglike>vv)[0]
	a= a[good_index]
	b = b[good_index]
	log_Ec =log_Ec[good_index]
	Ep = (2+a)/(a-b)*10**log_Ec
	Ep_max = Ep.max()
	Ep_min = Ep.min()
	Ep_mean = 0.5*(Ep_max+Ep_min)
	Ep_sigma = 0.5*(Ep_max-Ep_min)
	dx = Ep_mean - Epm
	sigmal = Ep_sigma - dx
	sigmah = Ep_sigma + dx
	
	Ep_list.append([Epm,sigmal,sigmah])
	
	fig,ax = plt.subplots()
	fit.plot_model(A,n=2,ax = ax,reference = True)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim(10**-2,)
	ax.legend()
	fig.savefig(savedir + 'B_plot_model_'+str(i)+'.png')
	plt.close(fig)
	
	#fig,ax = plt.subplots()
	fig = plt.figure(constrained_layout=True,figsize=(8,11))
	fig.suptitle('Multinest fit of band model')
	gs = GridSpec(4, 1, figure=fig)
	ax1 = fig.add_subplot(gs[0:2])
	fit.plot_data(A,n=0,ax = ax1,reference = True)
	ax1.set_xscale('log')
	ax1.set_xlim(8,30000)
	ax1.set_yscale('log')
	#ax1.set_ylim(10**-2,)
	ax1.legend()
	ax1.tick_params(labelbottom=False)
	ax1.set_ylabel('spectrum (rate/kev)')
	
	ax2 = fig.add_subplot(gs[2])
	fit.plot_data_in_model(A,ax=ax2)
	ax2.set_xscale('log')
	ax2.set_xlim(8,30000)
	ax2.legend()
	ax2.tick_params(labelbottom=False)
	
	ax3 = fig.add_subplot(gs[3])
	fit.plot_data_in_reference(ax=ax3)
	ax3.set_xscale('log')
	ax3.set_xlim(8,30000)
	ax3.legend()
	ax3.set_xlabel('spectrum energy (kev)')
	fig.savefig(savedir + 'A_plot_data_'+str(i)+'.png')
	plt.close(fig)
	
	fit.plot_corner(A)
	plt.savefig(savedir + 'C_plot_corner_'+str(i)+'.png')
	plt.close()


myfile.printdatatofile(savedir+'D_Ep.txt',data = np.array(Ep_list).T,format = ['.6f']*3)
myfile.printdatatofile(savedir+'D_E0.txt',data = np.array(E0_list).T,format = ['.6f']*3)
'''
myfile.printdatatofile(savedir+'D_kt1.txt',data = np.array(kt1).T,format = ['.6f']*3)
myfile.printdatatofile(savedir+'D_kt2.txt',data = np.array(kt2).T,format = ['.6f']*3)
myfile.printdatatofile(savedir+'D_a.txt',data = np.array(a_list).T,format = ['.6f']*3)
'''
