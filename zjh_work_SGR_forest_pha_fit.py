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
band_prior_list = [
	Prior([-5,4]),
	Prior([-3,8]),
	Prior([-5,2]),
	Prior([-1.3,10])
]
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

model_bb2_pl_prior_list = [
	Prior([-2,10]),
	Prior([1,4]),
	Prior([-2,10]),
	Prior([-2,3]),
	Prior([-5,5]),
	Prior([-3,5])
]

def model_bb_pl(e,cube):
	K1 = cube[0]
	kt1 = cube[1]
	#K2 = cube[2]
	#kt2 = cube[3]
	K3 = cube[2]
	a = cube[3]
	#print('kt1:',10**kt1,'kt2:',10**kt2,'K1',10**K1,'K2',10**K2)
	#if kt1<kt2:
	#	return [np.nan]
	return bb(e,10**K1,10**kt1)+ 10**K3*e**(-a)  #+bb(e,10**K2,10**kt2)

model_bb_pl_prior_list = [
	Prior([-2,10]),
	Prior([-1,4]),
	#Prior([-2,10]),
	#Prior([-2,3]),
	Prior([-5,5]),
	Prior([-3,5])
]


marker = 'tg_20200427_18330585'

datatop = '/home/laojin/my_work/' + marker + '/'
savedir = '/home/laojin/my_work/' + marker + '/'

if os.path.exists(savedir)  ==False:
	os.makedirs(savedir)


detectors = ['b1','n9','n6']
#detectors = ['b1','n9']
effective_band = [[600,30000],[8,800],[8,800],[8,800]]
time = []
SNR = []
spc = []
spc_bs = []

kt = []
a_ = []

kt1 = []
kt2 = []
a_list = []
for dete in detectors:

	spectrum_read = savedir + 'Z_pha_'+ dete +'.hdf5'
	detector = h5py.File(spectrum_read,'r')
	time.append(detector['time'][()])
	SNR.append(detector['SNR'][()])
	spc.append(detector['spc'][()])
	spc_bs.append(detector['spc_bs'][()])

time = time[0]
SNR = np.vstack(SNR).T

#dt = time[1]-time[0]
dt = 0.5

t_index = []
for index,snri in enumerate(SNR):
	snr_index = np.where(snri>=5)[0]
	if snr_index.size >= 2:
		t_index.append(index)
Ep_list = []
E0_list = []
new_time = []


spectrum_save = savedir + 'spectral_fit_BB_pl3/'
if os.path.exists(spectrum_save) ==False:
	os.makedirs(spectrum_save)
for i in t_index:
	print('time:',i)
	spc_list = []
	for index_,dete in enumerate(detectors):
		spci = spc[index_][i]
		bs_i = spc_bs[index_][i]
		bs_i[bs_i<0] = 0
		spci_err = np.sqrt(spci/dt)
		spci_err[spci_err<=0] =1
		rsp_link = datatop + 'B__'+dete + '.rsp'
		spectrum = Spectrum(spci,rsp_link,bs_spectrum=bs_i,spectrum_err=spci_err,time = dt,effective_band = effective_band[index_],spectrum_name = dete)
		spc_list.append(spectrum)

	#parameters = ['log K','a','b','log Ec']  #band
	#parameters_bb2_pl = ['log K1','log kt1','log K2','log kt2','log K3','a']#BB2_pl
	parameters_bb_pl = ['log K','log kt','log K2','a']#BB2_pl
	#fit = Fit(spc_list,band,band_prior_list,parameters,stats='pstat',reference=False)
	fit = Fit(spc_list,model_bb_pl,model_bb_pl_prior_list,parameters_bb_pl,stats='pstat',reference=False)
	mul_dir = spectrum_save+'A_n'+str(i)+'_out/'
	if os.path.exists(mul_dir) ==False:
		os.makedirs(mul_dir)
	A = fit.run(outputfiles_basename=mul_dir+'A_n'+str(i)+'_',resume = False, verbose = True)

	equ_w = A.get_equal_weighted_posterior()
	c = A.get_stats()['modes'][0]


	mean = np.array(c['mean'])
	sigma = np.array(c['sigma'])
	best = np.array(c['maximum'])
	dx = mean - best
	sigmal = sigma - dx
	sigmah = sigma + dx

	print('aaaaaaaaaaa',[time[i],best[1],sigmal[1],sigmah[1]])
	kt.append([time[i],best[1],sigmal[1],sigmah[1]])
	#kt2.append([time[i],best[3],sigmal[3],sigmah[3]])
	a_.append([time[i],best[3],sigmal[3],sigmah[3]])


	'''
	kt1.append([time[i],best[1],sigmal[1],sigmah[1]])
	kt2.append([time[i],best[3],sigmal[3],sigmah[3]])
	a_list.append([time[i],best[5],sigmal[5],sigmah[5]])
	
	
	
	logkM,am,bm,logEcM = np.array(c['maximum'])
	mean = np.array(c['mean'])[3]
	sigma = np.array(c['sigma'])[3]
	Epm = (2+am)/(am-bm)*10**logEcM
	dx = mean - logEcM
	sigmal = sigma - dx
	sigmah = sigma + dx
	Ec = 10**logEcM
	E0_list.append([time[i],Ec,Ec-10**(logEcM-sigmal),10**(logEcM+sigmah)-Ec])

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

	Ep_list.append([time[i],Epm,sigmal,sigmah])
	'''


	fig,ax = plt.subplots()
	fit.plot_model(A,n=2,ax = ax,reference = True)
	ax.set_ylabel(r'$KeV^2(Photons/cm^2/s/KeV)$')
	ax.set_xlabel('Energy KeV')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim(10**-2,)
	ax.legend()
	fig.savefig(spectrum_save + 'B_plot_model_'+str(i)+'.png')
	plt.close(fig)

	#fig,ax = plt.subplots()
	fig = plt.figure(constrained_layout=True,figsize=(8,11))
	#fig.suptitle('Multinest fit of band model. chi2 p = '+str(p))
	gs = GridSpec(4, 1, figure=fig)
	ax1 = fig.add_subplot(gs[0:2])
	fit.plot_data(A,n=0,ax = ax1,reference = True)
	ax1.set_xscale('log')
	ax1.set_xlim(8,30000)
	ax1.set_yscale('log')
	#ax1.set_ylim(10**-2,)
	ax1.legend()
	ax1.tick_params(labelbottom=False)
	ax1.set_ylabel('spectrum (rate/KeV)')

	ax2 = fig.add_subplot(gs[2])
	fit.plot_data_in_model(A,ax=ax2)
	ax2.set_xscale('log')
	ax2.set_ylabel('Normalized rate/KeV')
	ax2.set_xlim(8,30000)
	ax2.legend()
	ax2.tick_params(labelbottom=False)

	ax3 = fig.add_subplot(gs[3])
	fit.plot_data_in_reference(ax=ax3)
	ax3.set_xscale('log')
	ax3.set_xlim(8,30000)
	ax3.legend()
	ax3.set_ylabel('Normalized rate/KeV')
	ax3.set_xlabel('Energy (kev)')
	fig.savefig(spectrum_save + 'A_plot_data_'+str(i)+'.png')
	plt.close(fig)

	fit.plot_corner(A)
	plt.savefig(spectrum_save + 'C_plot_corner_'+str(i)+'.png')
	plt.close()


'''
myfile.printdatatofile(savedir+'D_Ep.txt',data = np.array(Ep_list).T,format = ['.6f']*4)
myfile.printdatatofile(savedir+'D_E0.txt',data = np.array(E0_list).T,format = ['.6f']*4)


myfile.printdatatofile(savedir+'D_kt1.txt',data = np.array(kt1).T,format = ['.6f']*4)
myfile.printdatatofile(savedir+'D_kt2.txt',data = np.array(kt2).T,format = ['.6f']*4)
myfile.printdatatofile(savedir+'D_a.txt',data = np.array(a_list).T,format = ['.6f']*4)
'''
myfile.printdatatofile(savedir+'D_kt.txt',data = np.array(kt).T,format = ['.6f']*4)
myfile.printdatatofile(savedir+'D_a.txt',data = np.array(a_).T,format = ['.6f']*4)
