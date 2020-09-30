from Fermi_tool.burst import *
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import Data_analysis.file as myfile
from Data_analysis import ch_to_energy
from Fermi_tool.lag import get_band,Lag_plot,get_lag
from Fermi_tool.lag import Prior,Lag_fit

savetop = '/home/laojin/my_work/lag/result6/'
#savetop = '/home/laojin/my_work/lag/bn150118409_1/'
data_top = '/media/laojin/Elements/trigdata/'
#sample_link = '/home/laojin/result/lag_samples2.txt'
sample_link = '/home/laojin/result/lag_samples7.txt'
#sample_link = '/home/laojin/result/idea.txt'

binsize = 0.2
e_band = get_band([10,30000],20,ovelap=0.5)
#bins = np.arange(-10,300,binsize)

cont_e = 200


def model(e,para):
	
	dt = -10**para[0]*(e**(-para[1]) - e[0]**(-para[1]))
	
	return dt


model_bb2_pl_prior_list = [
	
	Prior([-1.5, 1.5]),
	Prior([-0.2, 1.5])
]
parameters = ['logt', 'B']



def model2(e,para):
	d = 6
	dt = 10**para[0]*(e**para[1]-e[0]**para[1]) - (d-3)*(e**(d-4) - e[0]**(d-4))*para[2]
	return dt

model_lv_pl_prior_list = [
	Prior([-10, 2]),
	Prior([0, 1.1]),
	Prior([-1, 1])
]
parameters_lv = ['log t', 'a','K']



#name,t_start,t_stop,year,ni = myfile.readcol(sample_link)
name,t_start,t_stop,year,ni,bi,w_start,w_stop = myfile.readcol(sample_link)

NaI = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
BGO = ['b0','b1']

logt = []
logt_erl = []
logt_erh = []
B = []
B_erl = []
B_erh = []
K = []
K_erl = []
K_erh = []

for i in range(len(name)):
	
	bins = np.arange(w_start[i],w_stop[i],binsize)
	sample_link = data_top + str(year[i]) + '/' + name[i] + '/'
	savedir  = savetop + str(year[i]) + '/' + name[i] + '/'
	save_all = savetop + str(year[i]) + '/A_all/'
	if os.path.exists(save_all) == False:
		os.makedirs(save_all)
	if os.path.exists(savedir) == False:
		os.makedirs(savedir)
	files = get_file([sample_link,name[i]],NaI,BGO)
	
	
	hl_ni = files[ni[i]]
	trigtime_ni = hl_ni[0].header['TRIGTIME']
	time_ni = hl_ni[2].data.field(0)
	ch_ni = hl_ni[2].data.field(1)
	ch_n_ni = hl_ni[1].data.field(0)
	e1_ni = hl_ni[1].data.field(1)
	e2_ni = hl_ni[1].data.field(2)
	
	t_ni = time_ni - trigtime_ni
	t_ni,energy_ni = ch_to_energy(t_ni,ch_ni,ch_n_ni,e1_ni,e2_ni)
	
	ni_index = np.where(energy_ni<=cont_e)[0]
	t_ni = t_ni[ni_index]
	energy_ni = energy_ni[ni_index]
	
	hl_bi = files[bi[i]]
	trigtime_bi = hl_bi[0].header['TRIGTIME']
	time_bi = hl_bi[2].data.field(0)
	ch_bi = hl_bi[2].data.field(1)
	ch_n_bi = hl_bi[1].data.field(0)
	e1_bi = hl_bi[1].data.field(1)
	e2_bi = hl_bi[1].data.field(2)
	
	t_bi = time_bi - trigtime_bi
	t_bi,energy_bi = ch_to_energy(t_bi,ch_bi,ch_n_bi,e1_bi,e2_bi)
	bi_index = np.where(energy_bi>cont_e)[0]
	t_bi = t_bi[bi_index]
	energy_bi = energy_bi[bi_index]
	t = np.concatenate((t_ni,t_bi))
	energy = np.concatenate((energy_ni,energy_bi))
	sort_index = np.argsort(t)
	t = t[sort_index]
	energy = energy[sort_index]
	
	results = get_lag([t,energy],e_band,bins,wind=[t_start[i],t_stop[i]],sigma=4,plot_savedir = savedir+'A_check/')
	#fit = Lag_fit(model,model_bb2_pl_prior_list,parameters,result = results)
	fit = Lag_fit(model2,model_lv_pl_prior_list,parameters_lv,result = results)
	mul_dir = savedir+'A_n_out/'
	if os.path.exists(mul_dir) ==False:
		os.makedirs(mul_dir)
	analy = fit.run(outputfiles_basename=mul_dir+'A_n_',resume = False, verbose = True)
	
	paras,errorl,errorh = analy.get_best_fit()
	logt.append(paras[0])
	logt_erl.append(errorl[0])
	logt_erh.append(errorh[0])
	B.append(paras[1])
	B_erl.append(errorl[1])
	B_erh.append(errorh[1])
	K.append(paras[2])
	K_erl.append(errorl[2])
	K_erh.append(errorh[2])
	
	
	analy.plot_corner()
	plt.savefig(savedir+'A_corner.png')
	plt.close()
	
	lgplt = Lag_plot(results)
	
	lgplt.plot_lightcurve(sigma=4)
	plt.savefig(savedir+'B_lightcurve.png')
	plt.savefig(save_all+'B_'+name[i]+'_lightcurve.png')
	plt.close()
	
	fig = plt.figure(constrained_layout=True,figsize = (5,5))
	gs = GridSpec(2, 1, figure=fig)
	ax1 = fig.add_subplot(gs[0])
	lgplt.plot_lag(ax = ax1)
	analy.plot_fit(ax=ax1)
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax1.set_xlim(20,30000)
	ax1.tick_params(labelbottom = False)
	ax1.set_ylabel('lag (s)')
	ax1.legend()
	
	ax2 = fig.add_subplot(gs[1])
	analy.plot_different(ax = ax2)
	ax2.set_xlim(20,30000)
	ax2.set_xlabel('energy Kev')
	ax2.set_ylabel('residual')
	ax2.set_xscale('log')
	ax2.legend()
	plt.savefig(savedir + 'C_lag.png')
	plt.savefig(save_all+'C_'+name[i]+'_lag.png')
	plt.close()
	
C = {'name':name,
     'logt':np.array(logt),
     'logt_erl':np.array(logt_erl),
     'logt_erh':np.array(logt_erh),
     'a':np.array(B),
     'a_erl':np.array(B_erl),
     'a_erh':np.array(B_erh),
     'K':np.array(K),
     'K_erl':np.array(K_erl),
     'K_erh':np.array(K_erh)}

pdata = pd.DataFrame(C)
pdata.to_csv(savetop+'B_lag_para.csv',index=False)
