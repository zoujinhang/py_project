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

savetop = '/home/laojin/my_work/lag/result3/'
#savetop = '/home/laojin/my_work/lag/bn150118409_1/'
data_top = '/media/laojin/Elements/trigdata/'
#sample_link = '/home/laojin/result/lag_samples2.txt'
sample_link = '/home/laojin/result/lag_samples6.txt'
#sample_link = '/home/laojin/result/idea.txt'

binsize = 0.2
e_band = get_band([10,1000],10,ovelap=0.5)
bins = np.arange(-10,300,binsize)

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
name,t_start,t_stop,year,ni,w_start,w_stop = myfile.readcol(sample_link)

NaI = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
BGO = ['b0','b1']

logt = []
logt_erl = []
logt_erh = []
B = []
B_erl = []
B_erh = []

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
	hl = files[ni[i]]
	
	trigtime = hl[0].header['TRIGTIME']
	time = hl[2].data.field(0)
	ch = hl[2].data.field(1)
	ch_n = hl[1].data.field(0)
	e1 = hl[1].data.field(1)
	e2 = hl[1].data.field(2)
	t = time - trigtime
	t,energy = ch_to_energy(t,ch,ch_n,e1,e2)
	results = get_lag([t,energy],e_band,bins,wind=[t_start[i],t_stop[i]],sigma=4,plot_savedir = savedir+'A_check/')
	fit = Lag_fit(model,model_bb2_pl_prior_list,parameters,result = results)
	#fit = Lag_fit(model2,model_lv_pl_prior_list,parameters_lv,result = results)
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
	ax1.set_xlim(10,1000)
	ax1.tick_params(labelbottom = False)
	ax1.set_ylabel('lag (s)')
	ax1.legend()
	
	ax2 = fig.add_subplot(gs[1])
	analy.plot_different(ax = ax2)
	ax2.set_xlim(10,1000)
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
     'B':np.array(B),
     'B_erl':np.array(B_erl),
     'B_erh':np.array(B_erh)}

pdata = pd.DataFrame(C)
pdata.to_csv(savetop+'B_lag_para.csv',index=False)
