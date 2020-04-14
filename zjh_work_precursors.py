import numpy as np
import matplotlib.pyplot as plt
import os
from Calculate_duration import accumulate_counts,Plot
import Data_analysis.file as myfile
from Data_analysis import TD_baseline,ch_to_energy,background_correction,get_bayesian_duration,WhittakerSmooth
from astropy.io import fits
from astropy.stats import bayesian_blocks
import pandas as pd


#infor_link = '/home/laojin/my_work/zbb_Precursors/list.txt'
#infor_link = '/home/laojin/my_work/zbb_Precursors/list1.txt'
infor_link = '/home/laojin/my_work/zbb_Precursors/list3.txt'#bn150922234
#infor_link = '/home/laojin/my_work/zbb_Precursors/list2.txt'#bn191221802
savedir = '/home/laojin/my_work/zbb_Precursors/result/'

sample_list,dete_list,year_list,el,eu,tl,tu = myfile.readcol(infor_link)

data_top = '/media/laojin/Elements/trigdata/'
dt = 0.064
#start = -50#bn191221802
start = -20
stop = 10
for i in range(len(sample_list)):
	print(sample_list[i])
	datalink = data_top+str(year_list[i]) + '/'+sample_list[i] + '/'
	filename = myfile.findfile(datalink,'glg_tte_'+dete_list[i]+'_'+sample_list[i]+'_v*')[0]
	edges_file = 'GRB'+ sample_list[i][2:]+dete_list[i]+'_Block.txt'
	print(edges_file)
	hl = fits.open(datalink + filename)
	trigtime = hl[0].header['TRIGTIME']
	time = hl[2].data.field(0)
	ch = hl[2].data.field(1)
	ch_n = hl[1].data.field(0)
	e1 = hl[1].data.field(1)
	e2 = hl[1].data.field(2)
	t = time - trigtime
	t, e = ch_to_energy(t, ch, ch_n, e1, e2)
	e_index = np.where((e >= el[i]) & (e <= eu[i]))[0]
	t = t[e_index]
	edges_name = myfile.findfile('/home/laojin/my_work/zbb_Precursors/',edges_file)
	print(edges_name)
	if edges_name is not False and len(edges_name)>0:
		edges = myfile.readcol('/home/laojin/my_work/zbb_Precursors/'+edges_name[0])[0]
		edges = np.array(edges)
		dt = 0.01
		
		bins = np.arange(edges[0],edges[-1]+dt,dt)
		bin_n,bin_edges = np.histogram(t,bins = bins)
		rate = bin_n/dt
		t_c = (bin_edges[1:]+bin_edges[:-1])*0.5
		t_c,cs_rate,bs_rate = TD_baseline(t_c,rate)
		rate_sm = cs_rate+bs_rate.mean()
		bin_n_sm = np.round(rate_sm*dt)
	else:
		dt = 0.064
		bins = np.arange(start,stop+dt,dt)
		bin_n,bin_edges = np.histogram(t,bins = bins)
		rate = bin_n/dt
		t_c = (bin_edges[1:]+bin_edges[:-1])*0.5
		t_c,cs_rate,bs_rate = TD_baseline(t_c,rate)
		rate_sm = cs_rate+bs_rate.mean()
		bin_n_sm = np.round(rate_sm*dt)
		#edges = bayesian_blocks(t_c,bin_n_sm,fitness='events',gamma = np.exp(-3))#bn191221802
		edges = bayesian_blocks(t_c,bin_n_sm,fitness='events',gamma = np.exp(-5))
		dt = 0.01
		bins = np.arange(start,stop+dt,dt)
		bin_n,bin_edges = np.histogram(t,bins = bins)
		rate = bin_n/dt
		t_c = (bin_edges[1:]+bin_edges[:-1])*0.5
		t_c,cs_rate,bs_rate = TD_baseline(t_c,rate)
		rate_sm = cs_rate+bs_rate.mean()
		bin_n_sm = np.round(rate_sm*dt)
		
		
	print(edges)
	if os.path.exists(savedir) == False:
		os.makedirs(savedir)
	#result = background_correction(t_c,rate_sm,edges,degree = 7,plot_save=savedir + 'Z_'+sample_list[i]+'_check.png')#bn191221802
	#result = background_correction(t_c,rate_sm,edges,degree = 30,plot_save=savedir + 'Z_'+sample_list[i]+'_check.png')#bn101208498
	#result = background_correction(t_c,rate_sm,edges,degree = 7,plot_save=savedir + 'Z_'+sample_list[i]+'_check.png')
	result = background_correction(t_c,rate_sm,edges,degree = 40,plot_save=savedir + 'Z_'+sample_list[i]+'_check.png')#bn150922234
	c_rate = result['lc'][1]
	sigma = result['bkg'][2]
	re_rate = result['re_hist'][0]
	#startedges,stopedges = get_bayesian_duration(result,sigma = 2)#bn191221802
	#startedges,stopedges = get_bayesian_duration(result,sigma = 5)#bn101208498
	startedges,stopedges = get_bayesian_duration(result,sigma = 5)
	w = np.ones(len(t_c))
	for ij in range(len(startedges)):
		#index_w = np.where((t_c>=startedges[ij]-dt)&(t_c<=stopedges[ij]+dt))[0]#bn101208498
		index_w = np.where((t_c>=startedges[ij])&(t_c<=stopedges[ij]))[0]
		w[index_w] = 0
	
	#result1 = accumulate_counts(t_c,c_rate*dt,np.sqrt(bin_n),w,startedges,stopedges,txx = 0.9,it = 100,lamd = 40/dt)#bn191221802
	result1 = accumulate_counts(t_c,c_rate*dt,np.sqrt(bin_n),w,startedges,stopedges,txx = 0.9,it = 300,lamd = 10/dt)#bn150922234
	#result1 = accumulate_counts(t_c,c_rate*dt,np.sqrt(bin_n),w,startedges-dt,stopedges+dt,txx = 0.9,it = 300,lamd = 20/dt)#bn101208498
	#result1 = accumulate_counts(t_c,c_rate*dt,np.sqrt(bin_n),w,startedges,stopedges,txx = 0.9,it = 300,lamd = 70/dt)
	txx = result1['txx']
	txx_err1 ,txx_err2 = result1['txx_err']
	t1 = result1['t1']
	t1_err1,t1_err2 = result1['t1_err']
	t2 = result1['t2']
	t2_err1,t2_err2 = result1['t2_err']
	news = {'t90':txx,"t90_err-":txx_err1,"t90_err+":txx_err2,'t90_1':t1,"t90_1_err-":t1_err1,
	        "t90_1_err+":t1_err2,'t90_2':t2,"t90_2_err-":t2_err1,"t90_2_err+":t2_err2}
	df = pd.DataFrame(news,columns=['t90',"t90_err-","t90_err+",'t90_1',"t90_1_err-", "t90_1_err+",'t90_2',"t90_2_err-","t90_2_err+"])
	df.to_csv(savedir + 'C_'+sample_list[i]+'_T90.csv',index=False,float_format='%.2f')
	result1['time_edges'] = [startedges,stopedges]
	result1['t_c'] = t_c
	result1['rate'] = c_rate
	result1['bs'] = WhittakerSmooth(c_rate,w,lambda_=100/dt)
	result1['good'] = True
	result1['sigma'] = sigma
	result1['bayesian_edges'] = [edges]
	result1['bayesian_rate'] = [np.concatenate((re_rate[:1], re_rate))]
	if os.path.exists(savedir) == False:
		os.makedirs(savedir)
	myplt = Plot(result1)
	plt.title('GRB'+sample_list[i][2:])
	myplt.plot_light_curve(sigma=4.5)
	plt.xlim(tl[i],tu[i])
	plt.savefig(savedir + 'A_'+sample_list[i]+'_lightcurve.png')
	plt.savefig(savedir + 'E_'+sample_list[i]+'_lightcurve.eps')
	plt.close()
	
	for ij in range(len(result1['txx'])):
		plt.title(sample_list[i])
		myplt.plot_distribution('90',num = ij)
		plt.savefig(savedir + 'D_'+sample_list[i]+'_distribution_'+str(ij)+'.png')
		plt.close()
		
	plt.figure(figsize = (10,10))
	plt.subplot(2,1,1)
	plt.title('GRB'+sample_list[i][2:])
	myplt.plot_Txx1('90')
	plt.xlim(tl[i],tu[i])
	plt.subplot(2,1,2)
	myplt.plot_Txx2('90')
	plt.xlim(tl[i],tu[i])
	plt.savefig(savedir + 'B_'+sample_list[i]+'_txx.png')
	plt.savefig(savedir + 'F_'+sample_list[i]+'_txx.eps')
	plt.close()











