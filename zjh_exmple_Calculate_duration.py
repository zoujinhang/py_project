
import numpy as np
import matplotlib.pyplot as plt
import os
import Data_analysis.file as myfile
from Data_analysis import ch_to_energy
from Calculate_duration import Plot,get_txx,save_result
from astropy.io import fits


infor_link = '/home/laojin/my_work/zbb_Precursors/list.txt'
savedir = '/home/laojin/my_work/text_for_tool/result/'
data_top = '/media/laojin/Elements/trigdata/'


sample_list,dete_list,year_list,el,eu,tl,tu = myfile.readcol(infor_link)
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
	if os.path.exists(savedir) == False:
		os.makedirs(savedir)
	#计算
	result = get_txx(t,binsize = 0.01,time_edges=[start,stop],background_degree=7,sigma = 5,txx = 0.9,it = 300,prior = 5,plot_check=savedir + 'Z_'+sample_list[i]+'_check.png',hardnss=100)
	#----------------------------------------------------------------------------
	#保存结果
	save_result(result,savedir + 'C_'+sample_list[i]+'_T90.csv',float_format='%.3f')
	#----------------------------------------------------------------------------
	#绘图
	
	myplt = Plot(result)
	plt.title('GRB'+sample_list[i][2:])
	myplt.plot_light_curve(sigma=4.5)#光变
	plt.xlim(tl[i],tu[i])
	plt.savefig(savedir + 'A_'+sample_list[i]+'_lightcurve.png')
	plt.close()

	for ij in range(len(result['txx'])):
		plt.title(sample_list[i])
		myplt.plot_distribution('90',num = ij)#画出mcmc采样分布
		plt.savefig(savedir + 'D_'+sample_list[i]+'_distribution_'+str(ij)+'.png')
		plt.close()
	#--------------------------------------------------------------------------
	#画出txx
	plt.figure(figsize = (10,10))
	plt.subplot(2,1,1)
	plt.title('GRB'+sample_list[i][2:])
	myplt.plot_Txx1('90')#
	plt.xlim(tl[i],tu[i])
	plt.subplot(2,1,2)
	myplt.plot_Txx2('90')
	plt.xlim(tl[i],tu[i])
	plt.savefig(savedir + 'B_'+sample_list[i]+'_txx.png')
	plt.close()
