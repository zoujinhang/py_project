
from style_Shao_spectral_lag import readcol,transpose,GRB
import os
from multiprocessing import Pool


sample_link = '/home/laojin/result/lag_samples6.txt'

information = readcol(sample_link)
information = transpose(information)
#name,t_start,t_stop,year,ni,w_start,w_stop = readcol(sample_link)


def work(infor):
	
	name,t_start,t_stop,year,ni,w_start,w_stop = infor
	savetop = '/home/laojin/my_work/lag/result5/'
	data_top = '/media/laojin/Elements/trigdata/'
	savedir  = savetop + str(year) + '/' + name + '/'
	if os.path.exists(savedir) == False:
		os.makedirs(savedir)
	sample = GRB(name,data_top,savedir)
	sample.get_lag(ni,[w_start,w_stop],[t_start,t_stop],binsize = 0.2,sigma = 4,Positive_lag_positive_value = True)
	
p = Pool(2)
p.map(work,information)
p.close()
p.join()








