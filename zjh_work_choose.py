

from Data_analysis import Clock
import Data_analysis.file as myfile
import numpy as np
import pandas as pd
import os
from Fermi_tool.daily import Database
import matplotlib
matplotlib.use('Agg')
import os
from Data_analysis.geometry import Geometry

topdir = '/media/laojin/Elements/daily/'
timestart = '2013-04-27T00:00:00'
timestop = '2013-04-27T00:59:00'


detector_list = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']

geometry = Geometry()
fermi = Database(topdir)
data = fermi.get_detector_data(timestart,timestop)

t_all = np.array([])

for ni_index in detector_list:

	ni = data[ni_index]
	event = ni['']



'''
topdir = '/home/laojin/my_work/SGRJ1935_search/SER/results/'
time_file_linke = '/home/laojin/source_search/SGRJ1935toa_time_range30.txt'

savedir = '/home/laojin/my_work/SGRJ1935_search/chooce2/'

if os.path.exists(savedir) == False:
	os.makedirs(savedir)

t_start_utc,t_stop_utc = myfile.readcol(time_file_linke)
name_list = os.listdir(topdir)
print(name_list)


Fermi_clock = Clock()
t_start_met = Fermi_clock.utc_to_met(t_start_utc)
t_stop_met = Fermi_clock.utc_to_met(t_stop_utc)

fig_links = []

bayesian_link = topdir + name_list[0] + '/SGRJ1935/A_save_bayesian.csv'
f_d = pd.read_csv(bayesian_link)

for i in f_d.index:
	fig_links.append(topdir + name_list[0] + '/SGRJ1935/good/A_bayes_'+str(i)+'.png')

for dirname in name_list[1:]:
	bayesian_link = topdir + dirname + '/SGRJ1935/A_save_bayesian.csv'
	f_di = pd.read_csv(bayesian_link)
	for i in f_di.index:
		fig_links.append(topdir + dirname + '/SGRJ1935/good/A_bayes_'+str(i)+'.png')
	f_d = f_d.append(pd.read_csv(bayesian_link))
fig_links = np.array(fig_links)
print(f_d)
print(f_d.index)
print('---------------------------------------')
index_i = (f_d['start_met']>=t_start_met[0])&(f_d['start_met']<=t_stop_met[0])
chooce_fig_links = list(fig_links[index_i])
print('index_i',index_i)
chooce_f_d = f_d[index_i]
for i in range(1,len(t_start_met)):
	index_i = (f_d['start_met']>=t_start_met[i])&(f_d['start_met']<=t_stop_met[i])
	chooce_fig_links = chooce_fig_links + list(fig_links[index_i])
	chooce_f_d = chooce_f_d.append(f_d[index_i])
print(chooce_f_d)
print(chooce_fig_links)
print('22',len(fig_links),len(chooce_fig_links))
for i,figlink in enumerate(chooce_fig_links):
	cmd = 'cp '+figlink + ' '+savedir + 'B_bayes_'+str(i)+'.png'
	os.system(cmd)

chooce_f_d.to_csv(savedir + 'A_time_list.csv',index=False)
'''












