import re
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Data_analysis import Clock
from astropy.time import Time


datadir = '/home/laojin/results/source_search/new_version/SGRJ1935/'
savedir = '/home/laojin/my_lat/daily_search/SGRJ1935_statistics/'
if os.path.exists(savedir) ==False:
	os.makedirs(savedir)

porise_list = os.listdir(datadir)

print(porise_list)

fermiclo = Clock()

time_start = []
time_stop = []

for proise in porise_list:

	sampledir = datadir + proise + '/'

	samplelike = os.listdir(sampledir)
	sample = []

	for samplelikei in samplelike:

		if os.path.isdir(sampledir + samplelikei):
			name_add =re.split('[_]',samplelikei)[-1]
			sampledatadir = sampledir + samplelikei + '/Z_candidate_'+name_add+'.csv'
			data_i = pd.read_csv(sampledatadir)
			dete_name = ['n0_SNR','n1_SNR','n2_SNR','n3_SNR','n4_SNR','n5_SNR','n6_SNR','n7_SNR','n8_SNR','n9_SNR','na_SNR','nb_SNR']
			time = data_i['met'].values
			data = data_i[dete_name].values

			find_start = True
			t_start = 0
			for index,d_i in enumerate(data):
				indexi = np.where(d_i>=5)[0]
				if find_start and indexi.size>=2:
					t_start = time[index]
					find_start = False
				elif indexi.size < 2 and find_start == False:
					time_start.append(t_start)
					time_stop.append(time[index])
wt = 0.05
time_start = np.array(time_start)
time_stop = np.array(time_stop)

frb200428 = '2020-04-28T14:34:33.04672'
mdj_frb200428 = Time(frb200428).mjd

higt_de = '2020-04-27T18:26:20.00'
mjd_hight_de = Time(higt_de).mjd


higt_de1 = '2020-04-27T18:33:04'
mjd_hight_de11 = Time(higt_de1).mjd


print(time_stop)
time_start = np.floor(time_start / wt) * wt
time_stop = np.ceil(time_stop / wt) * wt

time_start,un_index = np.unique(time_start,return_index=True)
time_stop = time_stop[un_index]

sort_index = np.argsort(time_start)
time_start = time_start[sort_index]
time_stop = time_stop[sort_index]
duration = time_stop - time_start


mdjstart = fermiclo.met_to_utc(time_start).mjd
print('the number of candidates:',len(mdjstart))
t_range = np.arange(mdjstart.min(),mdjstart.max()+1,1)

n = np.histogram(mdjstart,bins = t_range)[0]
n = np.concatenate((n[:1],n))
plt.step(t_range,n,color = 'k',label = 'Counts of SGRJ1958 bursts')
plt.axvline(x = mdj_frb200428,color = 'r',label = 'FRB200428 burst time')
#plt.axvline(x = mjd_hight_de,color = 'orange',label = 'a very high burst density on 2020 April 27')
plt.xlabel('Time (mjd)')
plt.ylabel('Counts')
plt.legend()
plt.savefig(savedir + 'A_counts.png')
plt.close()


plt.step(t_range,n,color = 'k',label = 'Counts of SGRJ1958 bursts')
plt.axvline(x = mdj_frb200428,color = 'r',label = 'FRB200428 burst time')
plt.axvline(x = mjd_hight_de,color = 'orange',label = 'A very high burst density on 2020 April 27')
plt.axvline(x = mjd_hight_de11,color = 'g',label = 'The long-GRB-like candidate')
plt.xlabel('Time (mjd)')
plt.ylabel('Counts')
plt.legend()
plt.xlim(58950,59000)
plt.savefig(savedir + 'A_counts2.png')
plt.close()



wit_t = time_start[1:]-time_start[:-1]
print(len(wit_t))
logbin = np.logspace(np.log10(wit_t.min()),np.log10(wit_t.max()),20)
print(logbin)
n = np.histogram(wit_t,bins = logbin)[0]

n = np.concatenate((n[:1],n))

plt.step(logbin,n)
plt.xscale('log')
plt.ylabel('Counts')
plt.xlabel('Waiting time (s)')
plt.savefig(savedir + 'A_writting.png')
plt.close()

plt.plot(mdjstart[1:],wit_t,'.--',color = 'k',label='Waiting time evolution')
plt.axvline(x = mdj_frb200428,color = 'r',label = 'FRB200428 burst time')
plt.xlabel('Time (mjd)')
plt.xlim(58966,58970)
plt.ylabel('Log waiting time (s)')
plt.yscale('log')
plt.legend()
plt.savefig(savedir + 'A_writting_evolution.png')
plt.close()

plt.plot(mdjstart,duration,'.',color = 'k',label='duration evolution')
plt.axvline(x = mdj_frb200428,color = 'r',label = 'FRB200428 burst time')
plt.xlabel('Time (mjd)')
plt.ylabel('Log duration (s)')
plt.yscale('log')
plt.legend()
plt.savefig(savedir + 'A_duration_evolution.png')
plt.close()

