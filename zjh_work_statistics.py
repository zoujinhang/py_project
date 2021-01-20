import re
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Data_analysis import Clock
from astropy.time import Time


datadir = '/home/laojin/results/source_search/new_version/SGRJ1935/'
savedir = '/home/laojin/my_lat/daily_search/SGRJ1935_statistics2/'
if os.path.exists(savedir) ==False:
	os.makedirs(savedir)

porise_list = os.listdir(datadir)

fermiclo = Clock()

time_start = []
time_stop = []

print(porise_list)
candidate_time = []

for i in porise_list:
	if os.path.isdir(datadir + i):
		year = i[3:7]
		month = i[7:9]
		day = i[9:11]
		hour = i[12:14]
		minute = i[14:16]
		second = i[16:18]
		dot = i[18:]
		utc = year + '-' + month + '-' + day + 'T' + hour + ':' + minute + ':' + second + '.' + dot
		candidate_time.append(utc)

candidate_time = Time(candidate_time)

wt = 0.05
time_start = np.array(time_start)
#time_stop = np.array(time_stop)

frb200428 = '2020-04-28T14:34:33.04672'
mdj_frb200428 = Time(frb200428).mjd

higt_de = '2020-04-27T18:26:20.00'
mjd_hight_de = Time(higt_de).mjd

dust_scattering = '2020-04-27T19:41:00'
mjd_dust_scattering = Time(dust_scattering).mjd

higt_de1 = '2020-04-27T18:33:04'
mjd_hight_de11 = Time(higt_de1).mjd

SGR_frb_list1 = [59148.495776886891690,
		 59148.497543767036404,
		 59148.497656593630381,
		 59148.498672654941888,
		 59148.499949716475385,
		 59148.500626730659860,
		 59148.500664259125188,
		 59148.506190907726705,
		 59148.512356245613773,
		 59148.513484796014382,
		 59148.513709780192585,
		 59148.514010301434610,
		 59148.514537224531523,
		 59148.515100902397535,
		 59148.515890604372544,
		 59148.527694602802512,
		 59148.529199002165115]
SGR_frb_list1 = Time(SGR_frb_list1,format='mjd',scale = 'utc').mjd

SGR_frb_list2 = [59151.443853909790050,
		 59151.445696001595934,
		 59151.446222920138098,
		 59151.446748591675714,
		 59151.446936563937925,
		 59151.448027064776397,
		 59151.448139534106303,
		 59151.448177019337891,
		 59151.449004031317600,
		 59151.450018983705377,
		 59151.450583180405374,
		 59151.451861956149514,
		 59151.452125283438363,
		 59151.452538169345644,
		 59151.452613169378310,
		 59151.453026703813521,
		 59151.453402887280390,
		 59151.453440349752782,
		 59151.455170049594017,
		 59151.455319690125179,
		 59151.455395680030051,
		 59151.455582587332174,
		 59151.455846210446907,
		 59151.456673950604454,
		 59151.457650275442575,
		 59151.458251991498400,
		 59151.458289754336874,
		 59151.458477722058888,
		 59151.458928630214359,
		 59151.459718657592020,
		 59151.460132105559751,
		 59151.460169247176964,
		 59151.460620683261368,
		 59151.461410610514577,
		 59151.462124336241686,
		 59151.462462472416519,
		 59151.463027265308483,
		 59151.464154088564101,
		 59151.464342934639717,
		 59151.464642896098667,
		 59151.464680877397768,
		 59151.465395765932044,
		 59151.465508041845169,
		 59151.465658330911538,
		 59151.467161851905985,
		 59151.467612534783257,
		 59151.468177204791573,
		 59151.468816239954322,
		 59151.469192025200755,
		 59151.472538534573687,
		 59151.472876017665840,
		 59151.472952141826681,
		 59151.473140196008899,
		 59151.474266787154193,
		 59151.475244852794276,
		 59151.475620763194456,
		 59151.478026728560508,
		 59151.478515442802745,
		 59151.479530472563056,
		 59151.480282277436245,
		 59151.481522225381923,
		 59151.482988468225813,
		 59151.483514940759051,
		 59151.483665420972102,
		 59151.484078802946897,
		 ]
SGR_frb_list2 = Time(SGR_frb_list2,format='mjd',scale = 'utc').mjd

SGR_frb_list3 = [59152.462500209643622,
		59152.464042387771769,
		59152.464530030229071,
		59152.466222547489451,
		59152.467049372877227,
		59152.467274386639474,
		59152.468967538778088,
		59152.471561651633237,
		59152.472913798119407,
		59152.473477828702016,
		59152.473628368075879,
		59152.474191672750749,
		59152.474418823258020,
		59152.475356996132177,
		59152.475695444059966,
		59152.478102415225294,
		59152.478441011066025,
		59152.478853805943800,
		59152.479305788168858,
		59152.480432800293784,
		59152.481335422155098,
		59152.481447743579338,
		59152.481636125448858,
		59152.484116772269772,
		59152.484154419063998,
		59152.484418262902182,
		59152.488440323279065,
		59152.490695749271254,
		59152.192951564281656]


SGR_frb_list3 = Time(SGR_frb_list3,format='mjd',scale = 'utc').mjd


#print(time_stop)
#time_start = np.floor(time_start / wt) * wt
#time_stop = np.ceil(time_stop / wt) * wt

#time_start,un_index = np.unique(time_start,return_index=True)
#time_stop = time_stop[un_index]

#sort_index = np.argsort(time_start)
#time_start = time_start[sort_index]
#time_stop = time_stop[sort_index]
#duration = time_stop - time_start


#mdjstart = fermiclo.met_to_utc(time_start).mjd
mdjstart = candidate_time.mjd
mdjstart = np.sort(mdjstart)
print('the number of candidates:',len(mdjstart))


t_range = np.arange(int(mdjstart.min())-2,int(mdjstart.max())+5,1)
n = np.histogram(mdjstart,bins = t_range)[0]
n = np.concatenate((n[:1],n))

t_range_frb_list1 = np.arange(int(SGR_frb_list1.min())-1,int(SGR_frb_list1.max())+1,0.05)
n_frb_list1 = np.histogram(SGR_frb_list1,bins = t_range_frb_list1)[0]
n_frb_list1 = np.concatenate((n_frb_list1[:1],n_frb_list1))

t_range_frb_list2 = np.arange(int(SGR_frb_list2.min())-1,int(SGR_frb_list2.max())+1,0.05)
n_frb_list2 = np.histogram(SGR_frb_list2,bins = t_range_frb_list2)[0]
n_frb_list2 = np.concatenate((n_frb_list2[:1],n_frb_list2))

t_range_frb_list3 = np.arange(int(SGR_frb_list3.min())-1,int(SGR_frb_list3.max())+1,0.05)
n_frb_list3 = np.histogram(SGR_frb_list3,bins = t_range_frb_list3)[0]
n_frb_list3 = np.concatenate((n_frb_list3[:1],n_frb_list3))

plt.step(t_range,n/1,color = 'k',label = 'Outbreak rate of SGRJ1958 bursts')
plt.step(t_range_frb_list1,n_frb_list1/0.05,label = 'Outbreak rate of FAST detection 1')
plt.step(t_range_frb_list2,n_frb_list1/0.05,label = 'Outbreak rate of FAST detection 2')
plt.step(t_range_frb_list3,n_frb_list1/0.05,label = 'Outbreak rate of FAST detection 3')
plt.axvline(x = mdj_frb200428,color = 'r',label = 'SGR-FRB200428 burst time')
#plt.axvline(x = mjd_hight_de,color = 'orange',label = 'a very high burst density on 2020 April 27')
plt.xlabel('Time (mjd)')
plt.ylabel('Outbreak rate   number/day')
plt.legend()
plt.savefig(savedir + 'A_counts.png')
plt.close()

t_range = np.arange(int(mdjstart.min())-2,int(mdjstart.max())+5,0.05)
n = np.histogram(mdjstart,bins = t_range)[0]
n = np.concatenate((n[:1],n))
plt.step(t_range,n/0.05,color = 'k',label = 'Outbreak rate of SGRJ1958 bursts')



plt.axvline(x = mdj_frb200428,color = 'r',label = 'SGR-FRB200428 burst time')
#plt.axvline(x = mjd_hight_de,color = 'orange',label = 'A very high burst density on 2020 April 27')
plt.axvline(x = mjd_hight_de11,color = 'g',label = 'The long-GRB-like candidate')
plt.axvline(x = mjd_dust_scattering,color = 'orange',label = 'dust scattering halo')
plt.xlabel('Time (mjd)')
plt.ylabel('Outbreak rate   number/day')
plt.legend()
plt.xlim(58966,58968)
plt.savefig(savedir + 'A_counts2.png')
plt.close()



#wit_t = time_start[1:]-time_start[:-1]
#print(len(wit_t))
#logbin = np.logspace(np.log10(wit_t.min()),np.log10(wit_t.max()),20)
#print(logbin)
#n = np.histogram(wit_t,bins = logbin)[0]

#n = np.concatenate((n[:1],n))

#plt.step(logbin,n)
#plt.xscale('log')
#plt.ylabel('Counts')
#plt.xlabel('Waiting time (s)')
#plt.savefig(savedir + 'A_writting.png')
#plt.close()

#plt.plot(mdjstart[1:],wit_t,'.--',color = 'k',label='Waiting time evolution')
#plt.axvline(x = mdj_frb200428,color = 'r',label = 'FRB200428 burst time')
#plt.xlabel('Time (mjd)')
#plt.xlim(58966,58970)
#plt.ylabel('Log waiting time (s)')
#plt.yscale('log')
#plt.legend()
#plt.savefig(savedir + 'A_writting_evolution.png')
#plt.close()

#plt.plot(mdjstart,duration,'.',color = 'k',label='duration evolution')
#plt.axvline(x = mdj_frb200428,color = 'r',label = 'FRB200428 burst time')
#plt.xlabel('Time (mjd)')
#plt.ylabel('Log duration (s)')
#plt.yscale('log')
#plt.legend()
#plt.savefig(savedir + 'A_duration_evolution.png')
#plt.close()

