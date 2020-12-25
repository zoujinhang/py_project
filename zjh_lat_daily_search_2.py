
from daily_search.satellite import Detectors,Geometry,Locate,Sky_map
from daily_search.data_base import Database
from daily_search.perception import Event
from daily_search.background import get_background_f
#from Fermi_tool.daily import Database
from Data_analysis import ch_to_energy,TD_bs,TD_baseline
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt

import os

Fermi_dete = Detectors()

print(Fermi_dete('n0'))
topdir = '/media/laojin/TOSHIBA_EXT/daily/'
trigtime = 608633290
source = SkyCoord(ra=6.100 ,dec= -32.000,frame ='icrs',unit='deg')
name = 'GRB200415367'

#source = SkyCoord(ra=293.729,dec= 21.3864,frame ='icrs',unit='deg')      #you can change these things
#name = 'SGRJ1935'


#name = 'bn180720598'
#source = SkyCoord(ra=0.53,dec= -2.93,frame ='icrs',unit='deg')
#trigtime = 553789304.654332

#source = SkyCoord(ra=54.5068,dec= -26.9467,frame ='icrs',unit='deg')      #you can change these things
#name = 'GRB190114C'
#trigtime = 569192227

#name = 'bn190110726'
#source = SkyCoord(ra=276.939,dec= -53.643,frame ='icrs',unit='deg')
#trigtime = 568833894

savedir = '/home/laojin/my_lat/location/'+name+'/'
if os.path.exists(savedir)==False:
	os.makedirs(savedir)


#topdir = '/home/laojin/daily_data/'
#topdir = '/media/laojin/TOSHIBA_EXT/daily/'
#topdir = '/media/laojin/My_Passport/Fermi_GBM_Dailiy_Data/'
#timestart = '2020-04-28T00:19:40.00'
#timestop = '2020-04-28T00:19:47.00'

fermi = Database(topdir)
clock = fermi.clock

#timestart = '2018-07-20T14:19:23'
#timestop = '2018-07-20T14:29:37'
#met_start = clock.utc_to_met(timestart)+1
#met_stop = clock.utc_to_met(timestop)-1



met_start = trigtime - 100
met_stop = trigtime + 151

timestart = clock.met_to_utc(met_start)
timestop = clock.met_to_utc(met_stop)
num = -1
#met_start=609725986.924219+1
#met_stop=609725991.924219-1

binszie = 0.05
binsize_baseline = 1
binsize_lightcurve = 0.05

bins_baseline = np.arange(met_start,met_stop,binsize_baseline)
bins_lightcurve = np.arange(met_start,met_stop,binsize_lightcurve)
bins_lightcurve_c = 0.5*(bins_lightcurve[1:]+bins_lightcurve[:-1])


bins = np.arange(met_start,met_stop,binszie)
bin_c = 0.5*(bins[1:]+bins[:-1])

indexi = np.where((bin_c>=trigtime-1)&(bin_c<=trigtime+6))[0]
t_v = bin_c[indexi]

#during = t_v.max()-t_v.min()


detector_list = fermi.detector
pos = fermi.get_poshist_data(timestart,timestop)
data = fermi.get_detector_data(timestart,timestop)

dete_rate = []
dete_bs = []
t_m = t_v.mean()


for ni in fermi.detector:

	ni_ch_n = data[ni]['ch_E']
	ni_event = data[ni]['events']
	time = ni_event['TIME'].values
	ch = ni_event['PHA'].values
	lightcurve = Event(ni_ch_n,time,ch)
	t_cc,rate_b = lightcurve(bins = bins_baseline,energy_band=[33,84],channel = True)
	rate_b = np.concatenate((rate_b[:1],rate_b,rate_b[-1:]))
	t_cc = np.concatenate(([met_start],t_cc,[met_stop]))
	bs_f = get_background_f(t_cc,rate_b)
	t_cc,rate_c = lightcurve(bins = bins_lightcurve,energy_band=[33,84],channel= True)
	bs = bs_f(t_cc)



	#cs,bs = TD_bs(bin_c,rate)
	dete_rate.append(rate_c)
	dete_bs.append(bs)
	plt.plot(bin_c,rate_c,label='lightcurve')
	plt.plot(bin_c,bs,label = 'bs')
	plt.axvline(x=t_v[0],color = 'g')
	plt.axvline(x=t_v[num],color = 'r')
	plt.legend()
	plt.savefig(savedir + 'Z_'+ni+'_lightcurve.png')
	plt.close()

print('rate\n',dete_rate)
print('bs \n',dete_bs)

dete_rate = (np.vstack(dete_rate).T)[indexi]
dete_bs = (np.vstack(dete_bs).T)[indexi]

sigma = np.sqrt(dete_bs/binszie)

SNR = (dete_rate - dete_bs)/sigma

index_list = []

for i in range(t_v.size):
	index = np.where(SNR[i]>=5)[0]
	if len(index) >= 3:
		index_list.append(i)

#during = binsize_lightcurve*len(index_list)
during = 1
dete_rate_good = dete_rate[index_list]

dete_bs_good = dete_bs[index_list]

m_rate = dete_rate_good.mean(axis = 0)
m_bs = dete_bs_good.mean(axis = 0)

detector_list = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
#detector_list = ['n7','n8']
Fermi = Geometry(pos)
print('time_band',Fermi.met_time_band)
print('get_pos',Fermi.get_pos(t_m))
print('get_earth_position',Fermi.get_earth_point(t_m))
print('qsj',Fermi.get_qsj(t_m))
print('detect',Fermi(t_m))
locate = Locate(Fermi)

print('m_rate',np.round(m_rate*during))
print('m_bs',np.round(m_bs*during))

#ra,dec,err_r,ra_rcat,dec_rcat,chi2 = locate.condidate(t_m,np.round(m_rate*during),np.round(m_bs*during),'case3',[44,300],detector_list = detector_list)



ra,dec,err_r,ra_rcat,dec_rcat,chi2,P,ra_rcat2,dec_rcat2 = locate.condidate(t_m,m_rate,m_bs,during,'case4',[44,300],detector_list = detector_list)

print('chi2 max',chi2.max())
print('chi2 min',chi2.min())
print('location',ra,dec,err_r)
print('during',during)
print('P max',P.max())

P2 = P/P.sum()
P_sort = np.sort(-P)
c = (chi2-chi2.min())/(chi2.max()-chi2.min())
index_x = np.where(c<1)[0]

sort_p = np.sort(P2)[::-1]

p_sum = 0
h_1sigma = None
h_2sigma = None
h_3sigma = None
for i in -P_sort:
	p_sum = i+p_sum
	if p_sum >= 0.6826 and h_1sigma is None:
		h_1sigma = i
	if p_sum >= 0.9545 and h_2sigma is None:
		h_2sigma = i
	if p_sum >= 0.9974 and h_3sigma is None:
		h_3sigma = i
print('sigma list',[h_3sigma,h_2sigma,h_1sigma])

line = plt.tricontourf(ra_rcat[index_x],dec_rcat[index_x],chi2[index_x]-chi2.min(),[0,2.3,4.61,9.21])
print(line)
plt.savefig(savedir+'A_contour.png')
plt.close()

smp = Sky_map(figsize=(10,5))
smp.add_subplot(1,1,1)
smp.add_source(source,name)
smp.scatter(ra_rcat[index_x],dec_rcat[index_x],marker = ',',c = c[index_x],s = 5)
smp.tricontour(ra_rcat[index_x],dec_rcat[index_x],chi2[index_x]-chi2.min(),[0,2.3,4.61,9.21])
#smp.scatter(ra_rcat[index_x],dec_rcat[index_x],marker = ',',c = P2[index_x],s = 5)
#smp.tricontour(ra_rcat[index_x],dec_rcat[index_x],P2[index_x],[h_3sigma,h_2sigma,h_2sigma])
smp.plot_earth(t_v[num],Fermi)
smp.plot_detector(t_v[num],Fermi)
smp.plot(ra,dec,'*',markersize = 10,color = 'orange')
smp.plot(54.5068,-26.9467,',',markersize = 2,color = 'r')
smp.savefig(savedir+'A_sky_map_n_101.png')
smp.close()


smp = Sky_map(figsize=(10,5))
smp.add_subplot(1,1,1)
smp.add_source(source,name)
#smp.scatter(ra_rcat[index_x],dec_rcat[index_x],marker = ',',c = c[index_x],s = 5)
#smp.tricontour(ra_rcat[index_x],dec_rcat[index_x],chi2[index_x]-chi2.min(),[0,2.3,4.61,9.21])
#smp.scatter(ra_rcat2,dec_rcat2,marker = ',',c = P,s = 5)
smp.tricontour(ra_rcat2,dec_rcat2,P,[h_3sigma,h_2sigma,h_1sigma])
smp.plot_earth(t_v[num],Fermi)
smp.plot_detector(t_v[num],Fermi)
smp.plot(ra,dec,'*',markersize = 10,color = 'orange')
smp.plot(54.5068,-26.9467,',',markersize = 2,color = 'r')
smp.savefig(savedir+'A_sky_map_P2_101.png')
smp.close()