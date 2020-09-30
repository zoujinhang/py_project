
from Fermi_tool.burst import *
from Fermi_tool.spectrum import get_spectrum
import os
import numpy as np
import matplotlib.pyplot as plt
import Data_analysis.file as myfile
from Data_analysis import overlap_bins
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import h5py


savedir = '/home/laojin/my_lat/spectrum/bn150330828/'
data_top = '/media/laojin/Elements/trigdata/2015/bn150330828/'
sample_name = 'bn150330828'

#savedir = '/home/laojin/my_lat/spectrum/bn090809978/'
#data_top = '/media/laojin/Elements/trigdata/2009/bn090809978/'
#sample_name = 'bn090809978'

if os.path.exists(savedir)  ==False:
	os.makedirs(savedir)

dt = 0.064
NaI = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
BGO = ['b0','b1']
files = get_file([data_top,sample_name],NaI,BGO)

hl = files['n2']
#hl = files['n4']

trigtime = hl[0].header['TRIGTIME']
time = hl[2].data.field(0)
ch = hl[2].data.field(1)
ch_n = hl[1].data.field(0)
e1 = hl[1].data.field(1)
e2 = hl[1].data.field(2)
t = time - trigtime

hl1 = files['b0']
#hl1 = files['b0']

trigtime1 = hl1[0].header['TRIGTIME']
time1 = hl1[2].data.field(0)
ch1 = hl1[2].data.field(1)
ch_n1 = hl1[1].data.field(0)
e11 = hl1[1].data.field(1)
e21 = hl1[1].data.field(2)
t1 = time1 - trigtime1


hl2 = files['n1']
#hl2 = files['n5']
trigtime2 = hl2[0].header['TRIGTIME']
time2 = hl2[2].data.field(0)
ch2 = hl2[2].data.field(1)
ch_n2 = hl2[1].data.field(0)
e12 = hl2[1].data.field(1)
e22 = hl2[1].data.field(2)
t2 = time2 - trigtime2

#ch_index = np.where((ch>=3)&(ch<123))[0]
#ch_n1 = np.arange(3,123,1,dtype = int)
#t = t[ch_index]
#ch= ch[ch_index]
bins = np.arange(t[0],t[-1],dt)
bin_n,bin_edges = np.histogram(t,bins = bins)
t_c = (bin_edges[1:] + bin_edges[:-1])*0.5
rate = bin_n/dt
t_c,cs_rate,bs_rate = TD_baseline(t_c,rate)
rate_sm = cs_rate+bs_rate.mean()
bin_n_sm = np.round(rate_sm*dt)
edges = bayesian_blocks(t_c,bin_n_sm,fitness='events',p0 = 0.05)
result = background_correction(t_c,rate_sm,edges,degree = 7)

edges_index = np.where((edges>=90)&(edges<=170))[0]
#edges_index = np.where((edges>=-10)&(edges<=30))[0]
lim_edges = edges[edges_index]

re_rate = result['re_hist'][0]
re_rate = np.concatenate((re_rate[:1],re_rate))
t_c,rate_c = result['lc']

t_c_index = np.where((t_c>=90)&(t_c<=170))[0]


index_t_c = np.where(rate_c[t_c_index]>=2500)[0]
t_c_band = (t_c[t_c_index])[index_t_c]
print(t_c_band.min(),t_c_band.max())

#bins_arr = overlap_bins((-10,30),binsize = 1.,stepsize = 0.2)
bins_arr = overlap_bins((90,170),binsize = 1.,stepsize = 0.1)
#spc,spc_err = get_spectrum(t,ch,ch_n,lim_edges,1)
#spc1,spc_err1 = get_spectrum(t1,ch1,ch_n1,lim_edges,1)
#spc2,spc_err2 = get_spectrum(t2,ch2,ch_n2,lim_edges,1)

spc,spc_err = get_spectrum(t,ch,ch_n,bin_arr = bins_arr,bg_dt = 1.0)

all_rate_n4 = spc.sum(axis = 1)
all_rate_n4_err = np.sqrt((spc_err**2).sum(axis = 1))

f1 = h5py.File(savedir+'ZZ_spectrum_n2.hdf5','w')
#f1 = h5py.File(savedir+'ZZ_spectrum_n4.hdf5','w')
spc_ = f1.create_dataset('spc',spc.shape,dtype = np.float)
spc_[...] = spc
spc_arr_ = f1.create_dataset('spc_err',spc_err.shape,dtype = np.float)
spc_arr_[...] = spc_err
f1.close()

spc1,spc_err1 = get_spectrum(t1,ch1,ch_n1,bin_arr = bins_arr,bg_dt = 1.0)

all_rate_b0 = spc1.sum(axis = 1)
all_rate_b0_err = np.sqrt((spc_err1**2).sum(axis = 1))
f2 = h5py.File(savedir+'ZZ_spectrum_b0.hdf5','w')
#f2 = h5py.File(savedir+'ZZ_spectrum_b0.hdf5','w')
spc_ = f2.create_dataset('spc',spc1.shape,dtype = np.float)
spc_[...] = spc1
spc_arr_ = f2.create_dataset('spc_err',spc_err1.shape,dtype = np.float)
spc_arr_[...] = spc_err1
f2.close()

spc2,spc_err2 = get_spectrum(t2,ch2,ch_n2,bin_arr = bins_arr,bg_dt = 1.0)

all_rate_n5 = spc2.sum(axis = 1)
all_rate_n5_err = np.sqrt((spc_err2**2).sum(axis = 1))
f3 = h5py.File(savedir+'ZZ_spectrum_n1.hdf5','w')
#f3 = h5py.File(savedir+'ZZ_spectrum_n5.hdf5','w')
spc_ = f3.create_dataset('spc',spc2.shape,dtype = np.float)
spc_[...] = spc2
spc_arr_ = f3.create_dataset('spc_err',spc_err2.shape,dtype = np.float)
spc_arr_[...] = spc_err2
f3.close()

'''
for i in range(len(spc)):
	myfile.printdatatofile(savedir+'ZZ_spectrum_n2_'+str(i)+'.txt',data = [spc[i],spc_err[i]],format = ['.6f','.6f'])
	myfile.printdatatofile(savedir+'ZZ_spectrum_b0_'+str(i)+'.txt',data = [spc1[i],spc_err1[i]],format = ['.6f','.6f'])
	myfile.printdatatofile(savedir+'ZZ_spectrum_n1_'+str(i)+'.txt',data = [spc2[i],spc_err2[i]],format = ['.6f','.6f'])
'''

#lim_tc = 0.5*(lim_edges[1:]+lim_edges[:-1])
lim_tc = 0.5*(bins_arr[:,0]+bins_arr[:,-1])
Ep,errl,errh = myfile.readcol('/home/laojin/my_lat/spectrum/A_spectrum1/D_Ep.txt')
E0,errl0,errh0 = myfile.readcol('/home/laojin/my_lat/spectrum/A_spectrum1/D_E0.txt')

#Ep,errl,errh = myfile.readcol('/home/laojin/my_lat/spectrum/A_spectrum_bn090809978/D_Ep.txt')
#E0,errl0,errh0 = myfile.readcol('/home/laojin/my_lat/spectrum/A_spectrum_bn090809978/D_E0.txt')

#kt1,kt1errl,kt1errh = myfile.readcol('/home/laojin/my_lat/spectrum/A_spectrum_bb/D_kt1.txt')
#kt2,kt2errl,kt2errh = myfile.readcol('/home/laojin/my_lat/spectrum/A_spectrum_bb/D_kt2.txt')

lim_t_index = np.where((lim_tc>=t_c_band.min())&(lim_tc<=t_c_band.max()))[0]

print(lim_tc[lim_t_index])


Ep = np.array(Ep)
errl = np.array(errl)
errh = np.array(errh)
E0 = np.array(E0)
errl0 = np.array(errl0)
errh0 = np.array(errh0)




fig, ax1 = plt.subplots(figsize = (20,10),constrained_layout=True)
ax1.plot(t_c,rate_c,label = 'light curve')
#for i in lim_edges:
#	ax1.axvline(x=i,color = 'r')

ax1.step(edges,re_rate,color = 'k',label = 'block')
ax1.set_ylabel('rate')
ax1.set_xlabel('time (s)')
#ax1.set_xlim(-10,30)
#ax1.set_ylim(0,5000)
ax1.set_xlim(90,170)
ax2 = ax1.twinx()


print(len(lim_tc),len(Ep))
#index_tc = np.where((lim_tc>=124)&(lim_tc<=163))[0]
#print(index_tc)
#ax2.errorbar(lim_tc[index_tc],Ep[index_tc],yerr = [errl[index_tc],errh[index_tc]],color = 'g',ecolor = 'g',fmt = 'o',alpha = 0.3,label = 'band Ep')
#ax2.errorbar(lim_tc[index_tc],E0[index_tc],yerr = [errl0[index_tc],errh0[index_tc]],color = 'r',ecolor = 'r',fmt = 'o',alpha = 0.3,label = 'band E0')
ax2.errorbar(lim_tc[lim_t_index],Ep[lim_t_index],yerr = [errl[lim_t_index],errh[lim_t_index]],color = 'g',ecolor = 'g',fmt = 'o',alpha = 0.3,label = 'band Ep')
ax2.errorbar(lim_tc[lim_t_index],E0[lim_t_index],yerr = [errl0[lim_t_index],errh0[lim_t_index]],color = 'r',ecolor = 'r',fmt = 'o',alpha = 0.3,label = 'band E0')
#ax2.errorbar(lim_tc[:len(kt1)],10**np.array(kt1),yerr = [np.log(10)*10**np.array(kt1errl),np.log(10)*10**np.array(kt1errh)],color = '#f58220',ecolor = '#f58220',fmt = 'o',alpha = 0.3,label = 'bb+pl kt1')
#ax2.errorbar(lim_tc[:len(kt2)],10**np.array(kt2),yerr = [np.log(10)*10**np.array(kt2errl),np.log(10)*10**np.array(kt2errh)],color = '#ea66a6',ecolor = '#ea66a6',fmt = 'o',alpha = 0.3,label = 'bb+pl kt2')

#ax2.set_xlim(90,170)
#ax2.set_ylim(1,10**4)
ax2.set_yscale('log')
ax2.set_ylabel('Energy Kev')
ax2.legend()
fig.savefig(savedir +'ZZ_lock_lc.png')
plt.close(fig)

plt.errorbar(all_rate_n4[lim_t_index],Ep[lim_t_index],xerr = all_rate_n4_err[lim_t_index],yerr = [errl[lim_t_index],errh[lim_t_index]],fmt = '-o',alpha = 0.5,label = 'Flux Ep n2')
plt.errorbar(all_rate_n5[lim_t_index],Ep[lim_t_index],xerr = all_rate_n5_err[lim_t_index],yerr = [errl[lim_t_index],errh[lim_t_index]],fmt = '-o',alpha = 0.5,label = 'Flux Ep n1')
plt.legend()
plt.xlabel('flux counts/s')
plt.ylabel('Ep')
plt.yscale('log')
plt.xscale('log')
#plt.ylim(0,150)
#plt.plot(all_rate_n4[lim_t_index],Ep[lim_t_index])
plt.savefig(savedir + 'A_correlation.png')
plt.close()


'''
spc_bins = np.arange(t[0],t[-1],0.5)
spc_bins_c = (spc_bins[1:]+spc_bins[:-1])*0.5
spc,spc_err = get_spectrum(t,ch,ch_n1,spc_bins,1)

e_c = np.sqrt(e1*e2)
e_c_lim = e_c[ch_n1]
plt.title('0')
plt.errorbar(e_c_lim,spc[0],yerr = spc_err[0])
plt.xlabel('energy kev')
plt.ylabel('rate')
plt.xscale('log')
plt.savefig(savedir +'ZZ_spc_one.png')
plt.close()

nn = len(spc)
print('number of spc',nn)
duration = nn/24
fig,ax = plt.subplots()
n_ = 0
def make_frame(t):
	n = int(t*24)
	ax.clear()
	ax.set_title('time %.1f'% spc_bins_c[n])
	ax.errorbar(e_c_lim,spc[n],yerr = spc_err[n])
	ax.set_xlabel('energy kev')
	ax.set_ylabel('rate')
	ax.set_xscale('log')
	ax.set_ylim(-10,160)
	return mplfig_to_npimage(fig)

animation = VideoClip(make_frame, duration=duration)
animation.write_videofile(savedir+"ZZ_spc_evolution.mp4", fps=24,codec = 'h264')
animation.close()
'''





