
from Fermi_tool import *
import os
import numpy as np
import matplotlib.pyplot as plt
import Data_analysis.file as myfile
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage

savedir = '/home/laojin/my_lat/spectrum/bn150330828/'
data_top = '/media/laojin/Elements/trigdata/2015/bn150330828/'
sample_name = 'bn150330828'
if os.path.exists(savedir)  ==False:
	os.makedirs(savedir)

dt = 0.064
NaI = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
BGO = ['b0','b1']
files = get_file([data_top,sample_name],NaI,BGO)

hl = files['n2']
trigtime = hl[0].header['TRIGTIME']
time = hl[2].data.field(0)
ch = hl[2].data.field(1)
ch_n = hl[1].data.field(0)
e1 = hl[1].data.field(1)
e2 = hl[1].data.field(2)
t = time - trigtime

hl1 = files['b0']
trigtime1 = hl1[0].header['TRIGTIME']
time1 = hl1[2].data.field(0)
ch1 = hl1[2].data.field(1)
ch_n1 = hl1[1].data.field(0)
e11 = hl1[1].data.field(1)
e21 = hl1[1].data.field(2)
t1 = time1 - trigtime1


hl2 = files['n1']
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
lim_edges = edges[edges_index]

re_rate = result['re_hist'][0]
re_rate = np.concatenate((re_rate[:1],re_rate))
t_c,rate_c = result['lc']



spc,spc_err = get_spectrum(t,ch,ch_n,lim_edges,1)
spc1,spc_err1 = get_spectrum(t1,ch1,ch_n1,lim_edges,1)
spc2,spc_err2 = get_spectrum(t2,ch2,ch_n2,lim_edges,1)
for i in range(len(spc)):
	myfile.printdatatofile(savedir+'ZZ_spectrum_n2_'+str(i)+'.txt',data = [spc[i],spc_err[i]],format = ['.6f','.6f'])
	myfile.printdatatofile(savedir+'ZZ_spectrum_b0_'+str(i)+'.txt',data = [spc1[i],spc_err1[i]],format = ['.6f','.6f'])
	myfile.printdatatofile(savedir+'ZZ_spectrum_n1_'+str(i)+'.txt',data = [spc2[i],spc_err2[i]],format = ['.6f','.6f'])
lim_tc = 0.5*(lim_edges[1:]+lim_edges[:-1])
Ep,errl,errh = myfile.readcol('/home/laojin/my_lat/spectrum/A_spectrum/D_Ep.txt')
E0,errl0,errh0 = myfile.readcol('/home/laojin/my_lat/spectrum/A_spectrum/D_E0.txt')

fig, ax1 = plt.subplots()
ax1.plot(t_c,rate_c,label = 'light curve')
#for i in lim_edges:
#	ax1.axvline(x=i,color = 'r')
ax1.step(edges,re_rate,color = 'k',label = 'block')
ax1.set_ylabel('rate')
ax1.set_xlim(90,170)
ax2 = ax1.twinx()
ax2.errorbar(lim_tc,Ep,yerr = [errl,errh],color = 'g',ecolor = 'g',fmt = 'o')
ax2.errorbar(lim_tc,E0,yerr = [errl0,errh0],color = 'r',ecolor = 'r',fmt = 'o')
ax2.set_xlim(90,170)
ax2.set_ylim(30,10**3)
ax2.set_yscale('log')
ax2.set_ylabel('Ep Kev')
fig.savefig(savedir +'ZZ_lock_lc.png')
plt.close(fig)


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





