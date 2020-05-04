
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

edges_index = np.where((edges>=140)&(edges<=160))[0]
lim_edges = edges[edges_index]

re_rate = result['re_hist'][0]
re_rate = np.concatenate((re_rate[:1],re_rate))
t_c,rate_c = result['lc']

plt.plot(t_c,rate_c,label = 'light curve')
for i in lim_edges:
	plt.axvline(x=i,color = 'r')
plt.step(edges,re_rate,color = 'k',label = 'block')
plt.xlim(80,200)
plt.savefig(savedir +'ZZ_lock_lc.png')
plt.close()

spc,spc_err = get_spectrum(t,ch,ch_n,lim_edges,1)
myfile.printdatatofile(savedir+'ZZ_spectrum.txt',data = [spc[2],spc_err[2]],format = ['.6f','.6f'])


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





