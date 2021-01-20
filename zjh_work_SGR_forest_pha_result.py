
import Data_analysis.file as myfile
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes
import os
from daily_search.tool import overlap_bins
from daily_search.perception import Event
from daily_search.background import get_background_f,TD_baseline
from astropy.io import fits
import os
import pandas as pd
marker = 'tg_20200427_18330585'

datatop = '/home/laojin/my_work/' + marker + '/'
savedir = '/home/laojin/my_work/' + marker + '/'

if os.path.exists(savedir)  ==False:
	os.makedirs(savedir)

detectors = ['b1','n9','n6','n0']
binsize = 0.5
stepsize = 0.1
bins = overlap_bins([-15,30],binsize,stepsize)
lc_binsize = 0.064
bins_lc = np.arange(-17,31,lc_binsize)

time1,a_,a_l,a_h = myfile.readcol(datatop + 'D_a.txt')
time1,logkt,logkt_l,logkt_h = myfile.readcol(datatop + 'D_kt.txt')


plt.errorbar(time1,a_,yerr = [a_l,a_h],elinewidth=1,capsize=2,alpha=0.5,fmt = '.',label="alpha")
plt.savefig(savedir + 'Z_bb_pl_alpha.png')
plt.close()


plt.errorbar(time1,logkt,yerr = [logkt_l,logkt_h],elinewidth=1,capsize=2,alpha=0.5,fmt = '.',label="log kt")
plt.savefig(savedir + 'Z_bb_pl_logkt.png')
plt.close()


dete = 'n9'

link = datatop + marker +'_' + dete + '.fits'
hl = fits.open(link)
trigtime = hl[0].header['TRIGTIME']

time = hl[2].data.field(0)
ch = hl[2].data.field(1)
ch_n = hl[1].data.field(0)
e_ch = np.vstack([ch_n,ch_n]).T
e1 = hl[1].data.field(1)
e2 = hl[1].data.field(2)
t = time - trigtime
ch_E = pd.DataFrame({'CHANNEL':ch_n,'E_MIN':e1,'E_MAX':e2})
lightcurve = Event(ch_E,t,ch)

overlap_lc_t,overlap_lc_rates = lightcurve.ovelap_lightcurve(bins)
lc_t,lc_rates = lightcurve(bins = bins_lc)

bs_f = get_background_f(lc_t,lc_rates)
bs = bs_f(overlap_lc_t)

'''
fig = plt.figure()

host = HostAxes(fig, [0.15, 0.1, 0.65, 0.8])
par1 = ParasiteAxes(host, sharex=host)
par2 = ParasiteAxes(host, sharex=host)
host.parasites.append(par1)
host.parasites.append(par2)

host.axis["right"].set_visible(False)

par1.axis["right"].set_visible(True)
par1.axis["right"].major_ticklabels.set_visible(True)
par1.axis["right"].label.set_visible(True)

par2.axis["right2"] = par2.new_fixed_axis(loc="right", offset=(60, 0))

fig.add_axes(host)

p1, = host.plot(overlap_lc_t,overlap_lc_rates,color = 'k',label = 'overlap lc')
p2, = host.plot(overlap_lc_t,bs,color = 'orange',label = 'overlap bs')
p3, = par1.errorbar(time1,a_,yerr = [a_l,a_h],elinewidth=1,capsize=2,alpha=0.5,fmt = '.',label="alpha")
p4, = par2.errorbar(time1,logkt,yerr = [logkt_l,logkt_h],elinewidth=1,capsize=2,alpha=0.5,fmt = '.',label="log kt")

#host.set_xlim(0, 2)
#host.set_ylim(0, 2)
#par1.set_ylim(0, 4)
#par2.set_ylim(1, 65)

host.set_xlabel("Time (s)")
host.set_ylabel("count rate n/s")
par1.set_ylabel("alpha")
par2.set_ylabel("log kt (KeV)")

host.legend()

host.axis["left"].label.set_color(p1.get_color())
#host.axis["left"].label.set_color(p2.get_color())
par1.axis["right"].label.set_color(p3.get_color())
par2.axis["right2"].label.set_color(p4.get_color())

fig.savefig(savedir + 'Z_LC.png')
plt.close()
'''


fig, ax_f = plt.subplots()
ax_c = ax_f.twinx()
# automatically update ylim of ax2 when ylim of ax1 changes.
#ax_f.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
ax_f.plot(overlap_lc_t,overlap_lc_rates,color = 'k',label = 'overlap lc')
ax_f.plot(overlap_lc_t,bs,color = 'orange',label = 'overlap bs')
ax_c.errorbar(time1,a_,yerr = [a_l,a_h],elinewidth=1,capsize=2,alpha=0.5,fmt = '.',label="alpha (PL)")
#ax_f.set_xlim(0, 100)
#ax_f.set_title('Two scales: Fahrenheit and Celsius')
ax_f.set_ylabel('time (s)')
#ax_f.legend()
ax_c.set_ylabel('alpha')
ax_c.legend()
fig.savefig(savedir + 'Z_lc_alpha.png')
plt.close()

fig, ax_f = plt.subplots()
ax_c = ax_f.twinx()
# automatically update ylim of ax2 when ylim of ax1 changes.
#ax_f.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
ax_f.plot(overlap_lc_t,overlap_lc_rates,color = 'k',label = 'overlap lc')
ax_f.plot(overlap_lc_t,bs,color = 'orange',label = 'overlap bs')
ax_c.errorbar(time1,logkt,yerr = [logkt_l,logkt_h],elinewidth=1,capsize=2,alpha=0.5,fmt = '.',label="log kt (BB)")
#ax_f.set_xlim(0, 100)
#ax_f.set_title('Two scales: Fahrenheit and Celsius')
ax_f.set_xlabel('time (s)')
#ax_f.legend()
ax_c.set_ylabel('log kt')
ax_c.set_ylim(0.5,1.2)
ax_c.legend()
fig.savefig(savedir + 'Z_lc_logkt.png')
plt.close()
