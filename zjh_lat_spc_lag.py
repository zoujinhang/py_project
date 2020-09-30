
import  Data_analysis.file as myfile
from Data_analysis.Baseline import TD_baseline,TD_bs,WhittakerSmooth
from astropy.io import fits
import os
import numpy as np
import matplotlib.pyplot as plt
from Data_analysis import ch_to_energy
from Fermi_tool.lag import get_band,get_event_with_band,Lag_plot,get_lag
from scipy.interpolate import interp1d
from Fermi_tool.lag import Prior,Lag_fit



datalink = '/home/laojin/trigdata/bn190530430/glg_tte_n1_bn190530430_v00.fit'
savedir = '/home/laojin/my_lat/spc_lag/'

binsize = 0.2
e_band = get_band([10,1000],12,ovelap=0.1)

def model(e,para):
	
	dt = -10**para[0]*(e**(-para[1]) - e[0]**(-para[1]))
	
	return dt

model_bb2_pl_prior_list = [
	
	Prior([-1.5,1.5]),
	Prior([-0.2,1.5])
]
parameters = ['logt','B']

if os.path.exists(savedir) == False:
	os.makedirs(savedir)

hl = fits.open(datalink)
trigtime = hl[0].header['TRIGTIME']
ch_n = hl[1].data['CHANNEL']
e_min = hl[1].data['E_MIN']
e_max = hl[1].data['E_MAX']

time = hl[2].data['TIME']
ch = hl[2].data['PHA']
t = time - trigtime


print(e_band)
t,energy = ch_to_energy(t,ch,ch_n,e_min,e_max)
bins = np.arange(-50,100,binsize)


results = get_lag([t,energy],e_band,bins,sigma=4,plot_savedir = savedir+'A_check/')
fit = Lag_fit(model,model_bb2_pl_prior_list,parameters,result = results)
mul_dir = savedir+'A_n_out/'
if os.path.exists(mul_dir) ==False:
	os.makedirs(mul_dir)
analy = fit.run(outputfiles_basename=mul_dir+'A_n_',resume = False, verbose = True)
analy.plot_corner()
plt.savefig(savedir+'A_corner.png')
plt.close()


lgplt = Lag_plot(results)

lgplt.plot_lightcurve(sigma=4)
plt.savefig(savedir+'D_lightcurve.png')
plt.close()

ax = lgplt.plot_lag()
analy.plot_fit(ax = ax)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('energy Kev')
ax.set_ylabel('lag (s)')
ax.legend()
plt.savefig(savedir + 'E_lag.png')
plt.close()

'''
t_c = 0.5*(bins[1:]+bins[:-1])
event_band = get_event_with_band(t,energy,e_band)

t_i,e_i = event_band[0]
num_i = np.histogram(t_i,bins=bins)[0]
rate_i = num_i/binsize

cs,bs,sigma = TD_bs(t_i,rate_i,sigma=True,lambda_=200,it_ = 1)
n = len(t_c)
lag_t = np.arange(-n+1,n,1)*binsize
plt.plot(t_c,rate_i)
plt.plot(t_c,bs)
plt.plot(t_c,bs+sigma)
plt.savefig(savedir + 'B_lightcurve'+str(0)+'.png')
plt.close()


for index,(t_i,e_i) in enumerate(event_band[1:]):
	num_i = np.histogram(t_i,bins=bins)[0]
	rate_i = num_i/binsize
	cs0 = cs
	cs,bs,sigma = TD_bs(t_c,rate_i,sigma=True)
	
	
	nccf = np.correlate(cs0 / cs0.max(), cs / cs.max(), 'full')
	w = np.ones(len(nccf))
	smoo_nccf = WhittakerSmooth(nccf,w,0.15/binsize**1.5)
	intef = interp1d(lag_t,smoo_nccf,kind = 'quadratic')
	new_lat_t = np.arange(lag_t[0]+1,lag_t[-1]+0.001-1,0.001)
	#print(lag_t)
	#print(new_lat_t)
	plt.plot(lag_t,nccf)
	#plt.plot(lag_t,smoo_nccf)
	plt.plot(new_lat_t,intef(new_lat_t))
	plt.xlim(-10,10)
	plt.savefig(savedir + 'C_nccf'+str(index)+'.png')
	plt.close()
	
	plt.plot(t_c,rate_i)
	plt.plot(t_c,bs)
	plt.plot(t_c,bs+sigma)
	plt.savefig(savedir + 'B_lightcurve'+str(index+1)+'.png')
	plt.close()
num = np.histogram(t, bins=bins)[0]
rate = num/binsize

plt.plot(t_c,rate)
plt.savefig(savedir + 'A_lightcurve.png')
plt.close()


'''






