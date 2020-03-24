import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits

savedir = '/home/laojin/my_work/heb171223818/result/'
hebdatalink = '/home/laojin/my_work/heb171223818/data/HEB171223818_HE-Evt_detec09.fits'


binsize = 0.04


if os.path.exists(savedir) == False:
	os.makedirs(savedir)

hl =  fits.open(hebdatalink)
heb_time = hl[0].data.T[0]


heb_time  = np.sort(heb_time)
heb_time = heb_time-188681898.00
print(heb_time)
print(heb_time[-1]-heb_time[0])

heb_edges = np.arange(-10,10+binsize,binsize)
heb_n,heb_edges = np.histogram(heb_time,bins = heb_edges)
heb_c = (heb_edges[1:]+heb_edges[:-1])*0.5
heb_rate = heb_n/binsize


plt.plot(heb_c,heb_rate)
plt.xlim(-10,10)
plt.savefig(savedir+ 'A_heb_lightcurve.png')
plt.close()








