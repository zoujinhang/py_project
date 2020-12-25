
import numpy as np
from GECAM_geometry.satellite import sky_map

import Data_analysis.file as myfile
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import os
from scipy.optimize import curve_fit


def gausion(x,sigma,K):

	return K * np.exp(-0.5*(x/sigma)**2)




link1 = '/home/laojin/my_lat/skyp/tru_locat.txt'
link2 = '/home/laojin/my_lat/skyp/hitl_locat.txt'
savedir = '/home/laojin/my_lat/skyp/'

if os.path.exists(savedir) ==False:
	os.makedirs(savedir)

tru_ra,tru_dec = myfile.readcol(link1)

hitl_ra,hitl_dec = myfile.readcol(link2)

tru = SkyCoord(ra = tru_ra,dec = tru_dec,frame = 'icrs',unit = 'deg')
hitl = SkyCoord(ra = hitl_ra,dec = hitl_dec,frame = 'icrs',unit = 'deg')


sep = []

for i in range(len(tru_ra)):

	sep.append(tru[i].separation(hitl[i]).deg)

sep = np.array(sep)
#print(sep)

bins = np.linspace(0,30,29)
d_deg = bins[1:]-bins[:-1]

bin_n = np.histogram(sep,bins = bins)[0]
rate = bin_n/(d_deg*2*np.pi*bins[1:])

popt, pcov = curve_fit(gausion, bins[1:], rate)
print('sigma :',popt[0])


rate = np.concatenate((rate[:1],rate))

x = np.linspace(0,30,100)
plt.figure(figsize=(5,5))
plt.plot(bins,rate,label = 'data')
plt.plot(x,gausion(x,popt[0],popt[1]),label = 'gausion')
plt.xlabel(r'$\theta$ (degree)')
plt.ylabel(r'$number \cdot degree^{-2}$')
plt.savefig(savedir+'A_pdf.png')
plt.close()
