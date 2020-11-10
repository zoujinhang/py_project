import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import Data_analysis.file as myfile
import cartopy.crs as ccrs


filelink = '/home/laojin/my_lat/location/locrates_1deg_50_300_hard_n.dat'
alocdat_comp_link = '/home/laojin/my_lat/location/alocdat_comp.dat'
savedir = '/home/laojin/my_lat/location/'
if os.path.exists(savedir) ==False:
	os.makedirs(savedir)


response_res = 20


grid_spacing=190./(180/np.pi)/(response_res-1)
for i in range(20):
	print('grid_point',grid_spacing*i/np.pi*180)


scatterdata = myfile.readcol(alocdat_comp_link)
scatterdata = np.vstack(scatterdata)
re_scatterdata = scatterdata.reshape((16,236,19))
#print('scatterdata:\n',scatterdata)
print('scatterdata shape',scatterdata.shape)
#print('reshape scatterdata\n',re_scatterdata)
print('reshape scatterdata shape',re_scatterdata.shape)
print('re_scatterdata[0]',re_scatterdata[0])
print('re_scatterdata[0] shape',re_scatterdata[0].shape)


dat_ = myfile.readcol(filelink)
print('data\n',dat_[0])
dtorad = 180.0/np.arccos(-1.0)
print(dtorad)

print('sin',np.sin(90/dtorad))

ra = np.array(dat_[0])/60#/dtorad
dec = 90-np.array(dat_[1])/60#/dtorad

npoints = 0
for i in range(1,180):
	dist_at_zen = np.round(np.sin(i/dtorad)*360)
	npoints = npoints+dist_at_zen
print('look',1800/60/dtorad)
print(npoints)

#-------
#chi2 experiment



new_x = np.arange(-100,101,1)
new_f = np.zeros(201)
new_chi = np.zeros(201)

x = np.arange(0,1000,1)
y1 = 10*np.exp(-0.5*((x-500)/200)**2)

for i in range(201):
	y2 = 100*np.exp(-0.5*((x-400-i)/200)**2)+100
	f = (y1*(y2-100)/y2).sum()/(y1**2/y2).sum()
	chi = ((y2-100-f*y1)**2/(100+f*y1)).sum()
	new_f[i] = f
	new_chi[i] = chi
print('f\n',new_f)
plt.plot(new_x,new_f,label = 'f')
plt.plot(new_x,new_chi,label = 'chi')
min_index = np.argmin(new_chi)
min_new_x = new_x[min_index]
plt.axvline(x=min_new_x,color = 'r',label = str(min_new_x))
plt.legend()
plt.savefig(savedir + 'A_f_chi2.png')
plt.close()



#-------
plt.figure(figsize=(20,10))
ax = plt.axes(projection=ccrs.Mollweide(central_longitude=180),facecolor = '#f6f5ec')
ax.plot(ra,dec,color='blue', linewidth=2, marker='.', transform=ccrs.Geodetic(),alpha=0.2)


#----------------------------------------------------------------------

ax.set_global()
ax.invert_xaxis()
lats_y = np.arange(-75,76,15)
lats_x = np.zeros(lats_y.size)
lons_x = np.arange(0,360,30)
lons_y = np.zeros(lons_x.size)

#print(Nightshade(date, alpha=0.2))
xticks = list(range(-180, 180, 30))
yticks = list(range(-90, 90, 15))
ax.gridlines(xlocs=xticks, ylocs=yticks)
y = np.linspace(-90,90,100)
y_z = np.zeros(y.size)
x = np.linspace(0,360,100)
x_z = np.zeros(x.size)
proj_xyz = ax.projection.transform_points(ccrs.Geodetic(),np.array([0,30]), np.array([0,0]))
#ax.text(proj_xyz[0][0],proj_xyz[0][1],'240+15',size = 20,ha = 'right')
dx_ = np.abs(proj_xyz[0][0]-proj_xyz[1][0])
#print(proj_xyz)
lats_y_ticke = ax.projection.transform_points(ccrs.Geodetic(),lats_x+180+179.99, lats_y*1.1)
lats_y_x = lats_y_ticke[:,0]*0.86
lats_y_y = lats_y_ticke[:,1]
#print(lats_y_ticke)
#print(lats_y_x)
#print(lats_y_y)
for index,i in enumerate(lons_x):

	#ax.plot(y_z+i,y,'-',color = '#d3d7d4',linewidth = 1,transform = ccrs.Geodetic())
	ax.text(i,lons_y[index],r'$%d^{\circ}$'%i,transform = ccrs.Geodetic(),size = 20)

for index,i in enumerate(lats_y):

	#ax.plot(x,x_z+i,'-',color = '#d3d7d4',linewidth = 1,transform = ccrs.Geodetic())
	ax.text(lats_y_x[index]+dx_,lats_y_y[index],r'$%d^{\circ}$'%i,size = 20,ha = 'right',va = 'center')
plt.savefig(savedir+'A_map.png')
plt.close()


