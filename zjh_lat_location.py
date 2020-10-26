import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import Data_analysis.file as myfile
import cartopy.crs as ccrs


filelink = '/home/laojin/my_lat/location/locrates_1deg_50_300_hard_n.dat'

savedir = '/home/laojin/my_lat/location/'
if os.path.exists(savedir) ==False:
	os.makedirs(savedir)


c = np.fromfile(filelink,dtype=np.float)
dat_ = myfile.readcol(filelink)
print(c[0])
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


