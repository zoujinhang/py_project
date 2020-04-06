from matplotlib.patches import Polygon
import pyproj
from spherical_geometry.polygon import SphericalPolygon
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.patches import Polygon
from matplotlib import patches

savedir = '/home/laojin/my_lat/basemap/'
if os.path.exists(savedir) == False:
	os.makedirs(savedir)


poly = SphericalPolygon.from_cone(0,0,15,steps=50)
x,y = [p for p in poly.to_radec()][0]

poly1 = SphericalPolygon.from_cone(0,90,15,steps=50)
x1,y1 = [p for p in poly1.to_radec()][0]

xx = np.concatenate((x,x1))
yy = np.concatenate((y,y1))
print(x,y)



fig = plt.figure()
ax = fig.add_subplot(1,1,1)
e1 = patches.Ellipse((0, 0), 15,15,angle=0, linewidth=2, fill=False, zorder=2)


map = Basemap(projection='moll',lat_0=0,lon_0=60,resolution = 'l',area_thresh=1000.0,celestial=False,ax = ax)
xx,yy = map(xx,yy)
new_xy = list(zip(xx,yy))
print(new_xy)
poly1 = Polygon(new_xy,facecolor='coral',edgecolor='coral',linewidth=0)
ax.add_patch(poly1)
ax.add_patch(e1)
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
#map.fillcontinents(color='coral',lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
#map.drawmapboundary(fill_color='aqua')
# draw lat/lon grid lines every 30 degrees.

map.drawmeridians(np.arange(-180,180,30),dashes=[1,0],color='#d9d6c3')
map.drawparallels(np.arange(-90,90,30),dashes=[1,0], labels=[1,0,0,1], color='#d9d6c3')
# make up some data on a regular lat/lon grid.
nlats = 73; nlons = 145; delta = 2.*np.pi/(nlons-1)
lats = (0.5*np.pi-delta*np.indices((nlats,nlons))[0,:,:])
lons = (delta*np.indices((nlats,nlons))[1,:,:])
wave = 0.75*(np.sin(2.*lats)**8*np.cos(4.*lons))
mean = 0.5*np.cos(2.*lats)*((np.sin(2.*lats))**2 + 2.)
# compute native map projection coordinates of lat/lon grid.
#x, y = map(lons*180./np.pi, lats*180./np.pi)
# contour data over the map.
#cs = map.contour(x,y,wave+mean,15,linewidths=1.5)
plt.title('contour lines over filled continent background')
plt.savefig(savedir + 'A_basemap.png')
plt.close()

