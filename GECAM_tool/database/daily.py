


from astropy.io import fits
from astropy.time import Time
from ..file import findfile
import numpy as np
import pandas as pd
from ..clock import Clock

def time_slic(time_start,time_stop,format = None,scale='utc',H = 3):

	if isinstance(time_start,Time) ==False:
		time_start = Time(time_start, format=format, scale=scale)
	if isinstance(time_stop,Time) ==False:
		time_stop =Time(time_stop,format=format,scale=scale)
	re = pd.date_range(time_start.fits,time_stop.fits,freq = str(H)+'H')
	re = pd.DataFrame({'time':re})
	stp = pd.date_range(time_stop.fits,freq = 'H',periods=1)
	stp = pd.DataFrame({'time':stp})

	t_ar = pd.concat([re,stp])
	t_ar.drop_duplicates('time','first',inplace=True,ignore_index=True)
	t_ar.sort_values(by = 'time',inplace=True,ignore_index=True)
	time_list = t_ar['time'].values
	t_ar1 = pd.DataFrame({'time_start':time_list[:-1],'time_stop':time_list[1:]})
	return t_ar1.values


class Daily_database(object):

	def __init__(self,topdir,clock = None):

		self.clock = clock
		self.topdir = topdir
		if self.clock is None:
			self.clock = Clock()

	def get_satellite_position(self,time_start,time_stop,format = None,scale='utc'):


		names = ['A','B']
		dir_list=['GECAM_A','GECAM_B']
		col_name= ['time','q1','q2','q3','q4','x','y','z','xw','yw','zw']

		if isinstance(time_start,Time) == False:
			time_start = Time(time_start, format=format, scale=scale)
		if isinstance(time_stop,Time) == False:
			time_stop =Time(time_stop,format=format,scale=scale)

		met_start = self.clock.utc_to_met(time_start)
		met_stop = self.clock.utc_to_met(time_stop)
		date_time_arr = pd.date_range(time_start.fits,time_stop.fits,freq = 'H')

		return_value = {}
		for i,dir_ in enumerate(dir_list):
			q1, q2, q3, q4 = [], [], [], []
			x, y, z = [], [], []
			xw, yw, zw = [], [], []
			time = []
			for t_i in range(date_time_arr.shape[0]):
				date_t = date_time_arr[t_i]
				year = '%d' % date_t.year
				month = '%.2d' % date_t.month
				day = '%.2d' % date_t.day
				hour = '%.2d' % date_t.hour
				link = self.topdir + year + '/' + month + '/' + day + '/'+dir_+'/posatt/'
				fm = '_posatt_' + year[-2:] + month + day + '_'+hour+ '_v'
				name_list= findfile(link,fm)
				if len(name_list)>0:
					hl = fits.open(link+name_list[0])
					data = hl[1].data
					time_ = data.field(0)
					q1_,q2_,q3_,q4_ = data.field(1),data.field(2),data.field(3),data.field(4)
					x_,y_,z_ = data.field(9),data.field(10),data.field(11)
					xw_,yw_,zw_ = data.field(15),data.field(16),data.field(17)
					time.append(time_)
					q1.append(q1_),q2.append(q2_),q3.append(q3_),q4.append(q4_)
					x.append(x_),y.append(y_),z.append(z_)
					xw.append(xw_),yw.append(yw_),zw.append(zw_)
					hl.close()
				else:
					print('warning! the file under line is lost!')
					print(fm)
			if len(time)>0:
				time = np.concatenate(time)
				q1,q2,q3,q4 = np.concatenate(q1),np.concatenate(q2),np.concatenate(q3),np.concatenate(q4)
				x,y,z = np.concatenate(x),np.concatenate(y),np.concatenate(z)
				xw,yw,zw = np.concatenate(xw),np.concatenate(yw),np.concatenate(zw)
				inde_sort = np.argsort(time)
				new_data = np.vstack([time,q1,q2,q3,q4,x,y,z,xw,yw,zw]).T
				new_data = new_data[inde_sort]
				_,uni_index = np.unique(new_data[:,0],return_index=True)
				new_data = new_data[uni_index]
				index_ = np.where((new_data[:,0]>=met_start-2)&(new_data[:,0]<=met_stop+2))[0]
				new_data = new_data[index_]
				pd_data = pd.DataFrame(data=new_data,dtype=np.float128,columns=col_name)
				return_value[names[i]] = pd_data
			else:
				return_value[names[i]] = None

		return return_value






