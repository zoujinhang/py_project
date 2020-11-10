
import numpy as np
from astropy.io import fits
import pandas as pd
from ..file import findfile


class Burst_database(object):

	def __init__(self,topdir):

		self.topdir = topdir


	def get_satellite_position(self,sample_name):

		'''

		sample_name: sample_name
		return : {
			'A':pandas format:'time','q1','q2','q3','q4','x','y','z','xw','yw','zw',
			'B':pandas format:'time','q1','q2','q3','q4','x','y','z','xw','yw','zw'
			}
		'''
		dir_list = [self.topdir + sample_name + '/GECAM_A/',self.topdir + sample_name + '/GECAM_B/']
		names = ['A','B']
		fm = '_posatt_' + sample_name + '_v'
		col_name= ['time','q1','q2','q3','q4','x','y','z','xw','yw','zw']
		return_value = {}
		for i,dir_ in enumerate(dir_list):

			name_list = findfile(dir_,fm)
			if len(name_list) >0 :
				hl = fits.open(dir_+name_list[0])
				data = hl[1].data
				time = data.field(0)
				q1,q2,q3,q4 = data.field(1),data.field(2),data.field(3),data.field(4)
				x,y,z = data.field(9),data.field(10),data.field(11)
				xw,yw,zw = data.field(15),data.field(16),data.field(17)
				hl.close()
				inde_sort = np.argsort(time)
				new_data = np.vstack([time,q1,q2,q3,q4,x,y,z,xw,yw,zw]).T
				new_data = new_data[inde_sort]
				_,uni_index = np.unique(new_data[:,0],return_index=True)
				new_data = new_data[uni_index]
				pd_data = pd.DataFrame(data=new_data,dtype=np.float128,columns=col_name)
				return_value[names[i]] = pd_data
			else:
				return_value[names[i]] = None
				print('file lost!!')
		return  return_value







