from astropy.time import Time
import numpy as np



class Time_triansform():
	def __init__(self,time_origin):
		leap_second = [41499.0,41499.0,41864.0,42413.0,42778.0,43144.0,43509.0,43874.0,44239.0,44786.0,45151.0,45516.0,46247.0,47161.0,47892.0,48257.0,48804.0,49169.0,49534.0,50083.0,50630.0,51179.0,53736.0,54832.0,56109.0,57204.0,57754.0]
		add_second = np.arange(len(leap_second))+1
		index_r = np.argsort(-add_second)
		self.leap_second = leap_second[index_r]
		self.add_second = add_second[index_r]
		self.time_origin = Time(time_origin).mjd
	def utc_to_met(self,utc):
		mjd_ = Time(utc).mjd
		add_s = 0
		for index,value in enumerate(self.leap_second):
			if mjd_ >= value:
				add_s = self.add_second[index]
				break
		dt = (mjd_ - self.time_origin)*86400
		dt  = dt + add_s
		return dt






