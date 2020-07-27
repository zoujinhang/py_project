from zjh_download import download_all_in_one_path
import pandas as pd

databasetop = ''
timestart = '2020-03-01'
timestop = '2020-06-30'

def time_list(time_start,time_stop):
	
	list = pd.date_range(time_start,time_stop,freq = 'D')
	returnlist = []
	
	for date in list:
		year = '%d' % date.year
		month = '%.2d' % date.month
		day = '%.2d' % date.day
		list = year+'/'+month+'/'+day +'/'
		returnlist.append(list)
	return returnlist

datelist = time_list(timestart,timestop)

for date in datelist:
	targetdir = '/fermi/data/gbm/daily/' + date + 'current/'
	resultdir = databasetop + date
	download_all_in_one_path(targetdir,resultdir,num = 15)
	




