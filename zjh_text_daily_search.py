

from daily_search import Database
from daily_search import Search
import matplotlib.pyplot as plt
import os



datatop = '/media/laojin/TOSHIBA_EXT/daily/'
savetop = '/home/laojin/my_lat/daily_search/'

t_start = '2020-04-28T00:00:00'
t_stop = '2020-04-28T01:00:00'

fermi_data = Database(datatop)

data = fermi_data.get_detector_data(t_start,t_stop)
#print('data:',data['n0']['events'])
pd_pos_data = fermi_data.get_poshist_data(t_start,t_stop)

search = Search(data,pd_pos_data)

print(search.candidate)



