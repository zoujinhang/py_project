
import numpy as np




def get_energy_of_ch(time,e1,e2):
	'''

	:param time: time
	:param e1:
	:param e2:
	:return:
	'''
	numb = len(time)
	energy_random_arr = np.random.random_sample(numb)
	energy_array = e1 + (e2-e1)*energy_random_arr
	return energy_array

def ch_to_energy(time,ch,ch_n,e1,e2):
	'''

	:param time:
	:param ch:
	:param ch_n:
	:param e1:
	:param e2:
	:return:
	'''
	new_t = []
	new_energy = []
	new_ch = []
	for index,channel in enumerate(ch_n):
		ch_t_index = np.where(ch == channel)
		ch_t = time[ch_t_index]
		ch_n1 = ch[ch_t_index]
		energy_array = get_energy_of_ch(ch_t,e1[index],e2[index])
		new_t.append(ch_t)
		new_energy.append(energy_array)
		new_ch.append(ch_n1)
	new_t = np.concatenate(new_t)
	new_energy = np.concatenate(new_energy)
	new_ch = np.concatenate(new_ch)
	index_all = np.argsort(new_t)
	new_t = new_t[index_all]
	new_energy = new_energy[index_all]
	new_ch = new_ch[index_all]
	return new_t,new_energy,new_ch


def overlap_bins(range_ ,binsize = 1.0 ,stepsize = 0.5):

	start = np.arange(range_[0] ,range_[-1]-binsize ,stepsize)
	stop = start + binsize

	return np.vstack([start ,stop]).T

def overlap_curve(t ,bins_arr):

	n = np.zeros(bins_arr.shape[0])
	t = np.array(t)
	for index,edges in enumerate(bins_arr):
		nn_ = np.where((t>=edges[0])&(t<edges[-1]))[0]
		n[index] = nn_.size
	return n


