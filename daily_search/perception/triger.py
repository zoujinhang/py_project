
import numpy as np
from ..background import TD_baseline
from .lightcurve import Event
from astropy.stats import bayesian_blocks
from .bayesian_duration import background_correction,get_bayesian_duration


def bayesian_trig(data,name,windowlist,detector):

	'''
	bayesian blocks trigger.

	data: obtain from database.
	windowlist: list :[[start,stop],[start,stop],...]

	return: edges_list,new_window_list,name_list

	'''

	#sigma = 3
	#strong = 7
	#pf_strong = 10
	#wt = 0.1
	#binsize_baseline = 1
	binszie_lightcurve = 0.05
	#energy_band = [[5,50],[50,300],[300,900]]
	#name = data.keys()
	new_window_list  = []
	new_edges_list = []
	new_name_list = []
	lc_window_index_list = []
	for index0,(start,stop) in enumerate(windowlist):
		during = stop - start
		#print('lc durtion',during)
		if during >= 10*binszie_lightcurve:
			#bins_baseline = np.arange(start, stop, binsize_baseline)
			#bins_baseline_c = 0.5 * (bins_baseline[1:] + bins_baseline[:-1])
			#bb_cc = np.concatenate(([start], bins_baseline_c, [stop]))
			bins_lightcurve = np.arange(start, stop, binszie_lightcurve)
			#bins_lightcurve_c = 0.5 * (bins_lightcurve[1:] + bins_lightcurve[:-1])
			new_start = []
			new_stop = []
			trig_name = []
			for detei in name:
				ni = data[detei]
				ch_E = ni['ch_E']
				t = ni['events']['TIME'].values
				ch = ni['events']['PHA'].values
				lightcurve = Event(ch_E, t, ch)
				#rate_b = lightcurve(bins = bins_baseline,energy_band=[8,105],channel=True)[1]
				#rate_b = np.concatenate((rate_b[:1],rate_b,rate_b[-1:]))
				#bs_f_5_900 = get_background_f(bb_cc,rate_b)
				t_cc,rate_c = lightcurve(bins = bins_lightcurve,energy_band=[8,105],channel=True)
				#bs5_900 = bs_f_5_900(t_cc)
				bs5_900 = TD_baseline(t_cc,rate_c)
				n_c5_900 = np.round(rate_c*binszie_lightcurve)
				rate_c_n = rate_c-bs5_900 + bs5_900.mean()
				index = n_c5_900>0
				try:
					edges = bayesian_blocks(t_cc[index],n_c5_900[index],fitness='events',gamma = np.exp(-4))
				except:
					continue
				edges = np.round(edges / binszie_lightcurve) * binszie_lightcurve
				edges = np.unique(edges)
				edges = np.sort(edges)

				if len(edges)>=4:
					result = background_correction(t_cc,rate_c_n,edges,degree = 7.5)
					startedges,stopedges = get_bayesian_duration(result,sigma = 4,max_snr=False)
					if startedges.size == stopedges.size:
						if startedges.size > 0:
							if startedges.size == 1:
								during1 = stopedges - startedges
								if during1[0] >= 10:
									edges_index = np.where(
										(edges >= startedges[0]) & (
												edges <= stopedges[0]))[
										0]
									if edges_index.size <= 2:
										print('edges_index.size <= 2')
										continue
									elif edges_index.size > 2:
										ne_ed = edges[edges_index]
										d_ed = ne_ed[1:] - ne_ed[:-1]
										indei2 = np.where(d_ed <= 3)[0]
										if indei2.size == 0:
											print('indei2.size == 0')
											continue

								t_index = \
								np.where(t_cc < startedges[0])[0]
								rate_c1 = rate_c[t_index].mean()
								t_index = np.where(t_cc > stopedges[0])[0]
								rate_c2 = rate_c[t_index].mean()
								if abs(rate_c1 - rate_c2) > 600:
									print('abs(rate_c1 - rate_c2) > 600')
									continue

							if startedges.size > 1:
								startedges, stopedges = new_start_stop(startedges,stopedges,cafe=0.3)
							new_start.append(startedges)
							new_stop.append(stopedges)
							trig_name.append([detei]*len(startedges))
						elif during >= 6:

							edges_during = edges[1:] - edges[:-1]
							if edges_during[0] > 1:
								edges[0] = edges[0] + 1
							else:
								edges = edges[1:]
							if edges_during[-1] >1:
								edges[-1] = edges[0] - 1
							else:
								edges = edges[:-1]
							edges = np.round(edges / binszie_lightcurve) * binszie_lightcurve
							edges = np.unique(edges)
							edges = np.sort(edges)
							if len(edges) >= 4:
								#bins_baseline1 = np.arange(start+1, stop-1, binsize_baseline)
								#bins_baseline_c1 = 0.5 * (
								#		bins_baseline1[1:] + bins_baseline1[:-1])
								#bb_cc1 = np.concatenate(
								#	([start+1], bins_baseline_c1, [stop-1]))
								bins_lightcurve1 = np.arange(start+1, stop-1, binszie_lightcurve)
								#rate_b = lightcurve(bins = bins_baseline1,energy_band=[8,105],channel=True)[1]
								#rate_b = np.concatenate((rate_b[:1],rate_b,rate_b[-1:]))
								#bs_f_5_900 = get_background_f(bb_cc1,rate_b)
								t_cc,rate_c = lightcurve(bins = bins_lightcurve1,energy_band=[8,105],channel=True)
								#bs5_900 = bs_f_5_900(t_cc)
								bs5_900 = TD_baseline(t_cc,rate_c)
								rate_c_n = rate_c-bs5_900 + bs5_900.mean()
								result = background_correction(t_cc,rate_c_n,edges,degree = 7.5)
								startedges,stopedges = get_bayesian_duration(result,sigma = 4,max_snr=False)

								if startedges.size == stopedges.size:
									if startedges.size > 0:
										if startedges.size > 1:
											startedges, stopedges = new_start_stop(
												startedges, stopedges,
												cafe=0.5)
										new_start.append(startedges)
										new_stop.append(stopedges)
										trig_name.append([detei] * len(startedges))
					elif during >= 6:
						edges_during = edges[1:] - edges[:-1]
						if edges_during[0] > 1:
							edges[0] = edges[0] + 1
						else:
							edges = edges[1:]
						if edges_during[-1] >1:
							edges[-1] = edges[0] - 1
						else:
							edges = edges[:-1]
						edges = np.round(edges / binszie_lightcurve) * binszie_lightcurve
						edges = np.unique(edges)
						edges = np.sort(edges)
						if len(edges) >= 4:
							#bins_baseline1 = np.arange(start+1, stop-1, binsize_baseline)
							#bins_baseline_c1 = 0.5 * (
							#	bins_baseline1[1:] + bins_baseline1[:-1])
							#bb_cc1 = np.concatenate(
							#	([start + 1], bins_baseline_c1, [stop - 1]))
							bins_lightcurve1 = np.arange(start+1, stop-1, binszie_lightcurve)
							#rate_b = lightcurve(bins = bins_baseline1,energy_band=[8,105],channel=True)[1]
							#rate_b = np.concatenate((rate_b[:1],rate_b,rate_b[-1:]))
							#bs_f_5_900 = get_background_f(bb_cc1,rate_b)
							t_cc,rate_c = lightcurve(bins = bins_lightcurve1,energy_band=[8,105],channel=True)
							#bs5_900 = bs_f_5_900(t_cc)
							bs5_900 = TD_baseline(t_cc,rate_c)
							rate_c_n = rate_c-bs5_900 + bs5_900.mean()
							result = background_correction(t_cc,rate_c_n,edges,degree = 7.5)
							startedges,stopedges = get_bayesian_duration(result,sigma = 4,max_snr=False)

							if startedges.size == stopedges.size:
								if startedges.size > 0:
									if startedges.size > 1:
										startedges, stopedges = new_start_stop(
											startedges, stopedges, cafe=0.5)
									new_start.append(startedges)
									new_stop.append(stopedges)
									trig_name.append([detei] * len(startedges))
			#print('bayesian over!')
			if len(new_start) >= 2: # the trigger of bayesian blocks is good when trigger number is 2.
				time_edges,namelist = time_overlap(new_start,new_stop,trig_name,binszie_lightcurve)
				new_time_edges = []
				new_namelist = []
				lc_window_index = []
				for index_1,edges1 in enumerate(time_edges):
					if len(namelist[index_1])==2:
						#print('bayes_namelist',namelist[index_1])
						if detector.detector_association(namelist[index_1],n = 2,m = 2):
							new_time_edges.append(edges1)
							new_namelist.append(namelist[index_1])
							lc_window_index.append(index0)
					elif len(namelist[index_1])>2:
						if detector.detector_association(namelist[index_1],n = 3,m = 3):
							new_time_edges.append(edges1)
							new_namelist.append(namelist[index_1])
							lc_window_index.append(index0)
					#else:
					#	new_time_edges.append(edges1)
					#	new_namelist.append(namelist[index_1])
				new_window_list = new_window_list + get_windows(new_time_edges,start,stop,dt = 0,no_overlap=False)
				new_edges_list = new_edges_list + new_time_edges
				new_name_list = new_name_list + new_namelist
				lc_window_index_list = lc_window_index_list + lc_window_index

	return new_edges_list,new_window_list,new_name_list,lc_window_index_list

def new_start_stop(start,stop,cafe = 1.0):

	a = start[1:]
	b = stop[:-1]
	dt = a-b
	index_ = np.where(dt>cafe)[0]
	a = a[index_]
	b = b[index_]
	return np.concatenate((start[:1],a )),np.concatenate((b,stop[-1:]))


def time_overlap(start,stop,namelist,dt):

	'''
	this function will combine the periods overlap with others.

	start : list : [array(),array(),...]
	stop : list : [array(),array(),...]
	namelist: list :[list(),list(),...]

	return : time list[[t_start,t_stop],[t_start,t_stop],...], name list[list(),list(),...]
	'''
	all_start = np.concatenate(start)
	all_stop = np.concatenate(stop)
	#all_namelist = np.concatenate(namelist)
	smi_t = np.arange(all_start.min()-3,all_stop.max()+3,dt)
	zeo = np.zeros((smi_t.size,len(start)),dtype=int)
	name_all = []
	for indexi,namelist in enumerate(namelist):
		name_all.append(namelist[0])
		start_ar = start[indexi]
		stop_ar = stop[indexi]
		for i in range(start_ar.size):
			ind = np.where((smi_t>=start_ar[i])&(smi_t<=stop_ar[i]))[0]
			zeo[ind,indexi] = 1
	name_all = np.array(name_all)
	new_list = []
	new_namelist = []
	find_start = True
	t_start = all_start.min()
	name = None
	for indexi,zeoi in enumerate(zeo):
		ind = np.where(zeoi == 1)[0]
		if name is not None:
			if ind.size>len(name):
				name = name_all[ind]
		if find_start and ind.size>=2:
			t_start = smi_t[indexi]
			name = name_all[ind]
			find_start = False
		elif ind.size < 2 and find_start == False:
			new_list.append([t_start,smi_t[indexi]])
			new_namelist.append(name)
			find_start = True

	'''
	index_ = np.argsort(all_start)
	all_start = all_start[index_]
	all_stop = all_stop[index_]
	#all_during = all_stop-all_start
	all_namelist = all_namelist[index_]
	#print('all_namelist',all_namelist)
	all_dt = (all_start[1:]-all_start[:-1])
	one_list_start = [all_start[0]]
	one_list_stop = [all_stop[0]]
	one_list_name = [all_namelist[0]]
	#one_list_during = [all_during[0]]
	new_list = []
	new_namelist = []
	for i in range(len(all_dt)):

		if all_dt[i]<=dt*2:
			one_list_start.append(all_start[i+1])
			one_list_stop.append(all_stop[i+1])
			one_list_name.append(all_namelist[i+1])
			#one_list_during.append(all_during[i+1])
		else:
			trun_n = 0
			for j in range(len(one_list_start)):
				if all_start[i+1] < one_list_stop[j] and all_stop[i+1] > one_list_start[j] :#  and one_list_during[j] >= 2 and all_during[i+1] >=2:
					trun_n = trun_n+1
			if len(one_list_start) == trun_n:
				one_list_start.append(all_start[i+1])
				one_list_stop.append(all_stop[i+1])
				one_list_name.append(all_namelist[i+1])
				#one_list_during.append(all_during[i+1])
			else:
				if len(one_list_start) == 1:
					t_start = one_list_start[0]
					t_stop  = one_list_stop[0]
				elif len(one_list_start) == 2:
					dt_ = abs(one_list_start[1] - one_list_start[0])
					if dt_<=dt*3:
						t_start = np.mean(one_list_start)
						t_stop = max(one_list_stop)
					else:
						t_start = max(one_list_start)
						t_stop = max(one_list_stop)
				else:
					mean_start = np.mean(one_list_start)
					one_arr_start = np.array(one_list_start)
					std_ = (one_arr_start).std()
					while std_> 3*dt and one_arr_start.size >= 2:
						dt_ = np.abs(one_arr_start-mean_start)
						index_dt = np.where(dt_<dt_.max())[0]
						one_arr_start = one_arr_start[index_dt]
						mean_start = (one_arr_start).mean()
						std_ = (one_arr_start).std()
					t_start = mean_start
					t_stop = max(one_list_stop)
				new_list.append([t_start,t_stop])
				new_namelist.append(one_list_name)
				one_list_start = [all_start[i+1]]
				one_list_stop = [all_stop[i+1]]
				one_list_name = [all_namelist[i+1]]


	if len(one_list_start) == 1:
		t_start = one_list_start[0]
		t_stop  = one_list_stop[0]
	elif len(one_list_start) == 2:
		dt_ = abs(one_list_start[1] - one_list_start[0])
		if dt_<=dt*3:
			t_start = np.mean(one_list_start)
			t_stop = max(one_list_stop)
		else:
			t_start = max(one_list_start)
			t_stop = max(one_list_stop)
	else:
		mean_start = np.mean(one_list_start)
		one_arr_start = np.array(one_list_start)
		std_ = (one_arr_start).std()
		while std_> 3*dt and one_arr_start.size >= 2:
			dt_ = np.abs(one_arr_start-mean_start)
			index_dt = np.where(dt_<dt_.max())[0]
			one_arr_start = one_arr_start[index_dt]
			mean_start = (one_arr_start).mean()
			std_ = (one_arr_start).std()
		t_start = mean_start
		t_stop = max(one_list_stop)
	new_list.append([t_start,t_stop])
	#print('new_list',new_list)
	new_namelist.append(one_list_name)
	#print('new_namelist',new_namelist)
	'''

	return new_list,new_namelist

def check_edges(edges):

	new_edges = []
	for (start,stop) in edges:
		try:
			if stop-start >0:
				new_edges.append([start,stop])
		except:
			print('')
			continue
	return new_edges

def try_to_trig(data,name,detector):

	'''

	try to find the window edges.

	data: obtain from database.
	detector: the class Detector()
	return : the window of candidate : list[[start,stop],[start,stop],...]

	'''
	sigma = 3.0       #---the frist standerd.
	#strong = 7
	#pf_strong = 10
	#binsize_baseline = 1
	binszie_lightcurve = 0.05
	#name = data.keys()
	name = np.array(list(name))
	edges = get_bins(data,name)
	#print('get edges:',edges)
	window_list = []
	for start,stop in edges:
		during = stop-start
		if during < 2:
			continue
		#print(start,stop)
		#bins_baseline = np.arange(start,stop,binsize_baseline)
		#bins_baseline_c = 0.5*(bins_baseline[1:]+bins_baseline[:-1])
		#bb_cc = np.concatenate(([start],bins_baseline_c,[stop]))
		bins_lightcurve = np.arange(start,stop,binszie_lightcurve)
		bins_lightcurve_c = 0.5*(bins_lightcurve[1:]+bins_lightcurve[:-1])

		SNR5_900_list = []
		for detei in name:

			ni = data[detei]
			ch_E = ni['ch_E']
			t = ni['events']['TIME'].values
			ch = ni['events']['PHA'].values
			lightcurve = Event(ch_E,t,ch)
			#rate_b = lightcurve(bins = bins_baseline,energy_band=[8,105],channel=True)[1]
			#rate_b = np.concatenate((rate_b[:1],rate_b,rate_b[-1:]))
			#bs_f_5_900 = get_background_f(bb_cc,rate_b)
			t_cc,rate_c = lightcurve(bins = bins_lightcurve,energy_band=[8,105],channel=True)
			#bs5_900 = bs_f_5_900(t_cc)
			bs5_900 = TD_baseline(t_cc,rate_c)
			scale = np.sqrt(bs5_900/binszie_lightcurve)
			SNR5_900_list.append((rate_c-bs5_900)/scale)

		SNR5_900_list = np.vstack(SNR5_900_list).T
		#print('SNR5_900_list:',SNR5_900_list)
		time_index_list = []

		for index,SNR_i in enumerate(SNR5_900_list):
			good_index = np.where(SNR_i>=sigma)[0]
			#strong_index = np.where(SNR_i>strong)[0]
			#pf_strong_index = np.where(SNR_i>pf_strong)[0]
			#nnn = len(good_index)
			if len(good_index)>=3:
				#time_index_list.append(index)
				namelist = name[good_index]
				if detector.detector_association(namelist,n = 3,m = 3): #---the frist standerd.
					time_index_list.append(index)

		if len(time_index_list)>0:
			time_index_list = np.array(time_index_list)
			time_index_list = get_subsection_index(time_index_list,binszie_lightcurve,distinguish = 10)
			lc_t_list = []
			for index_i in time_index_list:

				ti_ = bins_lightcurve_c[index_i]
				lc_t_list.append([ti_.min(),ti_.max()])
			window_list = window_list + get_windows(lc_t_list,start,stop,dt = 1)
	return window_list

def get_windows(lc_t_list,start,stop,dt = 0,no_overlap = True):

	'''
	expand the time from candidate .
	'''

	window_list = []
	for ind,lc_ti in enumerate(lc_t_list):

		lc_dt = lc_ti[-1]-lc_ti[0]

		if lc_dt<=2:
			add_t = 10
		elif lc_dt<=20:
			add_t = 20
		elif lc_dt<=50:
			add_t = 30
		elif lc_dt<=80:
			add_t = 50
		elif lc_dt<=100:
			add_t = 60
		else:
			add_t = 80
		range_t_min = lc_ti[0]-add_t
		if ind == 0:
			if range_t_min < start+dt:
				range_t_min = start+dt
		else:
			if no_overlap and range_t_min < lc_t_list[ind-1][-1]:
				range_t_min = lc_t_list[ind-1][-1]
			elif no_overlap == False and (range_t_min < lc_t_list[ind-1][-1] and range_t_min >= lc_t_list[ind-1][0]):
				range_t_min = lc_t_list[ind-1][0]-5
				if range_t_min < start+dt:
					range_t_min = start+dt
				if ind >= 2:
					if range_t_min < lc_t_list[ind-2][-1]:
						range_t_min = lc_t_list[ind-2][-1]



		range_t_max = lc_ti[-1] + add_t
		if ind == len(lc_t_list)-1:
			if range_t_max >stop-dt:
				range_t_max = stop-dt
		else:
			if no_overlap and range_t_max >lc_t_list[ind+1][0]:
				range_t_max = lc_t_list[ind+1][0]
			elif no_overlap == False and (range_t_max >lc_t_list[ind+1][0] and range_t_max <= lc_t_list[ind+1][-1]):
				range_t_max = lc_t_list[ind+1][-1]+5
				if range_t_max >stop-dt:
					range_t_max = stop-dt
				if ind <= len(lc_t_list)-1-2:
					if range_t_max >lc_t_list[ind+2][0]:
						range_t_max = lc_t_list[ind + 2][0]

		window_list.append([range_t_min,range_t_max])
	return window_list

def get_subsection_index(index,binsize,distinguish=1.1):
	'''

	:param index:
	:param binsize:
	:param distinguish:
	:return:
	'''

	if len(index) <= 1:
		return [index]
	else:
		index_list = []
		index_one = []
		d_index = index[1:]-index[:-1]
		d_index = np.concatenate(([0],d_index))
		for id,value in enumerate(d_index):
			if value*binsize < distinguish:
				index_one.append(index[id])
			else:
				index_list.append(np.array(index_one))
				index_one = [index[id]]
		index_list.append(np.array(index_one))
		return index_list



def get_bins(data,name,wt = 0.1):
	'''
	try to find the edges for every detector data

	data: obtain from database.
	return: the list of edges.

	'''
	#name = data.keys()
	n_dete = len(name)
	t_start = []
	t_stop = []

	for dete in name:
		t_edges = []
		ni = data[dete]
		t = ni['events']['TIME'].values
		if t.size == 0:
			return []
		indexlist = check_event(t,wt = wt)
		for index in indexlist:
			t_i = t[index]
			t_edges.append([t_i.min(),t_i.max()])

		t_edges = np.array(t_edges).T
		t_start.append(t_edges[0])
		t_stop.append(t_edges[1])


	all_t_start = np.ceil(np.concatenate(t_start) / wt) * wt
	all_t_stop = np.floor(np.concatenate(t_stop) / wt) * wt
	t_sim = np.arange(all_t_start.min()-2,all_t_stop.max()+2,wt)
	zer = np.zeros(t_sim.size,dtype=int)
	for i in range(n_dete):
		t_starti = t_start[i]
		t_stopi = t_stop[i]
		for j in range(t_starti.size):
			indexi = np.where((t_sim>=t_starti[j])&(t_sim<=t_stopi[j]))[0]
			zer[indexi] = zer[indexi]+1
	good_index = np.where(zer==n_dete)[0]
	index_list = get_subsection_index(good_index,1,distinguish = 2)
	bins_list = []
	for index in index_list:
		t_simi = t_sim[index]
		bins_list.append([np.ceil(t_simi.min() / wt) * wt + 1.0,np.floor(t_simi.max()/ wt) * wt - 0.5])
	return bins_list
'''
	t_start = np.ceil(np.concatenate(t_start) / wt) * wt
	t_stop = np.floor(np.concatenate(t_stop) / wt) * wt
	t_start = np.unique(t_start)
	start_c = np.zeros(len(t_start))+1
	t_stop = np.unique(t_stop)
	stop_c = np.zeros(len(t_stop))-1

	edges = np.concatenate([t_start,t_stop])
	c_ = np.concatenate([start_c,stop_c])
	edges_index = np.argsort(edges)
	edges = edges[edges_index]
	c_ = c_[edges_index]
	t_start = []
	t_stop = []
	start = edges[0]
	fand = 0
	for i in range(len(edges)):

		if c_[i] == 1:
			start = edges[i]
			fand = 1
		elif c_[i] == -1 and fand == 1:
			fand = 0
			stop = edges[i]
			t_start.append(start)
			t_stop.append(stop)
	bins_list = []
	n_name = len(name)
	for i in range(len(t_start)):
		n_ = 0
		for dete in name:
			ni = data[dete]
			t = ni['events']['TIME'].values
			t = t[np.where((t>=t_start[i])&(t<=t_stop[i]))[0]]
			indexlist = check_event(t,wt = wt)
			if len(indexlist)==1:
				n_ = n_ + 1
				#t_ = t[indexlist[0]]
				#if t_.min()<=t_start[i] and t_.max() >= t_stop[i]:
				#	n_ = n_ + 1
		#print('n_',n_)
		if n_ == n_name and t_stop[i]-t_start[i] > 0:

			bins_list.append([t_start[i],t_stop[i]])

	return bins_list
'''

def check_event(t,wt = 0.064):
	#print('len t',len(t))
	dt = t[1:]-t[:-1]
	dt = np.concatenate((dt[:1],dt))
	#index_d = np.where(dt>=wt)[0]
	#if len(index_d)>0
	#print('len dt',len(dt))
	#print('dt',dt)
	index_list = []
	index_ = []
	not_first = False
	#print('all:',dt.size)
	for index,dti in enumerate(dt):

		if dti < wt:
			not_first = True
			index_.append(index)
		else:
			#print('chack!')
			if not_first:
				#print('add',len(index_))
				index_list.append(np.array(index_))
				index_ = [index]

	#print('add',len(index_))
	index_list.append(np.array(index_))
	return index_list