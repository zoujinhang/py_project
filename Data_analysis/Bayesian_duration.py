import numpy as np
import matplotlib.pyplot as plt


def background_correction(t,rate,edges,degree = 50,plot_save = None):
	'''
	
	:param t:
	:param rate:
	:param edges:
	:param degree:
	:return:
	'''
	re_rate,re_sigma,index_list = re_histogram(t,rate,edges)
	dt = t[1]-t[0]
	binsize = edges[1:]-edges[:-1]
	sort_index = np.argsort(-binsize)
	sort_binsize = binsize[sort_index]
	sort_sigma = re_sigma[sort_index]
	sort_re_rate = re_rate[sort_index]
	sort_index_list = np.array(index_list)[sort_index]
	background_pool = rate[sort_index_list[0]]
	correction_t = [edges[0],edges[-1]+dt,t[sort_index_list[0]][0],t[sort_index_list[0]][-1]]
	mean1 = sort_re_rate[0]
	correction_rate = [rate[0],rate[-1],mean1,mean1]
	sigma1 = sort_sigma[0]
	binsize1 = sort_binsize[0]
	n = 1
	for i in range(1,len(sort_binsize)):
		if confidence_analysis(mean1,sigma1,binsize1,sort_re_rate[i],sort_sigma[i],sort_binsize[i],dt,degree):
			n = n+1
			correction_t.append(t[sort_index_list[i]][0])
			correction_t.append(t[sort_index_list[i]][-1])
			correction_rate.append(sort_re_rate[i])
			correction_rate.append(sort_re_rate[i])
			background_pool = np.concatenate((background_pool,rate[sort_index_list[i]]))
			mean1 = background_pool.mean()
			sigma1 = background_pool.std(ddof = 1)
			binsize1 = binsize1 + sort_binsize[i]
	print('The number of background blocks is %d.'%n)
	correction_t = np.array(correction_t)
	correction_rate = np.array(correction_rate)
	sort_t_index = np.argsort(correction_t)
	correction_t = correction_t[sort_t_index]
	correction_rate = correction_rate[sort_t_index]
	if plot_save is not None:
	
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax.step(edges,np.concatenate((re_rate[:1],re_rate)))
		ax.plot(correction_t,correction_rate)
		fig.savefig(plot_save)
		plt.close(fig)
	
	correction_b = np.interp(t,correction_t,correction_rate)
	new_rate = rate - correction_b + correction_b.mean()
	re_rate,re_sigma,index_list = re_histogram(t,new_rate,edges)
	result = {'lc':(t,new_rate),
	          're_hist':(re_rate,re_sigma,index_list),
	          'bkg':(mean1,sigma1,binsize1),
	          'edges':edges,
	          'dt':dt
	          }
	return result

def get_SNR(edges,re_rate,background_mean,background_sigma,dt,non_negative = True):
	'''
	
	:param edges:
	:param re_rate:
	:param background_mean:
	:param background_sigma:
	:param dt:
	:param non_negative:
	:return:
	'''
	binsize = edges[1:] - edges[:-1]
	sigma = background_sigma * np.sqrt(dt/binsize)
	SNR = (re_rate - background_mean)/sigma
	if non_negative:
		SNR[SNR<0] = 0 #we only care about the part where the SNR is greater than 0.
	return SNR

def re_histogram(t,rate,edges):
	'''
	This is a stable function, which can re-bin the data 't' and 'rate' with the new bin edges 'edges'.
	
	:param t: The time of light curve.
	:param rate: The counts rate of light curve.
	:param edges: the new bin-edges, whose interval can different from each other.
	:return: three array ,The rebined rate and sigma of each bins.
	'''
	
	index = np.where(t >= edges[0])[0]
	t = t[index]
	rate = rate[index]
	re_rate = np.zeros(len(edges)-1)
	re_sigma = np.zeros(len(edges)-1)
	edges_num = 1
	index_list = []
	ones = []
	for index1,value in enumerate(t):
		if value > edges[edges_num] or value == t[-1]:
			ones.append(index1)
			index_list.append(ones)
			ones = []
			edges_num = edges_num+1
			if edges_num == len(edges):
				break
		else:
			ones.append(index1)
	for index1,value in enumerate(index_list):
		one_rate = rate[value]
		index_list[index1] = np.array(value)+index[0]
		re_rate[index1] = one_rate.mean()
		if len(one_rate) == 1:
			re_sigma[index1] = one_rate.std(ddof = 0)
		else:
			re_sigma[index1] = one_rate.std(ddof = 1)
	return re_rate,re_sigma,np.array(index_list)

def confidence_analysis(mean1,sigma1,T1,mean2,sigma2,T2,dt,degree = 3):
	
	'''
	conditional dependence of confidence
	:param mean1:
	:param sigma1:
	:param T1:
	:param mean2:
	:param sigma2:
	:param T2:
	:param dt:
	:param degree:
	:return: True or False
	'''
	
	sigma1_2 = sigma1*np.sqrt(dt/T2)
	sigma2_1 = sigma2*np.sqrt(dt/T1)
	return (mean1>=mean2-degree*sigma2_1)&(mean1<=mean2+degree*sigma2_1)&(mean2>=mean1-degree*sigma1_2)&(mean2<=mean1+degree*sigma1_2)




def get_bayesian_duration(data,sigma = 5):
	'''
	
	:param data:
	:param sigma:
	:return:
	'''
	start = 0
	stop = 0
	pulse = 1
	cafe = 2
	fringe_up = 3
	fringe_down = 4
	#t,rate = data['lc']
	edges = data['edges']
	re_rate,re_sigma,index_list = data['re_hist']
	dt = data['dt']
	bkg,bkgsigma,bkgsize = data['bkg']
	SNR = get_SNR(edges,re_rate,bkg,bkgsigma,dt,non_negative=True)
	
	binstart = edges[:-1]
	#binstop = edges[1:]
	trait = []
	for index1,hight in enumerate(re_rate):
		if (index1 == 0):
			trait.append(start)
		elif (index1 == len(re_rate)-1):
			trait.append(stop)
		else:
			if (hight> re_rate[index1-1]) and (hight > re_rate[index1+1]):
				trait.append(pulse)
			elif (hight <= re_rate[index1-1]) and (hight <= re_rate[index1+1]):
				trait.append(cafe)
			elif (hight < re_rate[index1-1]) and (hight >= re_rate[index1+1]):
				trait.append(fringe_down)
			elif (hight >= re_rate[index1-1]) and (hight < re_rate[index1+1]):
				trait.append(fringe_up)
			else:
				trait.append(cafe)
	start_tag = False
	start_edges = []
	stop_edges = []
	
	for index,value in enumerate(trait):

		if (SNR[index]>sigma) and (start_tag == False):
			
			if value != fringe_down:
				start_edges.append(binstart[index])
				start_tag = True
			else:
				stop_edges.pop()
				start_tag = True
				
			
		elif (SNR[index] <= sigma) and start_tag:
			
			if value != fringe_up:
				stop_edges.append(binstart[index])
				start_tag = False
	if start_tag:
		start_edges.pop()
	if start_edges[0] == binstart[0]:
		start_edges = start_edges[1:]
		stop_edges = stop_edges[1:]
	return np.array(start_edges),np.array(stop_edges)

	



