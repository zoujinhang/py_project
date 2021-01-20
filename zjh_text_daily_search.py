
from daily_search import Database,time_slic
from daily_search import Search
import matplotlib
matplotlib.use('Agg')
import os
from astropy.coordinates import SkyCoord
from multiprocessing import Pool
from daily_search.file import readcol
from daily_search.file import printdatatofile
datatop = '/media/laojin/My_Passport/Fermi_GBM_Dailiy_Data/'
savetop = '/home/laojin/my_lat/daily_search/new_code12/'
if os.path.exists(savetop) ==False:
	os.makedirs(savetop)


n_max = 3
t_start = '2020-04-01T00:00:00'
t_stop = '2020-04-30T23:59:00'

t_ar = time_slic(t_start,t_stop,H=3)
print('The time interval, H =',3)
print(t_ar)
n = len(t_ar)
print('the number of time window:',n)
if n > n_max:
	n = n_max
print('multiprocessing number:',n)
source = SkyCoord(ra=[293.729,128.426],dec= [21.3864,27.712],frame ='icrs',unit='deg')      #you can change these things
name = ['SGRJ1935','sn2020adow']



def run(t_ar):

	t_start, t_stop = t_ar
	time_markker = str(t_start)[:13]+'_3H'
	fermi_data = Database(datatop)
	data = fermi_data.get_detector_data(t_start,t_stop)
	#print('data:',data['n0']['events'])
	pd_pos_data = fermi_data.get_poshist_data(t_start,t_stop)
	print(time_markker+' searching begin!')
	errs = None
	errs_list =[]
	try:
		search = Search(data,pd_pos_data,marker = time_markker)
	except:
		errs = 'search candidate error!'
	else:
		s_n = search.get_candidate_number()
		print(time_markker + ' searching over! The number of the candidates is ',s_n)
		if s_n >=1:
			save_candidate = savetop + 'all_candidates/'
			if os.path.exists(save_candidate) == False:
				os.makedirs(save_candidate)

			#search.save_triggers(save_candidate+'A_trig_summary.csv')
			save_perspec = search.get_triggers()
			utc_time = save_perspec['utc'].values

			#print('-----------------------------------')
			for i in range(s_n):
				utc_i = utc_time[i]
				year = utc_i[0:4]
				month = utc_i[5:7]
				date = utc_i[8:10]
				hour = utc_i[11:13]
				minute = utc_i[14:16]
				second = utc_i[17:19]
				dotb = utc_i[20:22]
				time_markker1 = 'tg_' + year + month + date + '_' + hour + minute + second + dotb
				savedir = save_candidate + time_markker1 + '/'
				if os.path.exists(savedir) == False:
					os.makedirs(savedir)
				search.save_candidate_poshist(i,savedir,time_markker1)
				search.save_candidate_TTE(i,savedir,time_markker1)
				search.save_candidate_signal(i,savedir+'Z_candidate.csv')
				search.save_candidate_5_50_signal(i,savedir + 'Z_candidate_5_50.csv')
				search.save_candidate_50_300_signal(i, savedir + 'Z_candidate_50_300.csv')
				search.plot_candidate_5_50_lc(i, savedir + 'X_lc_5_50.png')
				search.plot_candidate_50_300_lc(i, savedir + 'X_lc_50_300.png')
				search.plot_candidate_lc(i, savedir + 'X_lc.png')
				search.plot_candidate_5_50_SNR(i, savedir + 'Y_snr_5_50.png')
				search.plot_candidate_50_300_SNR(i, savedir + 'Y_snr_50_300.png')
				search.plot_candidate_SNR(i, savedir + 'Y_snr.png')
				try:
					search.plot_sky_map(i, savedir + 'W_skymap.png')
				except:
					errs_list.append(str(i) + ' plot sky map error!')
				try:
					search.save_location_50_300_chi2(i,savedir + 'V_location_chi2_50_300.csv')
				except:
					errs_list.append(str(i) + ' save location 50-300 chi2 error!')
				try:
					search.save_location_5_50_chi2(i, savedir + 'V_location_chi2_5_50.csv')
				except:
					errs_list.append(str(i) + ' save location 50-300 chi2 error!')

			for i in range(len(source)):

				try:
					track = search.track_one(source[i], name=name[i])
				except:
					errs_list.append(name[i] + ' track error!')
				else:
					n = track.get_candidate_number()
					print(time_markker + ' track the source of', name[i], '. The number of candidates is '+str(n))
					if n >=1:

						savedir = savetop + name[i]+'/'
						if os.path.exists(savedir) == False:
							os.makedirs(savedir)

						#track.save_triggers(savedir + 'A_trig_summary.csv')
						track_perspec = track.get_triggers()
						utc_time = track_perspec['utc'].values
						track_perspec.to_csv(savedir + 'A_trig_summary_'+ time_markker+'.csv')
						for i in range(n):
							utc_i = utc_time[i]
							year = utc_i[0:4]
							month = utc_i[5:7]
							date = utc_i[8:10]
							hour = utc_i[11:13]
							minute = utc_i[14:16]
							second = utc_i[17:19]
							dotb = utc_i[20:22]
							time_markker2 = 'tg_' + year + month + date + '_' + hour + minute + second + dotb
							savecandidate = savedir + time_markker2 +'/'
							if os.path.exists(savecandidate)==False:
								os.makedirs(savecandidate)
							track.save_candidate_5_50_signal(i,savecandidate+'Z_candidate_5_50.csv')
							track.save_candidate_50_300_signal(i,savecandidate+'Z_candidate_50_300.csv')
							track.save_candidate_signal(i,savecandidate+'Z_candidate.csv')

							track.plot_candidate_5_50_lc(i,savecandidate + 'X_lc_5_50.png')
							track.plot_candidate_50_300_lc(i,savecandidate + 'X_lc_50_300.png')
							track.plot_candidate_lc(i,savecandidate + 'X_lc.png')

							track.plot_candidate_5_50_SNR(i,savecandidate + 'Y_snr_5_50.png')
							track.plot_candidate_50_300_SNR(i,savecandidate + 'Y_snr_50_300.png')
							track.plot_candidate_SNR(i,savecandidate + 'Y_snr.png')
							try:
								track.plot_sky_map(i,savecandidate + 'W_skymap.png')
							except:
								errs_list.append(name[i] + ' ' + str(i) + ' plot sky map error!')
	if errs is not None:
		print('something wrong of ' + time_markker)
		savedir = savetop + 'wrong/'
		if os.path.exists(savedir) == False:
			os.makedirs(savedir)
		printdatatofile(savedir + 'A_wrong_' + time_markker + '.txt', data=[[errs]])
	if len(errs_list) > 0:
		print('something wrong of ' + time_markker)
		savedir = savetop + 'wrong/'
		if os.path.exists(savedir) == False:
			os.makedirs(savedir)
		printdatatofile(savedir + 'A_wrong_' + time_markker + '.txt', data=[errs_list])


pool = Pool(n)     # Parallel computing
pool.map(run,list(t_ar))
pool.close()
pool.join()




