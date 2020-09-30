
import re

def transpose(data):
	da = [[row[col] for row in data] for col in range(len(data[0]))]
	return da


def readcol(file_name):
	'''
	read txt
	:param file_name:
	:return:
	'''
	c = []
	f = open(file_name,'r')
	for line in f:
		a = line
		b = [i for i in re.split('[\t,\n,\s]',a) if i != '']
		#print(b)
		c = c + [b]
	f.close()
	c = [i for i in c if i != []]
	data_all = transpose(c)
	nl = len(data_all)
	for k in range(nl):
		for l in range(len(data_all[k])):
			try:
				data_all[k][l] = int(data_all[k][l])
			except ValueError:
				try:
					data_all[k][l] = float(data_all[k][l])
				except ValueError:
					data_all[k][l] = data_all[k][l]
	return data_all