
import numpy as np
import matplotlib.pyplot as plt
import pymultinest
from pymultinest import PlotMarginalModes
import os


savedir = '/home/laojin/my_lat/multinest_toul/'
if os.path.exists(savedir) ==False:
	os.makedirs(savedir)


x = np.arange(0,100,1)
y = 6*x+1 + np.random.randn(100)*5

plt.plot(x,y,'.')
plt.savefig(savedir + 'A_data.png')
plt.close()
parameters = ['a','b']
n_params = len(parameters)
def model(x,a,b):
	y = a*x+b
	return y
def prior(cube,ndim,nparams):
	cube[0] = 10*cube[0]
	cube[1] = 5*cube[1]-2.5

def loglike(cube,ndim,nparams):
	a = cube[0]
	b = cube[1]
	m_y = a*x+b
	loglike_ = -0.5*(((y-m_y)/5)**2).sum()
	return loglike_

pymultinest.run(loglike, prior, n_params, outputfiles_basename=savedir + 'A_model_', resume = False, verbose = False)
a_ = pymultinest.Analyzer(outputfiles_basename=savedir + 'A_model_', n_params = n_params)
b_ = a_.get_stats()['modes'][0]

print('mean',b_['mean'])
print('sigma',b_['sigma'])
print('maximum',b_['maximum'])
plt.plot(x,y,'.')
for (a,b) in a_.get_equal_weighted_posterior()[::100,:-1]:
	plt.plot(x, model(x, a, b), '-', color='blue', alpha=0.3, label='data')
plt.savefig(savedir+ 'B_1_posterior.png')
plt.close()

p = PlotMarginalModes(a_)
plt.figure(figsize=(5*n_params,5*n_params))
for i in range(n_params):
	for j in range(i,n_params):
		if(i == 0):
			plt.subplot(n_params, n_params,i*n_params+j+1)
			plt.title(parameters[j],size = 30)
			plt.tick_params(labelsize = 15)
			p.plot_marginal(j, with_ellipses = True, with_points = False, grid_points=50)
			#plt.errorbar(maximum[j],0,xerr=[[sigmal[j]],[sigmah[j]]],xlolims=lolims[j],xuplims=uplims[j],
			#	markersize=10,elinewidth = 3,
			#            fmt = 'o',ecolor = 'r',capthick=5,color = 'r')
			if(j == 0):
				plt.ylabel("Probability",size = 30)
				plt.xlabel(parameters[j],size = 30)
					#plt.tight_layout()
		else:
			plt.subplot(n_params, n_params,i*n_params+j+1)
			plt.tick_params(labelsize = 15)
			p.plot_conditional(j, i-1, with_ellipses = False, with_points = False, grid_points=30)
			#plt.errorbar(maximum[j], maximum[i-1], xerr=[[sigmal[j]],[sigmah[j]]], yerr=[[sigmal[i-1]],[sigmah[i-1]]],
			#            xlolims=lolims[j],xuplims=uplims[j],uplims=uplims[i-1],lolims=lolims[i-1],
			#             markersize=10,elinewidth = 3,
			#             fmt='o', ecolor='r', capthick=5, color='r')
			if(j == i):
				plt.xlabel(parameters[j],size = 30)
				plt.ylabel(parameters[i-1],size = 30)
plt.savefig(savedir + 'B_coner.png')
plt.close()



