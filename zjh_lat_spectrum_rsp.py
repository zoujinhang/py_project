from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

import sympy
from scipy import optimize
import Data_analysis.file as myfile
from Data_analysis.Baseline import WhittakerSmooth
import emcee
import corner
import json
import pymultinest
import os
from zjh_plotmarginalmodes import PlotMarginalModes


#rsp_link = '/home/laojin/my_lat/spectrum/response_matrix/glg_cspec_b0_bn150330828_v00.rsp'
rsp_link = '/home/laojin/my_lat/spectrum/response_matrix/glg_cspec_n2_bn150330828_v00.rsp'
savedir = '/home/laojin/my_lat/spectrum/response_matrix/'


#y_sp,y_sp_err = myfile.readcol(savedir +'ZZ_spectrum0.txt')
y_sp,y_sp_err = myfile.readcol(savedir +'ZZ_spectrum_n2_0.txt')



hl = fits.open(rsp_link)
matr = hl[2].data['MATRIX']#响应矩阵在这里
matr[-1] = np.zeros(128)
e_lo = hl[2].data['ENERG_LO']
e_hi = hl[2].data['ENERG_HI']
e_min = hl[1].data['E_MIN']
e_max = hl[1].data['E_MAX']
#print(matr)
#print(matr[0])
#print(len(matr))
#print(len(matr[0]))
a = np.ones(128)
a = np.mat(a).T
aa = np.vstack(matr)
R = np.mat(aa)
b = R*a
#print(R)
#print(np.array(b).T)
e_c0 = np.sqrt(e_min*e_max)
#y = e_c0**2*np.exp(-0.001*(e_c0-40)**2)
y = np.array(y_sp)/(e_max-e_min)
y_sp_err = np.array(y_sp_err)/np.sqrt(e_max-e_min)
#y = np.ones(128)+10
Y = np.mat(y)
Y_E = np.mat(y_sp_err)

def f(x):
	X = np.mat(x)
	return np.array(X*R)[0]
	
def log_xx(x):
	X = np.mat(x)
	M = (X*R-Y)/Y_E
	S = M*M.T
	return S[0,0]
def log_like(x):
	return -0.5*log_xx(x)
def log_fprime(x):
	X = np.mat(x)
	G = 2*((X*R-Y)/Y_E/Y_E)*R.T
	return np.array(G)[0]
def log_hess(x):
	return np.array(2*R*R.T)
e_c = np.sqrt(e_lo*e_hi)
x0 = np.zeros(140)
x_opt = optimize.fmin_ncg(log_xx,x0,fprime=log_fprime,fhess=log_hess)
print('---------------------------------')
print(x_opt)

plt.plot(e_c,x_opt)
plt.xscale('log')
plt.savefig(savedir+'A_fit_ncg.png')
plt.close()
plt.plot(e_c0,y)
plt.plot(e_c0,f(x_opt))
plt.xscale('log')
plt.savefig(savedir+'A_data_ncg.png')
plt.close()
bounds = [[0,+np.inf]]*140

xxx = optimize.fmin_l_bfgs_b(log_xx,x0,fprime = log_fprime,bounds = bounds,epsilon=1.e-10)[0]
#xxx = optimize.fmin_tnc(log_xx,x0,fprime = log_fprime,bounds = bounds)[0]
w = np.ones(140)
xxx1 = WhittakerSmooth(xxx,w,1)
print(bounds)
print('---------------------------------')
print(xxx)

plt.title('Unfolded spectrum')
plt.plot(e_c,xxx,'o',color = 'orange',label = 'Unfolded data',markersize = 2)
plt.plot(e_c,xxx1,'o',color = 'g',label = 'Unfolded data of smooth',markersize = 2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy (Kev)')
plt.ylabel('A(E)')
plt.legend()
plt.savefig(savedir+'A_fit_l_bfgs_b.png')
plt.close()

plt.title('Folded spectrum')
plt.errorbar(e_c0,y,yerr=y_sp_err,fmt = '.',elinewidth=2,capsize=2,alpha=0.3,label = 'spectrum data')
#plt.plot(e_c0,y)
plt.plot(e_c0,f(xxx),label = 'Folded data')
plt.plot(e_c0,f(xxx1),label = 'Folded data of smooth')
plt.xlabel('Energy (Kev)')
plt.ylabel('n /s/Kev')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig(savedir+'A_data_l_bfgs_b.png')
plt.close()

def log_prior(x):
	x_0 = x[x<-0.00001]
	if len(x_0)>0:
		return -np.inf
	else:
		return 0.0
def log_probability(x):
	lp = log_prior(x)
	if not np.isfinite(lp):
		return -np.inf
	return lp+log_like(x)

def band(e,k,a,b,Ec):
	val = Ec*(a-b)
	index1 = np.where(e<val)[0]
	e1 = e[index1]
	index2 = np.where(e>=val)[0]
	e2 = e[index2]
	a1 = k*(e1/100)**a*np.exp(-(e1/Ec))
	a2 = k*np.exp(b-a)*(val/100)**(a-b)*(e2/100)**b
	A = np.concatenate((a1,a2))
	return A

def band_model_log_like(cube,ndim,nparams):
	k = cube[0]
	a = cube[1]
	b = cube[2]
	Ec = cube[3]
	if a<b:
		return -np.inf
	spec = np.zeros(140)
	for i in range(140):
		e = np.linspace(e_lo[i],e_hi[i],5)
		rate = band(e,k,a,b,Ec)
		if True in np.isnan(rate):
			return -np.inf
		spec[i] = rate.mean()
	yu = f(spec)
	loglike = (-0.5*((yu[10:110] - y_sp[10:110])/y_sp_err[10:110])**2).sum()
	#loglike = log_like(spec)
	#print(loglike)
	return loglike

def prior(cube,ndim,nparams):
	cube[0] = 10000*cube[0]
	cube[1] = 30*cube[1]-20
	cube[2] = (cube[1]+30)*cube[2]-30
	cube[3] = 10**(20*cube[3]-10)
	
parameters = ['K','a','b','Ec']
n_params = len(parameters)
mul_dir = savedir+'out/'
if os.path.exists(mul_dir) ==False:
	os.makedirs(mul_dir)
pymultinest.run(band_model_log_like, prior, n_params, outputfiles_basename=mul_dir+'_1_',resume = True, verbose = True)
json.dump(parameters, open(mul_dir + 'params.json', 'w'))
a1 = pymultinest.Analyzer(outputfiles_basename=mul_dir + '_1_', n_params = n_params)
p = PlotMarginalModes(a1)
plt.figure(figsize=(5*n_params,5*n_params))
for i in range(n_params):
	for j in range(i,n_params):
		if(i == 0):
			plt.subplot(n_params, n_params,i*n_params+j+1)
			plt.title(parameters[j],size = 30)
			plt.tick_params(labelsize = 15)
			p.plot_marginal(j, with_ellipses = True, with_points = False, grid_points=50)
			#plt.errorbar(maximum[j],0,xerr=[[sigmal[j]],[sigmah[j]]],fmt = 'o',ecolor = 'r',capthick=2,color = 'r')
			if(j == 0):
				plt.ylabel("Probability",size = 30)
				plt.xlabel(parameters[j],size = 30)
			#plt.tight_layout()
		else:
			plt.subplot(n_params, n_params,i*n_params+j+1)
			plt.tick_params(labelsize = 15)
			p.plot_conditional(j, i-1, with_ellipses = False, with_points = False, grid_points=30)
			#plt.errorbar(maximum[j], maximum[i-1], xerr=[[sigmal[j]],[sigmah[j]]], yerr=[[sigmal[i-1]],[sigmah[i-1]]], fmt='o', ecolor='r', capthick=2, color='r')
			if(j == i):
				plt.xlabel(parameters[j],size = 30)
				plt.ylabel(parameters[i-1],size = 30)
				#plt.tight_layout()

plt.savefig(savedir + 'Z_marginals_multinest.png')
plt.close()
k0,a0,b0,Ec0 = a1.get_best_fit()['parameters']

for (k,a,b,Ec) in a1.get_equal_weighted_posterior()[::100,:-1]:
	plt.plot(e_c, band(e_c,k,a,b,Ec), '-', color='k', alpha=0.3, label='data')
print('best fit ',k0,a0,b0,Ec0)
plt.plot(e_c, band(e_c,k0,a0,b0,Ec0), '-', color='r')
plt.plot(e_c,xxx1,color = 'orange')
plt.xscale('log')
#plt.yscale('log')
#plt.ylim(10**-4,)
plt.savefig(savedir + 'A_multinest.png')
plt.close()
results = a1.get_stats()
print(results)

plt.plot(e_c0,y)
plt.plot(e_c0,f(xxx))
plt.plot(e_c0,f(xxx1))
plt.plot(e_c0,f(band(e_c,k0,a0,b0,Ec0)), color='r')
plt.xscale('log')
plt.savefig(savedir+'Z_data_multinest.png')
plt.close()

fig, ax = plt.subplots()
ax.errorbar(e_c0,y,yerr = y_sp_err)
fig.savefig(savedir+'Z_data_multinest2.png')
plt.close(fig)


'''
pos = xxx1 + 1e-1 * np.random.rand(300,140)
pos[pos<0]=0
print(pos)
nwalkers, ndim = pos.shape
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability)
sampler.run_mcmc(pos, 5000,progress=True)


samples = sampler.get_chain()
fig, axes = plt.subplots(10, figsize=(10, 20), sharex=True)
plt.subplots_adjust(left = 0.1,right = 0.9,top = 0.95,bottom = 0.05)
for i in range(ndim):
	ax = axes[i%10]
	ax.clear()
	ax.plot(samples[:, :, i], "k", alpha=0.3)
	ax.set_xlim(0, len(samples))
	if i%10 == 9:
		fig.savefig(savedir+'Z_emcee_'+str(i)+'.png')
fig.savefig(savedir+'Z_emcee_'+str(i)+'.png')
plt.close(fig)

flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
print(flat_samples.shape)
er1_list = []
er2_list = []
for i in range(ndim):
	mcmc = np.percentile(flat_samples[:, i], [16,84])
	er1,er2 = np.abs(mcmc-xxx1[i])
	er1_list.append(er1)
	er2_list.append(er2)

er = np.array([er1_list,er2_list])
print(er)
plt.errorbar(e_c,xxx1,yerr = er)
#plt.plot(e_c,xxx)
#plt.plot(e_c,xxx-er[0])
#plt.plot(e_c,xxx+er[1])
plt.xscale('log')
plt.yscale('log')
plt.savefig(savedir+'B_fit_l_bfgs_b.png')
plt.close()
'''

'''
x1, x2 = sympy.symbols("x_1, x_2")
f_sym = (x1-1)**4 + 5 * (x2-1)**2 - 2*x1*x2
fprime_sym = [f_sym.diff(x_) for x_ in (x1, x2)]

fhess_sym = [[f_sym.diff(x1_, x2_) for x1_ in (x1, x2)] for x2_ in (x1, x2)]
f_lmbda = sympy.lambdify((x1, x2), f_sym, 'numpy')
fprime_lmbda = sympy.lambdify((x1, x2), fprime_sym, 'numpy')
fhess_lmbda = sympy.lambdify((x1, x2), fhess_sym, 'numpy')
def func_XY_to_X_Y(f):
	return lambda X: np.array(f(X[0], X[1]))


f = func_XY_to_X_Y(f_lmbda)
fprime = func_XY_to_X_Y(fprime_lmbda)
fhess = func_XY_to_X_Y(fhess_lmbda)
print('---------------------')
print(f([1,1]))
print(fprime([1,1]))
print(fhess([1,1]))

'''

