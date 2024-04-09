# 23/11/02 #



#################### Library ####################



import math
import matplotlib.pyplot as plt
import scipy.optimize as sci



#################### Form ####################



PAD = 5
TIC = 15
LAB = 20
TIT = 30



sep = ''



#################### Main ####################



fig = plt.figure()



plt.subplot(1, 3, 1)



def stationary_sigma_equation(sigma, *par):


	h, mu = par


	arg = h**2*sigma/math.sqrt(1+sigma**2) + (1-h**2)*sigma


	return arg**2 + mu**2 - sigma**2



def S_map(sigma):


	arg = sigma/math.sqrt(1+sigma**2)


	return arg



def H_map(S, h):


	arg = h**2*S + (1-h**2)*S/math.sqrt(1-S**2)


	return arg



def S_derivative(sigma):


	arg = 1/pow(1+sigma**2, 3/2)


	return arg



def H_derivative(S, h):


	arg = h**2 + (1-h**2)/pow(1-S**2, 3/2)


	return arg



def M_derivative(H, mu):


	return H/math.sqrt(H**2+mu**2)



h = 0.9



b_lis = [0.75, 0.50, 0.25, 0]
b_lis = [0, 0.25, 0.50, 0.75]



lin = ['dotted', 'dashdot', 'dashed', 'solid']
lin = ['solid', 'dashed', 'dashdot', 'dotted']



mu_inf = 0
mu_sup = 3
mu_div = pow(10, 3)



mu_lis = [(mu_sup-mu_inf)/mu_div*i + mu_inf for i in range(1, mu_div)]



for i in range(len(b_lis)):


	b = b_lis[i]


	mu_lis_b = [(1-b)*mu for mu in mu_lis]


	d_sigma = []


	for j in range(len(mu_lis_b)):


		mu = mu_lis_b[j]


		sta_sigma = sci.fsolve(stationary_sigma_equation, x0=1, args=(h, mu))


		d_sigma.append(abs(M_derivative(H_map(S_map(sta_sigma), h), mu)*H_derivative(S_map(sta_sigma), h)*S_derivative(sta_sigma)))


	lab = sep.join(('$b=', str(b), '$'))


	plt.plot(mu_lis, d_sigma, color='black', label=lab, linestyle=lin[i], marker='o', linewidth=3, markersize=0)


	asy = (1-h**2)**2


	plt.plot([mu_inf, mu_sup], [asy, asy], color='red', linestyle=lin[i], marker='o', linewidth=3, markersize=0)



plt.plot([mu_inf, mu_sup], [0, 0], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.plot([mu_inf, mu_sup], [1, 1], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.xlim([mu_inf, mu_sup])
#plt.ylim()
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\delta/\omega$', size=LAB, labelpad=PAD)
plt.ylabel('$\left|d\sigma_{t+1}/\,d\sigma_t\\right|_{\sigma_*}$', rotation=90, size=LAB, labelpad=PAD)
#plt.ylabel('$\left|\\frac{d\sigma_{t+1}}{d\sigma_t}\\right|_{\sigma_*}$', rotation=90, size=LAB, labelpad=PAD)
plt.title('$A$', size=TIT)
plt.legend(loc='best', fontsize=LAB)
#plt.legend(shadow=False, frameon=False, prop={'size': LAB})
plt.xticks(size=TIC)
plt.yticks(size=TIC)



plt.subplot(1, 3, 2)



def lag_derivative(sigma, h):


	arg = 1 - h**2*sigma**2/(1+sigma**2)


	return arg



for i in range(len(b_lis)):


	b = b_lis[i]

	mu_lis_b = [(1-b)*mu for mu in mu_lis]


	d_lag = []


	for j in range(len(mu_lis_b)):


		mu = mu_lis_b[j]


		sta_sigma = sci.fsolve(stationary_sigma_equation, x0=1, args=(h, mu))


		d_lag.append(abs(lag_derivative(sta_sigma, h)))


	lab = sep.join(('$b=', str(b), '$'))


	plt.plot(mu_lis, d_lag, color='black', label=lab, linestyle=lin[i], marker='o', linewidth=2, markersize=0)


	asy = 1-h**2


	plt.plot([mu_inf, mu_sup], [asy, asy], color='red', linestyle=lin[i], marker='o', linewidth=3, markersize=0)



plt.plot([mu_inf, mu_sup], [0, 0], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.plot([mu_inf, mu_sup], [1, 1], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.xlim([mu_inf, mu_sup])
#plt.ylim([y_inf, y_sup])
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\delta/\omega$', size=LAB, labelpad=PAD)
plt.ylabel('$\left|d\Delta_{t+1}/\,d\Delta_t\\right|_{(\Delta_*, \sigma_*)}$', rotation=90, size=LAB, labelpad=PAD)
#plt.ylabel('$\left|\\frac{d\Delta_{n+1}}{d\Delta_n}\\right|_{(\Delta_*, \sigma_*)}$', rotation=90, size=LAB, labelpad=PAD)
plt.title('$B$', size=TIT)
plt.legend(loc='best', fontsize=LAB)
#plt.legend(shadow=False, frameon=False, prop={'size': LAB})
plt.xticks(size=TIC)
plt.yticks(size=TIC)



plt.subplot(1, 3, 3)



for i in range(len(b_lis)):


	b = b_lis[i]


	ratio_lis = []


	mu_lis_b = [(1-b)*mu for mu in mu_lis]



	for j in range(len(mu_lis_b)):


		mu = mu_lis_b[j]


		sta_sigma = sci.fsolve(stationary_sigma_equation, x0=1, args=(h, mu))


		ratio_lis.append(H_map(S_map(sta_sigma), h)/sta_sigma)


	lab = sep.join(('$b=', str(b), '$'))


	plt.plot(mu_lis, ratio_lis, color='black', label=lab, linestyle=lin[i], marker='o', linewidth=3, markersize=0)


	asy = 1-h**2


	plt.plot([mu_inf, mu_sup], [asy, asy], color='red', linestyle=lin[i], marker='o', linewidth=3, markersize=0)



plt.plot([mu_inf, mu_sup], [0, 0], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.plot([mu_inf, mu_sup], [1, 1], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.xlim([mu_inf, mu_sup])
#plt.ylim()
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\delta/\omega$', size=LAB, labelpad=PAD)
plt.ylabel('$\sigma_*^H/\sigma_*$', rotation=90, size=LAB, labelpad=PAD)
plt.title('$C$', size=TIT)
plt.legend(loc='best', fontsize=LAB)
#plt.legend(shadow=False, frameon=False, prop={'size': LAB})
plt.xticks(size=TIC)
plt.yticks(size=TIC)



fig.set_size_inches(22, 7, forward=True)
fig.savefig('figure-5.pdf', bbox_inches='tight')



#################### Bottom ####################
