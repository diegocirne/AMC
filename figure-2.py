# 23/11/02 #



#################### Library ####################



import math
import matplotlib.pyplot as plt
import scipy.optimize as sci



#################### Form ####################



PAD = 10
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



h_lis = [0.7, 0.8, 0.9, 1.0]
h_lis = [1.0, 0.9, 0.8, 0.7]



mu_inf = 0
mu_sup = 2
mu_div = pow(10, 3)



mu_lis = [(mu_sup-mu_inf)/mu_div*i + mu_inf for i in range(mu_div)]



lin = ['dotted', 'dashdot', 'dashed', 'solid']
lin = ['solid', 'dashed', 'dashdot', 'dotted']



for i in range(len(h_lis)):


	h = h_lis[i]


	sta_sigma = []


	for j in range(len(mu_lis)):


		mu = mu_lis[j]


		sta_sigma.append(sci.fsolve(stationary_sigma_equation, x0=1, args=(h, mu)))


	lab = sep.join(('$h=', str(h), '$'))


	plt.plot(mu_lis, sta_sigma, color='black', label=lab, linestyle=lin[i], marker='o', linewidth=3, markersize=0)


	asy = 1/math.sqrt(1-(1-h**2)**2)


	asy_lis = [asy*mu for mu in mu_lis]


	plt.plot(mu_lis, asy_lis, color='red', linestyle=lin[i], marker='o', linewidth=3, markersize=0)



plt.plot([mu_inf, mu_sup], [0, 0], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.xlim([mu_inf, mu_sup])
#plt.ylim()
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\delta/\omega$', size=LAB, labelpad=PAD)
plt.ylabel('$\sigma_*/\omega}$', rotation=90, size=LAB, labelpad=PAD)
plt.title('$A$', size=TIT)
plt.legend(loc='best', fontsize=LAB)
#plt.legend(shadow=False, frameon=False, prop={'size': LAB})
plt.xticks(size=TIC)
plt.yticks(size=TIC)



plt.subplot(1, 3, 2)



def stationary_lag(sigma, h):


	arg = 1 + 1/sigma**2


	return arg/h**2



mu_lis = [(mu_sup-mu_inf)/mu_div*i + mu_inf for i in range(1, mu_div)]



for i in range(len(h_lis)):


	h = h_lis[i]


	sta_lag = []


	for j in range(len(mu_lis)):


		mu = mu_lis[j]


		sta_sigma = sci.fsolve(stationary_sigma_equation, x0=1, args=(h, mu))
		sta_lag.append(stationary_lag(sta_sigma, h))


	lab = sep.join(('$h=', str(h), '$'))


	plt.plot(mu_lis, sta_lag, color='black', label=lab, linestyle=lin[i], marker='o', linewidth=3, markersize=0)


	asy = 1/h**2


	plt.plot([mu_inf, mu_sup], [asy, asy], color='red', linestyle=lin[i], marker='o', linewidth=3, markersize=0)


plt.plot([mu_inf, mu_sup], [0, 0], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



y_inf = -0.25
y_sup = 5



plt.xlim([mu_inf, mu_sup])
plt.ylim([y_inf, y_sup])
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\delta/\omega$', size=LAB, labelpad=PAD)
plt.ylabel('$\Delta_*/B$', rotation=90, size=LAB, labelpad=PAD)
plt.title('$B$', size=TIT)
plt.legend(loc='best', fontsize=LAB)
#plt.legend(shadow=False, frameon=False, prop={'size': LAB})
plt.xticks(size=TIC)
plt.yticks(size=TIC)



plt.subplot(1, 3, 3)



def critical_B(sigma, Wmax, h):


	arg = Wmax/math.sqrt(sigma*sigma+1)


	return h**2*sigma**2/math.sqrt(sigma**2+1)*math.sqrt(2*math.log(arg))



Wmax = 2.0



mu_inf = pow(10, -3)
mu_inc = pow(10, -3)



for i in range(len(h_lis)):


	h = h_lis[i]


	mu_lis = []
	sta_sigma_lis = []


	mu = mu_inf
	sta_sigma = sci.fsolve(stationary_sigma_equation, x0=1, args=(h, mu_inf))


	while sta_sigma < math.sqrt(Wmax*Wmax-1):


		mu_lis.append(mu)
		sta_sigma_lis.append(sta_sigma)


		mu += mu_inc
		sta_sigma = sci.fsolve(stationary_sigma_equation, x0=1, args=(h, mu))


	Bc = [critical_B(sta_sigma, Wmax, h) for sta_sigma in sta_sigma_lis]


	Bc[-1] = 0


	lab = sep.join(('$h=', str(h), '$'))


	plt.plot(mu_lis, Bc, color='black', label=lab, linestyle=lin[i], marker='o', linewidth=3, markersize=0)



gap = 1.0



plt.plot([mu_inf, mu_lis[-1]+gap], [0, 0], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.xlim([mu_inf, mu_lis[-1]+gap])
#plt.ylim()
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\delta/\omega$', size=LAB, labelpad=PAD)
plt.ylabel('$B_c\,/\omega$', rotation=90, size=LAB, labelpad=PAD)
plt.title('$C$', size=TIT)
plt.legend(loc='best', fontsize=LAB)
#plt.legend(shadow=False, frameon=False, prop={'size': LAB})
plt.xticks(size=TIC)
plt.yticks(size=TIC)



fig.set_size_inches(22, 7, forward=True)
fig.savefig('figure-2.pdf', bbox_inches='tight')




#################### Bottom ####################
