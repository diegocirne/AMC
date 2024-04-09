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



def stationary_sigma(mu):


	arg = mu*mu + math.sqrt(mu**4+4*mu**2)


	return math.sqrt(arg/2)



b_lis = [0.75, 0.50, 0.25, 0]
b_lis = [0, 0.25, 0.50, 0.75]



lin = ['dotted', 'dashdot', 'dashed', 'solid']
lin = ['solid', 'dashed', 'dashdot', 'dotted']



mu_inf = 0
mu_sup = 2
mu_div = pow(10, 3)



mu_lis = [(mu_sup-mu_inf)/mu_div*i + mu_inf for i in range(mu_div)]



for i in range(len(b_lis)):


	b = b_lis[i]


	sta_sigma = []
	
	
	mu_lis_b = [(1-b)*mu for mu in mu_lis]
	
	
	for j in range(len(mu_lis_b)):


		mu = mu_lis_b[j]


		sta_sigma.append(stationary_sigma(mu))


	lab = sep.join(('$b=', str(b), '$'))


	plt.plot(mu_lis, sta_sigma, color='black', label=lab, linestyle=lin[i], marker='o', linewidth=3, markersize=0)


	asy = 1-b


	asy_lis = [(1-b)*mu for mu in mu_lis]


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



def stationary_lag(sigma):


	arg = 1 + 1/sigma**2


	return arg



mu_lis = [(mu_sup-mu_inf)/mu_div*i + mu_inf for i in range(1, mu_div)]



for i in range(len(b_lis)):


	b = b_lis[i]


	mu_lis_b = [(1-b)*mu for mu in mu_lis]


	sta_lag = []


	for j in range(len(mu_lis_b)):


		mu = mu_lis_b[j]

		
		sta_sigma = stationary_sigma(mu)
		sta_lag.append((1-b)*stationary_lag(sta_sigma))


	lab = sep.join(('$b=', str(b), '$'))


	plt.plot(mu_lis, sta_lag, color='black', label=lab, linestyle=lin[i], marker='o', linewidth=3, markersize=0)


	asy = 1-b


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



def critical_B(sigma, Wmax):


	arg = Wmax/math.sqrt(sigma**2+1)


	return sigma**2/math.sqrt(sigma**2+1)*math.sqrt(2*math.log(arg))



Wmax = 2.0



mu_inf = pow(10, -3)
mu_inc = pow(10, -3)



for i in range(len(b_lis)):


	b = b_lis[i]


	if stationary_sigma((1-b)*mu_inf) < math.sqrt(Wmax*Wmax-1):


		mu_lis = [mu_inf]


		while stationary_sigma((1-b)*(mu_lis[-1]+mu_inc)) < math.sqrt(Wmax*Wmax-1):


			mu_lis.append(mu_lis[-1]+mu_inc)


		Bc = [critical_B(stationary_sigma((1-b)*mu), Wmax)/(1-b) for mu in mu_lis]
		
		
		Bc[-1] = 0


		lab = sep.join(('$b=', str(b), '$'))


		plt.plot(mu_lis, Bc, color='black', label=lab, linestyle=lin[i], marker='o', linewidth=3, markersize=0)



gap = 1.0



plt.plot([mu_inf, mu_lis[-1]+gap], [0, 0], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.xlim([mu_inf, mu_lis[-1]+gap])
#plt.ylim()
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\delta/\omega$', size=LAB, labelpad=PAD)
plt.ylabel('$B_c/\omega$', rotation=90, size=LAB, labelpad=PAD)
plt.title('$C$', size=TIT)
plt.legend(loc='best', fontsize=LAB)
#plt.legend(shadow=False, frameon=False, prop={'size': LAB})
plt.xticks(size=TIC)
plt.yticks(size=TIC)



fig.set_size_inches(22, 7, forward=True)
fig.savefig('figure-4.pdf', bbox_inches='tight')




#################### Bottom ####################
