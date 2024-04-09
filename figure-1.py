# 23/11/01 #



#################### Library ####################



import matplotlib.pyplot as plt
import math



#################### Form ####################



PAD = 10
TIC = 15
LAB = 20
TIT = 30



sep = ''



#################### Main ####################



fig = plt.figure()



plt.subplot(1, 3, 1)



def mean_fitness(sigma, lag):


	arg1 = math.sqrt(1+sigma**2)


	arg2 = lag**2/(1+sigma**2)


	return 1/arg1*math.e**(-arg2/2)



sigma_inf = pow(10, -2)
sigma_sup = pow(10, 2)
sigma_div = pow(10, 5)



lag_lis = [0, 1.0, 2.0, 3.0]



lin = ['solid', 'dashed', 'dashdot', 'dotted']



sigma_lis = [(sigma_sup-sigma_inf)/sigma_div*sigma+sigma_inf for sigma in range(sigma_div)]



for i in range(len(lag_lis)):


	lag = lag_lis[i]


	mean_W_lis = [mean_fitness(j, lag) for j in sigma_lis]


	lab = sep.join(('$\Delta_t\,/\omega=', str(lag), '$'))


	plt.plot(sigma_lis, mean_W_lis, color='black', label=lab, linestyle=lin[i], marker='o', linewidth=3, markersize=0)



plt.plot([sigma_inf, sigma_sup], [0, 0], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.plot([sigma_inf, sigma_sup], [1, 1], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.xlim([sigma_inf, sigma_sup])
#plt.ylim()
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\sigma_t\,/\omega$', size=LAB, labelpad=PAD)
plt.ylabel('$\overline{W}_t\,/W_{max}$', rotation=90, size=LAB, labelpad=PAD)
#plt.ylabel('$\\frac{\overline{W}_t}{W_{max}}$', rotation=0, size=25, labelpad=PAD)
plt.title('$A$', size=TIT)
plt.legend(loc='best', fontsize=LAB)
#plt.legend(shadow=False, frameon=False, prop={'size': LAB})
plt.xticks(size=TIC)
plt.yticks(size=TIC)



plt.subplot(1, 3, 2)



def estationary_mean_W(sigma, B):


	arg1 = math.sqrt(1+sigma**2)


	arg2 = B**2*(1+sigma**2)/sigma**4


	return 1/arg1*math.e**(-arg2/2)



B_lis = [0, 0.5, 1.0, 2.0]



lin = ['solid', 'dashed', 'dashdot', 'dotted']



for i in range(len(B_lis)):


	B = B_lis[i]


	est_mean_W_lis = [estationary_mean_W(sigma, B) for sigma in sigma_lis]


	lab = sep.join(('$B/\omega=', str(B), '$'))


	plt.plot(sigma_lis, est_mean_W_lis, color='black', label=lab, linestyle=lin[i], marker='o', linewidth=3, markersize=0)



plt.plot([sigma_inf, sigma_sup], [0, 0], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.plot([sigma_inf, sigma_sup], [1, 1], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.xlim([sigma_inf, sigma_sup])
#plt.ylim()
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\sigma_*/\omega$', size=LAB, labelpad=PAD)
plt.ylabel('$\overline{W}_*/W_{max}$', rotation=90, size=LAB, labelpad=PAD)
#plt.ylabel('$\\frac{\overline{W}_*}{W_{max}}$', rotation=0, size=25, labelpad=PAD)
plt.title('$B$', size=TIT)
plt.legend(loc='best', fontsize=LAB)
#plt.legend(shadow=False, frameon=False, prop={'size': LAB})
plt.xticks(size=TIC)
plt.yticks(size=TIC)



plt.subplot(1, 3, 3)



def stationary_sigma(mu):


	arg = mu*mu + math.sqrt(mu**4+4*mu**2)


	return math.sqrt(arg/2)



def critical_B(sigma, Wmax):


	arg = Wmax/math.sqrt(sigma*sigma+1)


	return sigma*sigma/math.sqrt(sigma*sigma+1)*math.sqrt(2*math.log(arg))



mu_inf = pow(10, -3)
mu_inc = pow(10, -4)



Wmax_lis = [1.1, 1.3, 1.5, 2.0]



for i in range(len(Wmax_lis)):


	Wmax = Wmax_lis[i]


	if stationary_sigma(mu_inf) < math.sqrt(Wmax*Wmax-1):


		mu_lis = [mu_inf]


		while stationary_sigma(mu_lis[-1]+mu_inc) < math.sqrt(Wmax*Wmax-1):


			mu_lis.append(mu_lis[-1]+mu_inc)


		Bc_lis = [critical_B(stationary_sigma(mu), Wmax) for mu in mu_lis]
		
		
		Bc_lis[-1] = 0


		lab = sep.join(('$W_{max}=', str(Wmax), '$'))


		plt.plot(mu_lis, Bc_lis, color='black', label=lab, linestyle=lin[i], marker='o', linewidth=3, markersize=0)



gap = 0.2



plt.plot([mu_inf, mu_lis[-1]+gap], [0, 0], color='blue', linestyle='solid', marker='o', linewidth=1, markersize=0)



plt.xlim([mu_inf, mu_lis[-1]+gap])
#plt.ylim()
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\delta/\omega$', size=LAB, labelpad=PAD)
plt.ylabel('$B_c\,/\omega$', rotation=90, size=LAB, labelpad=PAD)
#plt.ylabel('$\\frac{B_c}{\omega}$', rotation=0, size=25, labelpad=PAD)
plt.title('$C$', size=TIT)
plt.legend(loc='best', fontsize=LAB)
#plt.legend(shadow=False, frameon=False, prop={'size': LAB})
plt.xticks(size=TIC)
plt.yticks(size=TIC)



fig.set_size_inches(22, 7, forward=True)
fig.savefig('figure-1.pdf', bbox_inches='tight')




#################### Bottom ####################
