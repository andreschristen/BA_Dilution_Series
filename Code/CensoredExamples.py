#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 14:21:24 2018

@author: jac
"""


from numpy import array, ones, exp, arange

from scipy.stats import binom


from pylab import plot, figure, savefig, subplots, close, rc, tight_layout
rc('font', size=18)

from DilExp import DilExp

close(3)
close(4)



############################################################
### In the standard drop plate method, post are the same even with q=0.0
############################################################
close(1)
close(2)

m=10 # # of drops

fig, ax = subplots(num=2) ###Figure 2
de = DilExp(q=0.05) #Default values, drop plate system
de.Data( s=array([0]*m), y=array([0]*m)) ### All obs == to zero
figure(1)
de.PlotPost(CFU_X_Area=False, log10_p1=False)
ax.plot( de.x, de.fx, 'b-', label=r"$q = 0.05$", linewidth=2)

de = DilExp(q=0.0) #Default values, drop plate system
de.Data( s=array([0]*m), y=array([0]*m))
de.CalculatePost()
figure(1)
de.PlotPost(CFU_X_Area=False, log10_p1=False)
ax.plot( de.x, de.fx, 'g-', label=r"$q = 0.00$", linewidth=2)

ax.set_xlabel("CFU")
ax.set_ylabel(r"Probability")
ax.legend(loc='upper right', shadow=True)#, fontsize='large')

figure(2)
tight_layout()
savefig("Images/Censored1.jpg")



####################################################################
### Differences in using censored likelihood or ignoring censorship
###################################################################
figure(3)
de = DilExp() #Default values
m=5  ##### ALL observatiions are sturated
de.Data( s=array([6]*m), y=array([de.c]*m))
de.PlotPost( CFU_X_Area=False, log10_p1=False, plot_data=False, color='green')

censored_pmf = ones(de.k_range.shape)
for i, k in enumerate(de.k_range):
    #for y in de.y:
    #    if y < de.c: ###Not censored
    #        censored_pmf[i] *= binom.pmf( k=y, n=k, p=de.dil[0])
    #    else: ### Censored
    #        censored_pmf[i] *= sum(binom.pmf( k=arange(de.c,2*de.prior_range[-1]), n=k, p=de.dil[0]))
    #print( log(censored_pmf[i]), de.LogLikelihoodCensored(k))
    censored_pmf[i] = exp(de.LogLikelihoodCensored(k))
censored_pmf /= sum(censored_pmf)
plot( de.k_range, censored_pmf, 'b-')
plot( [de.k_range[-1]]*2, [0,censored_pmf[-1]], 'k--')
plot( de.y*de.mul, [0.0]*de.m, 'k*')
tight_layout()
savefig("Images/Censored2.jpg")




figure(4)
de = DilExp() #Default values
#####ONLY THREE non TNTC
de.Data( s=array([6]*m), y=array([de.c]*(m-3)+[25]*3))
de.PlotPost( CFU_X_Area=False, log10_p1=False, plot_data=False, color='green')

censored_pmf = ones(de.k_range.shape)
for i, k in enumerate(de.k_range):
    ##Censored likelihood prod_i^m P(Y_i > c | k)  
    for y in de.y:
        if y < de.c: ###Not censored
            censored_pmf[i] *= binom.pmf( k=y, n=k, p=de.dil[0])
        else: ### Censored
            censored_pmf[i] *= sum(binom.pmf( k=arange(de.c,2*de.prior_range[-1]), n=k, p=de.dil[0]))
censored_pmf /= sum(censored_pmf)
plot( de.k_range, censored_pmf, 'b-')
plot( [de.k_range[-1]]*2, [0,censored_pmf[-1]], 'k--')
plot( de.y*de.mul, [0.0]*de.m, 'k*')
tight_layout()
savefig("Images/Censored3.jpg")


### There is three orders of magnitud difference in time
#timeit de.LogLikelihoodCensored(de.k_range[250]) ## 424 µs ± 27.7 µs
#timeit de.LogLikelihood(de.k_range[250])  ### 8.94 µs ± 432

