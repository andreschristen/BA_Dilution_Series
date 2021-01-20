#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 14:37:49 2018

@author: jac
"""

from numpy import log10, array, arange, cumsum, histogram, where, linspace, diff

from scipy.stats import gaussian_kde, gamma, uniform


from pylab import figure, savefig, subplots, close, rc, tight_layout, plot
rc('font', size=18)

from DilExp import DilExp, MultiDilExp


def PlotMargPriorN( M, b, N, T=1000000):
    """Plot the marginal prior dist for N0, in hierarchycal model with 1 rep."""
    e = M*uniform.rvs(size=T) #gamma.rvs( a=1, scale=2*M, size=T) #M*uniform.rvs(size=T)
    a = gamma.rvs( a=1, scale=b, size=T)
    S = gamma.rvs( a=a, scale=e/a, size=T)#gamma.rvs( a=a*e**2, scale=1/(a*e), size=T)#gamma.rvs( a=a, scale=e/a, size=T)#
    N0 = 10**M * uniform.rvs(size=T)
    fig, (ax1, ax2, ax3) = subplots(nrows=3, ncols=1, num=2)
    ax1.hist( S, density=True, bins=linspace( 0, M, num=160))
    ax2.hist( 10**S - 1, density=True, bins=N)
    ax3.hist( log10(N0+1), density=True, bins=linspace( 0, M, num=160))
    ax1.set_xlabel(r"$log_{10}(CFU + 1)$")
    ax2.set_xlabel("CFU")
    ax3.set_xlabel(r"$log_{10}(CFU + 1)$")
    ax1.set_ylabel(r"Probability")
    tight_layout()
    return S



if __name__=='__main__':

    #### Level of detection if all dilution is platted. ####
    de = DilExp( alphap=1.0, q=0.05) #all dilution is platted
    m=1
    de.Data( s=array([0]*m), y=array([0]*m))
    de.CalculatePost()
    print("alpha_p = 1.0, all dilution platted, P[N^k_0 = ", de.k_range[:3], " | Y=0] = ", de.post_array[:3])  

    #### Level of detection hierachycal model with 1, 3 and 6 repetiions. ####

    close(1)
    step = 10
    N = arange( 0, 500, step, dtype=int) ###Plotting area, every 10 CFUs
    m=10  ###m drops
    b=500

    ##### Calculate also E for various repetitions ####
    fig, ax = subplots(num=1) ###Figure 1
    print("b = %d" % (b,))
    Ks = [1, 3, 12]
    LOD = [0]*len(Ks)
    marker = [ '*', '.', 'o', '>', '<']
    for i,K in enumerate(Ks):
        print("\n### Analizing LOD for K=%d #####" % (K,))
        md = MultiDilExp(K=K, b=b, Sc=1.0) ###Default values for drop plate
        for k in range(K):
            md.Data( k, s=array([0]*m), y=array([0]*m))
        md.RunTwalk(T=300000) 
        rt = histogram( 10**md.twalk.Output[:,2] - 1, bins=N, density=True  )
        step_array = rt[0]*step
        ax.plot( N[:-1], step_array, marker=marker[i],\
            ls='None', markersize=4, label=r"$K = %d$" % (K,))
        LOD[i] = N[where(cumsum(step_array) >= 0.95)[0][0]]

    for i,K in enumerate(Ks):
        print("95%% level of detection with K=%d rep.: %d" % ( K, LOD[i]))

    ax.set_xlabel("CFU")
    ax.set_ylabel(r"Probability")
    ax.legend(loc='upper right', shadow=True)#, fontsize='large')

    figure(1)
    tight_layout()
    savefig("Images/LevelOfDetection.jpg")

    close(2)
    figure(2)
    step=10
    N = arange( 0, 500, step, dtype=int) ###Plotting area, every 10 CFUs
    b=500
    M=10
    S = PlotMargPriorN( M=M, b=b, N=N)
    kde = gaussian_kde( S, bw_method=lambda kde: kde.n**(-1.0/((kde.d + 4)*1.6)))
    e = linspace( 0, M, 160)
    kde = kde.evaluate(e)
    C = sum(diff(e)*kde[:-1]) ### Keep the integral == 1
    close(3)
    figure(3)
    plot( e, kde/C, 'g-')

