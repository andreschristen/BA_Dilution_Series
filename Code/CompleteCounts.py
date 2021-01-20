#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 14:21:24 2018

@author: jac
"""

from time import time, localtime, strftime

from numpy import array, ones, exp, arange, log, cumsum

from scipy.stats import binom

from pylab import plot, figure, savefig, xlabel, ylabel, legend, subplots, close, rc, tight_layout
rc('font', size=18)

from DilExp import DilExp

close(5)


##################### Analyzyng ALL counts across dilutions vs
### only the first obs count
### We do not need to worry for TNTC ... the likelihood is flat
### for higher values


def PerformExpPost( N0, T, m=5, J=7, alpha0=1, alpha=10, alphap=1000, c=30, q=0.05, pr_a=10):
    """Expected posterior calculation of using all counts vs first count.
       We do not need to worry for TNTC, the likelihood is flat and does not provide info.
       Then se assume the *first* dilution is countable, for this, N0 needs to be
       0 to mul * c * alpha**-1 * alphap**-1, where mul = 1 to 2 perhaps.
       Otherwise the expected value for the first counts is (well) above c, the TNTC level.

       N0: number of CFUs in dilution 0.
       T: number iterations, repeated samples.
       
       
    """
    print("########## PerformExpPost ############")
    print("\nN_0=%d, T=%d, J=%d, alpha0=%d,alpha=%d, alphap=%d, c=%d, q=%f, pr_a=%f" %\
          ( N0, T, J, alpha0, alpha, alphap, c, q, pr_a))

    de = DilExp( J=J, alpha0=alpha0, alpha=alpha, alphap=alphap, c=c, q=q, pr_a=pr_a) 
    de2 = DilExp( J=J, alpha0=alpha0, alpha=alpha, alphap=alphap, c=c, q=q, pr_a=pr_a) 

    de.SimData(N0=N0,s=array([0]*m)) ##To make de.dil

    print("Expected values al dilution 0:", de.dil*N0)
    de.Data( s=array([0]*m), y=de.dil*N0, calc_post=True) ##Expected values
    k_range = de.k_range #Fix the range
    s = []
    for i in range(J):
        s += [i]*m
    s = array(s) #Full set of dilutions and drops
    for t in range(T):
        ### Simulate full data
        y2 = de2.SimData( N0=N0, s=s)
        while all(y2[:m] >= de.c): ##All aboce TNTC
            print("Rejected, all above %d, TNTC." % (de.c,))
            y2 = de2.SimData( N0=N0, s=s)
        print( t, y2[:m], y2[m:2*m], y2[2*m:3*m], y2[2*m:3*m], "..." )
    
        ##For de take the First counted dil only
        de.Data( s=array([0]*m), y=y2[:m]) 
        tmp = array([de.LogPost(k) for k in k_range])
        tmp += 200 - max(tmp)  ### Set the maximum to 200 for numerical stability
        tmp = exp(tmp)
        if t==0:
            de.post_array  = tmp/sum(tmp)
        else:
            de.post_array += tmp/sum(tmp)

        ##For de2 take all counted dil
        de2.Data( s=s, y=y2) ##ALL data
        tmp = array([de2.LogPost(k) for k in k_range])
        tmp += 200 - max(tmp)  ### Set the maximum to 200 for numerical stability
        tmp = exp(tmp)
        
        ### Average posteriors
        if t==0:
            de2.post_array  = tmp/sum(tmp)
        else:
            de2.post_array += tmp/sum(tmp)

    #### Plot expected posteriors
    cdf = cumsum(de.post_array/T)
    q1 = 0.0001  ### Range within these quatiles
    q2 = 1-q1
    j1 = 0
    j2 = len(k_range)-1
    for j,p in enumerate(cdf):
        if (p <= q1):
            j1 = j
        if (p > q2):
            j2 = j
            break

    print("First counted dilution, green, and all dilutions counted, blue.\n")
    plot( k_range[j1:(j2+1)],  de.post_array[j1:(j2+1)]/T, 'g.', markersize=1)
    plot( k_range[j1:(j2+1)], de2.post_array[j1:(j2+1)]/T, 'b.', markersize=1)
    xlabel(r"$N_0$")
    ylabel("Probability")
    tight_layout()
    
    return k_range, de.post_array/T, de2.post_array/T


if __name__=='__main__':

    T=120
    ### Time esrtimate, 11 secs for drop plate and 23 secs for plated, per iteration.
    est_processing_secs = (11*T*5 + 23*T*5)
    print("Estimated procesing time, %d hours and %d minutes." %\
          (est_processing_secs // 3600, (est_processing_secs % 3600) // 60))

    inittime1 = time()
    print("############ Start Drop plate: %s" % (strftime("(%a), %Y.%m.%d:%H:%M:%S", localtime(inittime1)),))
     
    for N0 in [ 500, 5000, 10000, 20000, 30000]:
        close(1)
        figure(1)
        k_range, post_first1, post_all1 = PerformExpPost( N0=N0, T=T) ##CBE drop plate design
        savefig("Images/AllCounts_DropPlate_N0=%d.jpg" % (N0,))

    finishtime1 = time()
    print("############# Start Drop plate: %s" % (strftime("(%a), %Y.%m.%d:%H:%M:%S", localtime(inittime1)),))
    print("############ Finish Drop plate: %s" % (strftime("(%a), %Y.%m.%d:%H:%M:%S", localtime(finishtime1)),))
    inittime2 = time()
    print("\n################# Start plated: %s" % (strftime("(%a), %Y.%m.%d:%H:%M:%S", localtime(inittime2)),))
      
    for N0 in [ 50, 500, 1000, 2000, 3000]:
        close(2)
        figure(2)
        k_range, post_first1, post_all1 = PerformExpPost( N0=N0, T=T,\
            m=2, J=7, alpha0=4, alpha=10, alphap=100, c=300) ##Interlab plated design
        savefig("Images/AllCounts_Interlab_N0=%d.jpg" % (N0,))

    finishtime2 = time()
    print("############# Start Drop plate: %s" % (strftime("(%a), %Y.%m.%d:%H:%M:%S", localtime(inittime1)),))
    print("############ Finish Drop plate: %s" % (strftime("(%a), %Y.%m.%d:%H:%M:%S", localtime(finishtime1)),))
    print("########### Processing: %d minutes and %d seconds." %\
      ((finishtime1-inittime1) // 60, (finishtime1-inittime1) % 60))
    print("\n################# Start plated: %s" % (strftime("(%a), %Y.%m.%d:%H:%M:%S", localtime(inittime2)),))
    print("################ Finish plated: %s" % (strftime("(%a), %Y.%m.%d:%H:%M:%S", localtime(finishtime2)),))
    print("########### Processing: %d minutes and %d seconds." %\
      ((finishtime1-inittime1) // 60, (finishtime2-inittime2) % 60))

