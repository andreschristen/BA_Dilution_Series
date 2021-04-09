#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 14:55:08 2018

@author: jac

Inter lab analhysis for Clorox Low 1 experiment
"""
from pickle import load

from numpy import diff, linspace, var, array, mean, sqrt, abs, append

from scipy.stats import gaussian_kde, uniform
from scipy.stats.mstats import mquantiles


from pandas import read_excel

from pylab import savefig, rc, tight_layout, close, figure, plot, xlabel, ylabel, boxplot, ylim, hlines
rc('font', size=18)

from InterLabDataAnalysis import InterLabMakeList, AnaBF

from DilExp import InterLabDilExp

def AnaLRs( Il, Il_control, mdlist, mdlist_control, col, alpha=0.001, N=200, Sc_unit=r"cm^2"):
    LRs = [mdlist_control[i].TwalkE() - mdlist[i].TwalkE() for i in range(len(mdlist))]
    cLR = Il_control.TwalkcE() - Il.TwalkcE()
    
    Dict = {} #Dictonary to hold all posterior expectations and variances
    ### cE is \mathcal{E}
    Dict['Exp_cE']  = array([ mean(Il_control.TwalkcE()), mean(Il.TwalkcE()),  mean(cLR)])
    Dict['Vars_cE'] = array([ var(Il_control.TwalkcE()), var(Il.TwalkcE()),  var(cLR)])

    Dict['Exp_E']   = array([mean(mdlist[i].TwalkE()) for i in range(len(mdlist))])
    Dict['Vars_E']  = array([var(mdlist[i].TwalkE()) for i in range(len(mdlist))])
    
    Dict['Exp_E_control']  = array([mean(mdlist_control[i].TwalkE()) for i in range(len(mdlist))])
    Dict['Vars_E_control'] = array([var(mdlist_control[i].TwalkE()) for i in range(len(mdlist))])

    Dict['Exp_LR']   = array([mean(LRs[i]) for i in range(len(mdlist))])
    Dict['Vars_LR']  = array([var(LRs[i]) for i in range(len(mdlist))])

    e0=15
    e1=0
    for i,LR in enumerate(LRs):    
        q = mquantiles( LR, [alpha,1-alpha])
        e0 = min(q[0],e0)
        e1 = max(q[1],e1)
    e = linspace( e0, e1, N)
    
    for i,LR in enumerate(LRs):
        print("Calculating kde for free LR_%d ..." % (i,))
        kde_instance = gaussian_kde( LR, bw_method=lambda kde: kde.n**(-1.0/((kde.d + 4)*1.6)))
        kde = kde_instance.evaluate(e)
        C = sum(diff(e)*kde[:-1]) ### Keep the integral == 1
        kde /= C
        plot( e, kde, '-', color=col[i], linewidth=1)
        print("Interlab, lab %d, P[ LR > 3 ] = %6.4f" % (i, sum(LR > 3)/len(LR)))
        print("Interlab, E[|\mathcal{E} - E_%d|] = %6.4f\n" % (i, mean(abs(cLR[:-1:3]-LR[:-1]))) )

    print("Calculating kde of \mathcal{E} ...")
    kde_instance = gaussian_kde( cLR, bw_method=lambda kde: kde.n**(-1.0/((kde.d + 4)*1.6)))
    kde = kde_instance.evaluate(e)
    C = sum(diff(e)*kde[:-1]) ### Keep the integral == 1
    kde /= C
    plot( e, kde, 'k-', linewidth=2)
    xlabel(r"$log_{10}\left(\frac{CFU_0 + 1}{CFU + 1}\right)$")
    ylabel("Density")
    print("Interlab, global, P[ LR > 3 ] = %6.4f" % (sum(cLR > 3)/len(cLR),))
    
    return Dict



if __name__ == '__main__':
    T=200000

    ### BF plot
    InterLab_CloroxLow1 = read_excel( './Data/InterLab_CloroxLow1.xls',\
            ['Lab5', 'Lab6', 'Lab8', 'Lab5Cntrl', 'Lab6Cntrl', 'Lab8Cntrl'])
    print("%10s, %2s, %16s, %16s, %16s, %16s" %\
              ( "experiment", "k" , "bbinom" , "binom" , "corr" , "BF"))
    rt_InterLab = []
    for lab in ['Lab5', 'Lab6', 'Lab8']:
        rt_InterLab += AnaBF( InterLab_CloroxLow1, lab, control=False)
        rt_InterLab += AnaBF( InterLab_CloroxLow1, lab, control=True)
    close(1)
    figure(1)
    # rt_CBE is created in CBEExamples
    rt_CBE = load(open("CBE_rt_CBE.pkl", "rb"))
    boxplot( [rt_CBE, rt_InterLab], positions=[1,2], labels=['CBE', 'Interlab'])
    plot( [1]*len(rt_CBE), rt_CBE, 'k.')
    plot( [2]*len(rt_InterLab), rt_InterLab, 'k.')
    hlines( 3, 0, 3, color='red')
    hlines( 1, 0, 3, color='black', linestyles='dashed')
    ylim((-0.5,4))
    ylabel("Bayes Factors")
    savefig("Images/BFs.jpg")
    all_bf = append( array(rt_CBE), array(rt_InterLab))
    print("BF outliers > 4: 70CTemp_10min, k=2, BF=15, CloroxLow1, Lab5, k=0, BF=15, lab6, k=0, BF=5.5")
    all_bf = append( array(rt_CBE), array(rt_InterLab))
    print("BFs < 1: %d, of %d." % (sum(all_bf < 1.0), len(all_bf)))


    mdlist = InterLabMakeList( InterLab_CloroxLow1, T=T, control=False, data_all=False)
    mdlist_control = InterLabMakeList( InterLab_CloroxLow1, T=T, control=True, data_all=False)

    print("\n\nAnalizing hierarchycal model of experiment Clorox low1 for all labs.")
    Il = InterLabDilExp(mdlist)
    Il.RunTwalk(T=3*T)
    close(1)
    figure(1)
    kde_list = Il.AnaEs(col=['blue','blue','blue'])
    tight_layout()
    savefig("Images/InterLab%s_PostEs.jpg" % ('CloroxLow1',))

    print("\n\nAnalizing hierarchycal model of controls for all labs.")
    Il_control = InterLabDilExp(mdlist_control)
    Il_control.RunTwalk(T=3*T)
    close(2)
    figure(2)
    kde_list_control = Il_control.AnaEs(col=['blue','blue','blue'])
    xlabel(r"$log_{10}\left( CFU_0 + 1\right)$")
    tight_layout()
    savefig("Images/InterLab%s_PostEs_control.jpg" % ('CloroxLow1',))

    close(3)
    figure(3)
    ExpVars = AnaLRs( Il, Il_control, mdlist, mdlist_control, col=['blue','blue','blue'])
    tight_layout()
    savefig("Images/InterLab%s_PostLRs.jpg" % ('CloroxLow1',))
    
    print("\nPost mean for \mathcal{E}: cntrl, Exp, LR:", ExpVars['Exp_cE'])
    print("Variance  for \mathcal{E}: cntrl, Exp, LR:", ExpVars['Vars_cE'])

    print("\nPost mean for control E:", ExpVars['Exp_E_control'], mean(ExpVars['Exp_E_control']))
    print("Variance  for control E:", ExpVars['Vars_E_control'], mean(ExpVars['Vars_E_control']))

    print("\nPost mean for E:", ExpVars['Exp_E'], mean(ExpVars['Exp_E']))
    print("Variance  for E:", ExpVars['Vars_E'], mean(ExpVars['Vars_E']))

    print("\nPost mean for LR:", ExpVars['Exp_LR'], mean(ExpVars['Exp_LR']))
    print("Variance  for LR:", ExpVars['Vars_LR'], mean(ExpVars['Vars_LR']))


