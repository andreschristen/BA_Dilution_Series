# -*- coding: utf-8 -*-
"""
Created on Mon 2018.02.28:00:00:00

@author: jac-nosm

Abstrac:
Dilution experiments are performed routinely in microbial laboratories.
However, their results are not properly analyzed, with only very basic hueristic
formulas.  We present a proper statistical analysis of these experiments and
include also a further analysis of the bacterial evolution under temperature stress. 
"""

# Packages

from numpy import array, zeros, exp

from DilExp import DilExp, MultiDilExp
from matplotlib import pyplot as plt


def TransformTNTC( count, c=300):
    if isinstance( count, int):
        return count
    else:
        return c

def AnaBF( spreadsheet, lab, T=100000, control=False, data_all=True):
    """Analyse BF for the repetitions in all data
       experiment in lab = 'Lab5', 'Lab6', 'Lab8' or
       if control 'Lab5Cntrl', 'Lab6Cntrl', 'Lab8Cntrl'
       not(data_all), that is, only include the first countable dilution.
    """
    J=7 ## Dilutions < J
    c=300 ##maximum count for the drop or plated volumen, TNTC threshold
    d = DilExp( J=J, alpha0=4, alpha=10, alphap=100, c=c, Sc=4.13, Sc_unit=r"cm^2", q=0.05, pr_a=10)
    """     alpha0, dilution factor for tube 1 from tube 0: 1ml from 4*10ml
            alpha, dilution factor for each tube =10, 1ml from 10ml
            alphap, dilution factor for the drop or plated volumen =100, 0.1ml from 10ml
            q, probibility of misscounting (=0.05)
            pr_a, prior mean of the count, given the dilution s to be counted (=20)
            Sc, coupon surface area (=4.13)
            Sc_unit, coupon surface area units (=r"cm^2").
    """
    if control:
        sheet = lab + 'Cntrl'
        #print("Analizing control experiment of %s" % (lab,))
    else:
        sheet = lab
        #print("Analizing experiment of %s" % (lab,))
    rt = []
    for k in range(3): #Repetitions 1, 2, 3
            s = zeros(2)
            y = zeros(2)
            plate1 = array(spreadsheet[sheet].loc[ :, 'Plate 1'])
            for j, count in enumerate(plate1[J*k:(J*(k+1))]):
                if TransformTNTC(count) < c:
                    s[0] = spreadsheet[sheet].loc[ j+J*k, 'Dilution']
                    y[0] = count
                    break
            plate2 = array(spreadsheet[sheet].loc[ :, 'Plate 2'])
            for j, count in enumerate(plate2[J*k:(J*(k+1))]):
                if TransformTNTC(count) < c:
                    s[1] = spreadsheet[sheet].loc[ j+J*k, 'Dilution']
                    y[1] = count
                    break
            d.Data( s=s, y=y, calc_post=True)
            rt += [d.BF]
            print("%10s, %2d, %16g, %16g, %16g, %16f" %\
                  ( sheet, k, d.k_betabinom, d.k_binom, exp(d.k_binom_c - d.k_betabinom_c), d.BF))
    return rt


def InterLabGenMultiDilExp( spreadsheet, lab, T=100000, control=False, data_all=True):
    """Create a MultiDilExp object for the repetitions of
       experiment in lab = 'Lab5', 'Lab6', 'Lab8' or
       if control 'Lab5Cntrl', 'Lab6Cntrl', 'Lab8Cntrl'
       if data_all include all the counts including TNTC
       if not(data_all) only include first countable dilution.
    """
    J=7 ## Dilutions < J
    c=300 ##maximum count for the drop or plated volumen, TNTC threshold
    md = MultiDilExp( K=3, J=J, alpha0=4, alpha=10, alphap=100, c=c, Sc=4.13, Sc_unit=r"cm^2", q=0.05, pr_a=10, b=500, M=10)
    """     alpha0, dilution factor for tube 1 from tube 0: 1ml from 4*10ml
            alpha, dilution factor for each tube =10, 1ml from 10ml
            alphap, dilution factor for the drop or plated volumen =100, 0.1ml from 10ml
            q, probibility of misscounting (=0.05)
            pr_a, prior mean of the count, given the dilution s to be counted (=20)
            Sc, coupon surface area (=4.13)
            Sc_unit, coupon surface area units (=r"cm^2").
    """
    if control:
        sheet = lab + 'Cntrl'
        print("Analizing control experiment of %s" % (lab,))
    else:
        sheet = lab
        print("Analizing experiment of %s" % (lab,))
    for k in range(3): #Repetitions 1, 2, 3
        if data_all:
            s = array([[d,d] for d in range(J)]).flatten() # Dilutions
            y = zeros(2*J)
            for d in range(J):
                y[2*d] = TransformTNTC(spreadsheet[sheet].loc[ d+J*k, 'Plate 1'])
                y[2*d + 1] = TransformTNTC(spreadsheet[sheet].loc[ d+J*k, 'Plate 2'])
            print("Repetition %d" % (k,))
            print(s)
            print(y)
            print("\n")
            md.Data( k=k, s=s, y=y)
        else:
            s = zeros(2)
            y = zeros(2)
            plate1 = array(spreadsheet[sheet].loc[ :, 'Plate 1'])
            for j, count in enumerate(plate1[J*k:(J*(k+1))]):
                if TransformTNTC(count) < c:
                    s[0] = spreadsheet[sheet].loc[ j+J*k, 'Dilution']
                    y[0] = count
                    break
            plate2 = array(spreadsheet[sheet].loc[ :, 'Plate 2'])
            for j, count in enumerate(plate2[J*k:(J*(k+1))]):
                if TransformTNTC(count) < c:
                    s[1] = spreadsheet[sheet].loc[ j+J*k, 'Dilution']
                    y[1] = count
                    break
            print("Repetition %d" % (k,))
            print(s)
            print(y)
            print("\n")
        md.Data( k=k, s=s, y=y)
    md.data_other = sheet
    md.RunTwalk(T=T)
    md.twalk.Ana()
    return md

def InterLabMakeList( spreadsheet, T=100000, control=False, data_all=False):
    mdlist = []
    for lab in ['Lab5', 'Lab6', 'Lab8']:
        mdlist += [InterLabGenMultiDilExp( spreadsheet, lab, T=T, control=control, data_all=data_all)]
    return mdlist

