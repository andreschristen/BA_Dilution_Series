#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:50:46 2019

An example of how to run the program entering the data manually.
All settings may be changed, changing the default values.

Data are taken from the CBE example, 70 C, 10 min.
See Excel file CBE_BiofilmHotWaterStudies.xls

Change data or default parameter values as required.

JA Christen, jac@cimat.mx

"""

from numpy import array
from pylab import figure, tight_layout

from DilExp import DilExp, MultiDilExp


### Single dilution experiment, default values, see DilExp?
dilexp_single = DilExp()
dilexp_single.Data( s=array([0]), y=array([ 1]))

### Repetition (tube) 2: All 10 drops were counted at dilution 0, y are the counts
dilexp_single.Data( s=array([0]*10), y=array([14, 11, 14, 16, 17, 22, 32, 18, 11, 21]), calc_post=True)

figure(1)
#dilexp_single.PlotPost(log10_p1=False) ###No MCMC
#dilexp_single.PlotPost(betabinom=True, log10_p1=False)
#print("bbinom %g, binom %g, BF: %f" % ( dilexp_single.k_betabinom, dilexp_single.k_binom, dilexp_single.k_betabinom/dilexp_single.k_binom))

tight_layout()


###Multidilution experiment: all 3 tubes of CBE example, 70 C, 10 min.
### Check default values with MultiDilExp?
md = MultiDilExp()

md.Data( k=0, s=array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0]),\
              y=array([2, 4, 3, 6, 1, 2, 5, 7, 3, 1]))
md.Data( k=1, s=array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),\
              y=array([14, 11, 14, 16, 17, 22, 32, 18, 11, 21])) #Same as dilexp_single
md.Data( k=2, s=array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0]),\
              y=array([0, 0, 0, 0, 0, 7, 5, 0, 2, 9]))

md.RunTwalk( T=100000 )
figure(2)
# To analyze the information of the t-walk MCMC
md.twalk.Ana()
tight_layout()

figure(3)
md.PlotResults( fignum=3 )
tight_layout()










