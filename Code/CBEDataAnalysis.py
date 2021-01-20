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
from time import sleep

from scipy.special import gammaln
from scipy.stats import expon, gaussian_kde
from scipy.stats.mstats import mquantiles

from numpy.random import choice
from numpy import append, array, where, arange, linspace, matrix
from numpy import log, log10, exp, var, sqrt, diff, cumsum, mean, inf

from pylab import plot, vlines, xlabel, ylabel, xlim, ylim, subplot, title
from pylab import figure, savefig, boxplot, close, hlines, legend, axis, tight_layout

from pytwalk import pytwalk

from DilExp import DilExp, MultiDilExp, GetKDE




#### CBE Data: define the range of the graphs to each experiment 
PlotPostXlim = {'65CTemp_15min':[4.0,6.5], '65CTemp_30min':[3.5,6.0], '65CTemp_45min':[3.0,6.0],\
                '65CTemp_60min':[2.0,6.0], '65CTemp_90min':[0.0,6.0],\
                '70CTemp_10min':[-1.0,6.0], '70CTemp_20min':[2.0,6.0], '70CTemp_30min':[-4.0,6.0],\
                '70CTemp_40min':[0.0,4.0], '70CTemp_60min':[-2.0,3.0],\
                '75CTemp_10min':[-2.2,5.0],'75CTemp_20min':[-2.0,5.0],\
                '80CTemp_1min': [ 1.5,4.5],'80CTemp_2min': [-0.1,4.6],\
                'RoomTemp_1min':[ 7.5,10.0],'RoomTemp_10min':[7.0,9.5],'RoomTemp_15min':[7.0,9.5],\
                '21CTemp_10min':[ 4.5,9.5],'21CTemp_11min':[4.5,9.5],'21CTemp_12min':[4.5,9.5],'21CTemp_13min':[4.5,9.5]}


def AnaBF(experiment, time, CBEData, D=[5,5]):
    """Calculate the Bayes factor of data from experiment, time for CBEData.
       Return the results."""
    d = DilExp() #Default values
    rt = []
    #D=[5,5] ## Two plates of five drops
    ### Enter data
    sel = where(array(CBEData[experiment].loc[ :, 'Time'] == time))[0]
    for k,i in enumerate(sel[::(len(D))]): 
        # To creat a list with all the information analized of an experiment (temperature_time)
        s, y = [], []
        for r in range(len(D)):
            s += [CBEData[experiment].loc[ i+r,'Dilution']]*D[r]
            y = append(y , CBEData[experiment].loc[ i+r,['Drop%d'%j for j in range(1,D[r]+1)]]) 
        d.Data( s=s, y=y, calc_post=True)
        rt += [d.BF]
        Z_beta_binom = exp(-d.k_betabinom_c)*d.k_betabinom
        Z_binom = exp(-d.k_binom_c)*d.k_binom
        print("%10s_%02dmin, %2d, %16g, %16g, %16f, %16f" %\
          (experiment, time, k, Z_beta_binom, Z_binom, Z_binom/(Z_beta_binom+Z_binom), Z_beta_binom/Z_binom))
    return rt

       

def GenMultiDilExp(experiment, time, CBEData, D=[5,5], betabinom=False):
    """Create a MultiDilExp object for CBE data with 'experiment' and 'time'."""
    md = MultiDilExp(betabinom=betabinom) #Default values
    #D=[5,5] ## Two plates of five drops
    ### Enter drop data
    sel = where(array(CBEData[experiment].loc[ :, 'Time'] == time))[0]
    for k,i in enumerate(sel[::(len(D))]): 
        # To creat a list with all the information analized of an experiment (temperature_time)
        s, y = [], []
        for r in range(len(D)):
            s += [CBEData[experiment].loc[ i+r,'Dilution']]*D[r]
            y = append(y , CBEData[experiment].loc[ i+r,['Drop%d'%j for j in range(1,D[r]+1)]]) 
        md.Data( k=k, s=s, y=y)
        ### Should be the same for i and i+1
        md.experiment = experiment
        md.data_other = CBEData[experiment].loc[ i, ['Temperature', 'Time', 'Tube', 'Plate']]
        md.PlotPostXlim = PlotPostXlim["%s_%dmin" % (experiment, time)]
    return md


def AnaTimeData( experiment, time, CBEData, T=100000, fig=[0,1], betabinom=False):
    """Function to analyze the information obteined of t-walk"""
    md = GenMultiDilExp(experiment=experiment, time=time, CBEData=CBEData, betabinom=betabinom)
    print(("\n\nExperiment: %s_%dmin" % (experiment, time)))
    # To run the t-walk
    md.RunTwalk(T=T)
    figure(fig[0])
    # To analyze the information of the t-walk
    md.twalk.Ana()
    savefig("Images/%s_%dmin_TS.png" % (experiment, time))

    figure(fig[1])
    print(" ")
    # Results graphics
    if betabinom:
        md.PlotResults(fignum=fig[1], color="green")
    else:
        md.PlotResults(fignum=fig[1], color="blue")
    xlim(md.PlotPostXlim)
    tight_layout()
    #savefig("Images/%s_%dmin_Results.jpg" % (experiment, time))
    return md


def PlotIndPost( md, fig=[0,1], K=[0,1], color="skyblue"):
    """Plot individual posteriors, without hierachycal model, in
       multi dilution experiment object md.  
       USE THIS, K=[0,1], TO CREATE FIG 4 OF THE PAPER, with '80CTemp', '2min'."""
    figure(fig[0]) ### Add to current joint posterior 
    for k in K:
        md.d[k].PlotPost(linestyle='dashdot', color = color)
    figure(fig[1])
    for k in K:
        plot( md.d[k].k_range[0:1000], md.d[k].post_array[0:1000], 'b-', linestyle='dashdot', color = "limegreen")
    xlabel(r"$N_0$")
    ylabel(r"Probability")
    


      
def PlotTemp(md_list):
    """ First function to analyze the results of the activation threshold"""
    close(3)
    figure(3)
    xlim((-1,62))
    hlines( 0, -1, 62)
    ylim((-0.3,10.2))
    vlines( 0, -0.2,10.2)
    for md in md_list: 
            md.PlotvPar( par=1, position=md.data_other['Time'], alpha=0.9,\
                color=color[md.experiment][0], linewidth=2)
            x, y = [], []
            for d in md.d:
                x += [[md.data_other['Time']]*len(d.y)]
                y += [log10(d.y*d.mult_factors/d.Sc + 1)]
            plot( x, y, ' ', marker=color[md.experiment][1], markersize=4,\
                 color=color[md.experiment][0])
    ylabel(md.unit)
    xlabel(r"$min$")
    savefig("Images/All.png")


#### CBE Data: define the colors for the activation threshold plot for each experiment
color = {'21CTemp':['blue', '.', r"$21 ^oC$"],'RoomTemp':['blue', '.', r"$22 ^oC$"], '65CTemp':['yellow', 'v', r"$65 ^oC$"],\
            '70CTemp':['orange', '^', r"$70 ^oC$"], '75CTemp':['red', '<', r"$75 ^oC$"],\
            '80CTemp':['darkred', '>', r"$80 ^oC$"]} 

def PlotTemp2( md_list, th=2):
    """ Second function (and main) to analyze the results of the activation threshold"""
    hold = {} #empty dictionary
    for md in md_list:
        hold[md.experiment] = [[0],[0]]+color[md.experiment]    
    for md in md_list:
            p = md.CalcProb(th=th)
            hold[md.experiment][0] += [md.data_other['Time']]
            hold[md.experiment][1] += [p]

    ### This 'plot' is only to get the label right for room temp.
    plot( 0, 0, '-', marker='.', markersize=1, color='blue', linewidth=2)
    for experiment in hold:
        plot( hold[experiment][0], hold[experiment][1], '-',\
             marker=hold[experiment][3], markersize=5,\
                color=hold[experiment][2], linewidth=2, label=hold[experiment][4])

    ylabel(r"$P( E <= %d | Y)$" % (th))
    xlabel(r"$min$")
    hlines( 0, -1, 91)
    vlines( 0, -0.1,1.1)

    xlim((-1,91))
    ylim((-0.1,1.1))
    axis([-1,20,-0.1,1.1])
    legend()
    return hold

    
def Main_Intra_All( CBEData, LabExps, fig=[0,1,2]):
    """ Main_Intra is a function to analize intra labs data.
    """

    md_list = []
    for l,experiment in enumerate(LabExps): #['RoomTemp', '65CTemp', '70CTemp', '75CTemp', '80CTemp']): 
        unique_times = CBEData[experiment]["Time"].unique()
        for time in unique_times:
            close(fig[0])
            close(fig[1])
            md_list += [AnaTimeData( experiment=experiment, time=time, CBEData=CBEData, T=50000, fig=fig[:2])]      
    #PlotTemp(md_list)
    return md_list






