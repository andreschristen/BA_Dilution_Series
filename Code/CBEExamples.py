# -*- coding: utf-8 -*-
"""
Created on Mon 2018.02.28:00:00:00

@author: Nigel Orlando Santillán Morales

Abstrac:
Examples of experiments from different labs.
"""

######################################################################################################################################################################################
from pickle import load, dump


from pandas import read_excel

from pylab import plot, rc, figure, close, savefig, xlim, ylabel, xlabel, tight_layout, legend
rc('font', size=18)

from CBEDataAnalysis import AnaTimeData, PlotIndPost, AnaBF

from DilExp import GetKDE

######################################################################################################################################################################################

"""  
EXAMPLES: Intra Lab Experiment.
"""



##################### CBE data #######################
### Three repetitions, ten plated drops in two petri dishes, default values.
# Experiment:       Times:
CBE =\
[['RoomTemp', [ 1, 10, 15], "blue"],\
 [ '65CTemp', [15, 30, 45, 60, 90], "yellow"],\
 [ '70CTemp', [10, 20, 30, 40, 60], "orange"],\
 [ '75CTemp', [10, 20], "red"],\
 [ '80CTemp', [ 1,  2], "firebrick"]]


def ExpPlots( CBEData, experiment, time, plot_ind=False, K=[0,1], betabinom=False):
    md = AnaTimeData( CBEData=CBEData, experiment=experiment, time=time, T=500000, fig=[0,1], betabinom=betabinom)
    figure(0)
    tight_layout()
    savefig("Images/%s_%dmin_TS.png" % (experiment, time))
    if plot_ind:
        if betabinom:
            PlotIndPost(md, fig=[1,2], K=K, color="green") #without hierachycal model
        else:
            PlotIndPost(md, fig=[1,2], K=K)
        figure(2)               
        tight_layout()
        savefig("Images/%s_%dmin_IndPosts.jpg" % (experiment, time))
    figure(1)
    tight_layout()
    if betabinom:
        savefig("Images/%s_%dmin_bbinom_Results.jpg" % (experiment, time))
    else:
        savefig("Images/%s_%dmin_Results.jpg" % (experiment, time))
    return md


if __name__ == '__main__':
    CBEData = read_excel( 'CBE_BiofilmHotWaterStudies.xls',\
                ['RoomTemp', '65CTemp', '70CTemp', '75CTemp', '80CTemp'])
    ### Data taken from spreadsheet: "Biofilm Hot Water Studies.xlsx"
    ### NOTE: The RoomTemp experiment corresponds to the control of the 80C experiment. 

    ### md = ExpPlots(...) is a MultiDilExp object, see DilExp
    ### md.d[0].y is the data for repetition 0 etc.
    ### The simulated values of E are available in md.TwalkE()

     ### Run all data BF's with beta-binomial
    rt_CBE =[]
    print("%16s, %2s, %16s, %16s, %16s, %16s" %\
              ( "experiment", "k" , "bbinom" , "binom" , "Prob binom", "BF"))
    for item in CBE:
        experiment= item[0]
        for time in item[1]:
            rt_CBE += AnaBF( CBEData=CBEData, experiment=experiment, time=time)
    dump( rt_CBE, open("CBE_rt_CBE.pkl", "wb")) #To be used by InterLabExamples in the BF plot

    All = {} ### Dictiionary to hold MCMC iterations of E for all experiments
    
    experiment='70CTemp'
    time=10
    close(1)
    close(2)
    md = ExpPlots( CBEData=CBEData, experiment=experiment, time=time,\
                  betabinom=False)#, plot_ind=True, K=[0,1,2])
    #md_bb = ExpPlots( CBEData=CBEData, experiment=experiment, time=time, betabinom=True)#, plot_ind=True, K=[0,1,2])
    All[experiment] = {}
    All[experiment][time.__str__() + 'min'] = md.TwalkE()

    ######  Figs 3 and 4 take some 5 min to run

    ### Plots for Fig. 3
    experiment='65CTemp'
    time=15
    close(1)
    close(2)
    md = ExpPlots( CBEData=CBEData, experiment=experiment, time=time)
    All[experiment] = {}
    All[experiment][time.__str__() + 'min'] = md.TwalkE()

    experiment='75CTemp'
    time=10
    close(1)
    close(2)
    md = ExpPlots( CBEData=CBEData, experiment=experiment, time=time)
    All[experiment] = {}
    All[experiment][time.__str__() + 'min'] = md.TwalkE()

    ### This last one is also neede for fig. 4
    experiment='RoomTemp'
    time=15
    close(1)
    close(2)
    mdControl = ExpPlots( CBEData=CBEData, experiment=experiment, time=time)
    All[experiment] = {}
    All[experiment][time.__str__() + 'min'] = mdControl.TwalkE()


    ### Plots for Fig. 4
    experiment='80CTemp'
    time=2
    close(1)
    close(2)
    md = ExpPlots( CBEData=CBEData, experiment=experiment, time=time, plot_ind=True)
    All[experiment] = {}
    All[experiment][time.__str__() + 'min'] = md.TwalkE()

    figure(1)
    xlim((-0.5,3.5))
    tight_layout()
    savefig("Images/%s_%dmin_Results.jpg" % (experiment, time))
    figure(2)
    xlim((0,900))
    tight_layout()
    savefig("Images/%s_%dmin_IndPosts.jpg" % (experiment, time))
    ### Log reduction wrt RoomTemp experiment
    close(1)
    figure(1)
    LR = mdControl.TwalkE() - md.TwalkE()
    e, kde = GetKDE( LR, alpha=0.0000001)
    plot( e, kde, 'k-')
    ylabel("Density")
    xlabel(r"$log_{10}\left(\frac{CFU_0 + 1}{CFU + 1}\right)$")
    xlim(( 4, 9))
    tight_layout()
    savefig("Images/%s_%dmin_LR.jpg" % (experiment, time))
    print("%s_%s, $P[ LR > 3 ] = %6.4f" % (experiment, time, sum(LR > 3)/len(LR)))

    ### Activation threshold figure: Takes longer, some 15 min
    ### We load the data from "CBEallE.pkl" below ########
    
    ### Remaining experiments:
    CBE_R=[\
         ['RoomTemp', [ 10, 15]],\
         [ '65CTemp', [ 30, 45, 60, 90]],\
         [ '70CTemp', [ 20, 30, 40, 60]],\
         [ '75CTemp', [ 20]],\
         [ '80CTemp', [  1]]]

    for ex in CBE_R:
        experiment = ex[0]
        for time in ex[1]:
            md = ExpPlots( CBEData=CBEData, experiment=experiment, time=time)
            All[experiment][time.__str__() + 'min'] = md.TwalkE()
            #print("All[%s][%s]" % ( experiment, time.__str__() + 'min'))
    dump( All, open("CBEallE.pkl", "wb"))
    CBE =\
    [['RoomTemp', [ 1, 10, 15], "blue"],\
     [ '65CTemp', [15, 30, 45, 60, 90], "yellow"],\
     [ '70CTemp', [10, 20, 30, 40, 60], "orange"],\
     [ '75CTemp', [10, 20], "red"],\
     [ '80CTemp', [ 1,  2], "firebrick"]]

    #All = load(open("CBEallE.pkl", "rb"))
    close(2)
    figure(2)
    eh = 2 ### Threshold for E
    for ex in CBE:
        experiment = ex[0]
        act_prob = []
        for time in ex[1]:
            act_prob += [sum(All[experiment][time.__str__() + 'min'] < eh)/len(All[experiment][time.__str__() + 'min'])]
        ex += [act_prob]
    for ex in CBE[1:]:
        experiment, time, color, act_prob = ex
        plot( time, act_prob, '-o', color=color, label=experiment[:2] + r" $^o$C")
    legend(fontsize=14) #(loc=( 61, 0.3))
    xlabel("min")
    ylabel("Act. Probability")
    tight_layout()
    savefig("Images/ActivationProbability.jpg")






