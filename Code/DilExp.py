#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 11:20:29 2018

@author: jac
"""

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

from scipy.special import gammaln, betaln
from scipy.stats import expon, gaussian_kde, binom, poisson
from scipy.stats.mstats import mquantiles

from numpy.random import choice
from numpy import append, array, where, arange, linspace, zeros, ceil
from numpy import log, log10, exp, var, sqrt, diff, cumsum, mean, inf

from pylab import plot, vlines, xlabel, ylabel, xlim, ylim
from pylab import figure

from pytwalk import pytwalk



"""
#### Example importing data from Excel file using ReadData
#### Read data, typical a single lab per file, in standarized format.
#### See for example DataInputExampleLab5.xls 
#### Only the first countable dilution is considered in this format.

from pandas import read_excel
from pylab import tight_layout

#### Read data:

spreadsheet = read_excel( 'DataInputExampleLab5.xls', ['Parameters', 'Data', 'Exp Info'])
print("Reading data from Lab5.")
Lab5 = ReadData(spreadsheet)

spreadsheet = read_excel( 'DataInputExampleLab6.xls', ['Parameters', 'Data', 'Exp Info'])    
print("\nReading data from Lab6.")
Lab6 = ReadData(spreadsheet)

### Acces individual MultiDilExp objects, ie. single experiment with repetitions.
Lab5['Expose_low1'].RunTwalk(T=100000)
figure(1)
Lab5['Expose_low1'].PlotResults()
tight_layout()

### Do a multilab analysis:
Lab5['CaviCide_high1'].RunTwalk(T=100000)
Lab6['CaviCide_high1'].RunTwalk(T=100000)
mdlist = [ Lab5['CaviCide_high1'], Lab6['CaviCide_high1']]
Il = InterLabDilExp(mdlist)
Il.RunTwalk(T=200000)  
figure(2)
kde_list = Il.AnaEs(col=['blue','blue'])

"""            
            

def GetKDE( data, alpha=0.001, N=200, cut_at_zero=True, bw_method=lambda kde: kde.n**(-1.0/((kde.d + 4)*1.6))):
    """ Kernel density estimation of data using guassian kernels.
        returns mesh using the quantiles alpha and 1-alpha, and its respective evaluation."""
    kde = gaussian_kde( data, bw_method=bw_method)
    q = mquantiles( data, [alpha,1-alpha])
    if (q[0] < 1.0) and cut_at_zero:
        q[0]=0.0
    e = linspace( q[0], q[1], N)
    kde = kde.evaluate(e)
    C = sum(diff(e)*kde[:-1]) ### Keep the integral == 1
    return e, kde/C

            
def ReadData(spreadsheet):
    """Read data, typical a single lab, in standarized format.
       See DataInputExample.xls .
       
       Note that the surface area is not read and is set to 1.
       It is not currently used.
    """
    ###  Read in all parameters
    All = {} #Dictionary
    for i, par in enumerate(spreadsheet['Parameters'].loc[:,'Parameter']):
        All[par] = spreadsheet['Parameters'].loc[i,'Value']
    All['Experiments'] = len(spreadsheet['Data']['Experiment'].unique())
    for exp_name in spreadsheet['Data']['Experiment'].unique():
        All[exp_name] =\
            MultiDilExp( K=All['K'], J=All['J'], alpha0=All['alpha0'], alpha=All['alpha'],\
            alphap=All['alphap'], c=All['c'], Sc=1, Sc_unit=All['Sc_unit'],\
            q=All['q'], M=All['M'], pr_a=10, b=500)
    m = All['m'] ## Number of drops or plates
    for row, exp_name in enumerate(spreadsheet['Data']['Experiment']):
        k = spreadsheet['Data']['Repetition'][row] - 1 #Repetition
        s = array([spreadsheet['Data']['Dilution'][row]]*m) #Dilution counted
        y = array([spreadsheet['Data'].loc[ row, "Count%d" % (c+1,)] for c in range(m)]) #Counts
        print("Experiment %s, loading repetition %d." % (exp_name,k), s, y)
        All[exp_name].Data( k, s, y)
    return All



class DilExp:
    
    def __init__( self, J=7, alpha0=1, alpha=10, alphap=1000, c=30, q=0.05,\
                 pr_a=10, Sc=1, Sc_unit=r"cm^2", betabinom=False):
        """Dilution CFU count experiment analysis:
            (default values for a typical drop plate method)
            J, number of tube dilutions j=0,1,...,J-1 (=7)
            alpha0, initial dilution factor for tube 1 from tube 0: eg. 1ml from 4*10ml
            alpha, dilution factor for each tube (=10)
            alphap, dilution factor for the drop or plated volumen (=1000)
            c, maximum count for the drop or plated volumen (=30)
            q, probibility of misscounting (=0.05)
            pr_a, prior mean of the count, given the dilution s to be counted (=20)
            Sc, coupon surface area (=4.13)
            Sc_unit, coupon surface area units (=r"cm^2")."""

        ### Design
        self.J = J
        self.alpha0 = alpha0
        self.alpha = alpha
        self.alphap = alphap
        self.c = c #
        self.q = q
        self.Sc = Sc/self.alpha0 ##Then all is standarized to all tubes with the same volumen.
        self.Sc_unit = Sc_unit

        ### Prior dist parameters, this is the distribution for n_s
        self.center = self.c/2
        self.pr_2a2 = 4*2*float(pr_a**2)
        self.prior_range = [ 0, int(round(self.c*2.0))]
        self.l_range = arange( self.prior_range[0], self.prior_range[1], dtype=int)  ### Range for n_s  
        self.betabinom = betabinom


    def Loghs( self, i):  
        """Log Prior for n_s. Uniform."""
        #if (isinstance( i, int)): #does not work with self.l_range[i] !
        if ((i >= self.prior_range[0]) and (i <= self.prior_range[-1])):
            return 1.0
        else:
            return -inf
        #else:
        #    return -inf


    def PlotPriorhs(self):
        y = exp(array([self.Loghs(i) for i in self.l_range]))
        y /= sum(y)
        plot( self.l_range, y, 'go')
        vlines( self.l_range, 0, y,  colors='k', linestyles='-')
        xlabel(r"$n_s$")
        ylabel("Probability")

    
    def Logh0( self, k):
        """LogPrior for n_0, given the dilution to be counted self.s."""
        return sum([self.Loghs(int(round(k*self.dil[i]))) \
                    for i in range(self.m)]) 

    #def Logh0( self, k):
    #    """LogPrior for n_0, given the dilution to be counted self.s."""
    #    return sum([-log(k + 1) for i in range(self.m)])  


    
    def PlotPriorh0( self, post_rng_indx):
        y = exp(array([self.Logh0(k) for k in post_rng_indx]))
        y /= sum(y)
        if (self.s > 0):
            plot( post_rng_indx, y, 'g-')
        else:
            plot( post_rng_indx, y, 'go')
            vlines( post_rng_indx, 0, y,  colors='k', linestyles='-')
        xlabel(r"$n_0$")
        ylabel("Probability")

    def SimData( self, N0, s):
        """Simulate data with true N0 and
            s = a list with the decreasing number of the counted dilution for each plate.
        """
        self.dil = float(self.alpha)**(-s) * float(self.alphap)**-1 * (1.0-self.q)
        return binom.rvs( n=N0, p=self.dil)

    def Data( self, s, y, calc_post=False):
        """s = a list with the decreasing number of the counted dilution for each plate
           y = a list with the counts, same length as s."""
        self.s = array(s, dtype=int)
        self.y = array(y, dtype=int)
        #print "s",s
        #print "y",y
        self.m = len(self.y)
        ### Some constants:
        ### mult_factors are the blow up factors that is
        ### approximately one expects self.y[i]*self.mult_factors[i] CFUs
        self.mult_factors = self.alpha0 * self.alpha**self.s * self.alphap ### Vector
        if all(self.y == self.c):
            self.mul = max(self.mult_factors)
        else:
            if all(self.y >= self.c):
                jk = 0 ### All TNTC, choose any count
            else:
                jk = where(self.y < self.c)[0][0] ##First dilution with a countable CFU
            self.mul = self.mult_factors[jk] #max(self.mult_factors)
        ### This is the binomial probability of each y[i] | n_0, q ~ Bi( n_0, dil[i])
        self.dil = float(self.alpha0)**-1 * float(self.alpha)**(-self.s) *\
                         float(self.alphap)**-1 * (1.0-self.q)
        #### ### Vector  min(dil) most be > 10**-15
        self.ymax = max(self.y)
        ### These constants are auxiliary to the calculations 
        self.C1 = log(self.dil) ### Vector
        self.C2 = log(1.0-self.dil) ### Vector
        self.C3 = sum(self.y*(self.C1-self.C2) - gammaln(self.y + 1))
        self.C2 = sum(self.C2)
        
        self.mn = mean(self.y)
        self.sd = sqrt(var(self.y))
        if calc_post:
            tmp = self.betabinom
            self.CalculatePost(betabinom=not(tmp))
            self.CalculatePost(betabinom=tmp)
            self.BF = exp(self.k_binom_c - self.k_betabinom_c)*self.k_betabinom/self.k_binom

    def LogLikelihood( self, k):
        """We dont use binom.logpmf to improve speed, instead we calculate the logpmf
           mannulay, without some constants.
           LogLikelihoodCensored uses the logpmf.
        """
        return self.m*gammaln(k+1) + k*self.C2 - sum(gammaln(k-self.y + 1)) + self.C3

    def LogLikelihoodCensored( self, k):
        """If we are dealing with saturated data, that is y=de.c
           then we will need to subtutute the above likelihood with this one."""
        llkh = 0.0
        for j, y in enumerate(self.y):
            if y < self.c: ###Not censored
                llkh += binom.logpmf( k=y, n=k, p=self.dil[j])
            else: ### Censored: likelihood prod_i^m P(Y_i > c | k)
                llkh += log(sum(binom.pmf( k=arange(self.c,2*self.prior_range[-1]), n=k, p=self.dil[j])))
        return llkh

    def LogLikelihoodBetaBinom( self, k):
        """This is a more general likelihood, a BetaBinomial likeliohood."""
        llkh = 0.0
        for j, y in enumerate(self.y):
            #n0 = ceil(1.0/(1-self.dil[j]))
            #la = 2 #without the restriction of la s* > 1
            la = 1/self.dil[j] + 1
            a = la*self.dil[j] #= 2 p=self.dil[j]
            b = la*(1-self.dil[j])
            # BetaBinom( y | k, a, b) 
            llkh += gammaln(k+1) - gammaln(y+1) - gammaln(k-y+1) +\
                        betaln( y+a, k-y+b) - betaln( a, b)
        return llkh
 
    def LogPost( self, k):
        if self.betabinom:
            return self.Logh0(k) + self.LogLikelihoodBetaBinom(k)
        else:
            return self.Logh0(k) + self.LogLikelihood(k)


    def CalculatePost( self, betabinom=False, k_range=0):    
        """ Calculate and normlaize posterior.
            k_range=0, automatic range calculation for k."""
        ### Range for k:
        if isinstance(k_range, int):
            if (max(self.s) == 0):  ### Range for n_0
                self.k_range = arange( max( self.ymax, self.l_range[0]*self.mul), self.l_range[-1]*self.mul, dtype=int)   
            else:
                self.k_range = arange( max( self.ymax, self.l_range[0]*self.mul), self.l_range[-1]*self.mul,\
                    self.alpha**(max(self.s)-1) * self.alphap, dtype=int)
        else:
            self.k_range = k_range

        self.betabinom = betabinom
        self.post_array = array([self.LogPost( k ) for k in self.k_range])
        c = 200 - max(self.post_array)
        self.post_array += c ### Set the maximum to 200 for numerical stability
        self.post_array = exp(self.post_array)
        self.K = 1.0/sum(self.post_array)
        if betabinom:
            self.k_betabinom = 1/self.K
            self.k_betabinom_c = c
            ### For model comparison puerpuses
            ### the actual normalization constant is exp(-self.k_betabinom_c)*self.k_betabinom
        else:
            self.k_binom = 1/self.K
            self.k_binom_c = c
            ### For model comparison puerpuses
            ### the actual normalization constant is exp(-self.k_binom_c)*self.k_binom
        self.post_array *= self.K
        self.cdf = cumsum(self.post_array)


    def SimPost( self, size, CFU_X_Area=True, log10_p1=True):
        ### Post. dist previously calculated
        if CFU_X_Area:
            Sc = float(self.Sc)
        else:
            Sc = 1

        if log10_p1:
            x  = log10((self.k_range +1)/Sc)
        else:
            x = self.k_range/Sc

        return choice( a=x, size=size, replace=True, p=self.post_array)


    def Post( self, k):
        return self.K*exp(self.LogPost(k))


    def PlotPost( self, betabinom=False, CFU_X_Area=True, log10_p1=True, plot_data=True, color='blue', linestyle='solid'):
        self.CalculatePost(betabinom=betabinom)
        
        ### Select range for posterior:
        #S = 0.0
        q1 = 0.0005  ### Range within these quatiles
        q2 = 1-q1
        j1 = 0
        j2 = len(self.k_range)-1
        for j,p in enumerate(self.cdf):
            if (p <= q1):
                j1 = j
            if (p > q2):
                j2 = j
                break
            #S += p
        post_rng_indx = arange( j1, j2+1, dtype=int)

        if CFU_X_Area:
            if self.Sc == 1:
                unit = r"CFU + 1"
                Sc = 1
            else:
                Sc = float(self.Sc)
                unit = r"\frac{CFU + 1}{%s}" % self.Sc_unit
        else:
            Sc = 1
            unit = r"CFU"

        if log10_p1:
            unit = r"$log_{10}\left(%s\right)$" % (unit,) 
            tmp  = log10((self.k_range[post_rng_indx]+1)/Sc)
            x = tmp[:-1]
            fx = self.post_array[post_rng_indx[:-1]]/diff(tmp)
            ylab = "Density"
        else:
            unit = r"$%s$" % (unit,)
            x = self.k_range[post_rng_indx]/Sc
            fx = self.post_array[post_rng_indx]
            ylab = "Probability"
        
        ### Plot posterior:
        if ((max(self.s) > 0) or log10_p1):
            plot( x, fx, '-', color=color, linestyle=linestyle)  ###Cont. approximation or density
        else:
            plot( x, fx, '.', color=color) ### Discrete
            #vlines( phi(self.k_range[post_rng_indx]), 0, self.post_array[post_rng_indx],  colors='k', linestyles='-')
        #self.PlotPriorh0(post_rng_indx)

        ### Plot blown up data
        if (not(log10_p1) and plot_data):
            plot( self.y*self.mul, [self.post_array[post_rng_indx[0]]]*self.m, '*', color='orange')
            ### Plot mean+- 2 sd of blown up data 
            plot( array([self.mn-3*self.sd, self.mn, self.mn+3*self.sd])*self.mul,\
                 [self.post_array[post_rng_indx[0]]]*3, 'go', linewidth=1)
            plot( array([self.mn-self.sd, self.mn, self.mn+self.sd])*self.mul,\
                 [self.post_array[post_rng_indx[0]]]*3, 'go', linewidth=1)
            #plot( array([self.mn-3*self.sd, self.mn, self.mn+3*self.sd])*self.mul/Sc,\
                  #    [self.post_array[post_rng_indx[0]]]*3, 'r*')
        xlabel(unit)
        ylabel(ylab)
        #title("Posterior, data J=%d, y=[%d,%d]." % (J, y[0], y[1]))
        self.x = x
        self.fx = fx ###Save the last plotted values
        return xlim() + ylim()

  
class MultiDilExp:
    
    def __init__( self, K=3, b=500, M=10, J=7, alpha0=1, alpha=10, alphap=1000,\
                 c=30, q=0.05, pr_a=10, Sc=1, Sc_unit=r"cm^2", betabinom=False):
        """Dilution CFU count experiment analysis on several tubes:
            (default values for a typical drop plate method)
            K, number of tubes (=3)
            b, prior parameter for the exponential distribution of par. A (=500)
            M, upper value for the U(0,M) prior for par. E. (=10)
            
            Parameters for each plate:
            J, number of tube dilutions (=7)
            alpha0, dilution factor for tube 1 from tube 0: eg. 1ml from 4*10ml
            alpha, dilution factor for each tube (=10)
            alphap, dilution factor for the drop or plated volumen (=1000)
            c, maximum count for the drop or plated volumen (=30)
            q, probibility of misscounting (=0.05)
            pr_a, prior mean of the count, given the dilution s to be counted (=20)
            Sc, coupon surface area (=4.13)
            Sc_unit, coupon surface area units (=r"cm^2")."""

        self.K = K
        self.n = self.K + 2 # Number of parameters, A,E + S's
        self.b = float(b)
        self.Sc = Sc
        self.M = float(M)
        self.MinE = log10(1/self.Sc) #Minimum value for E and S's: CFU=0
        self.Cte = -log(self.M) - log(b) 
        ### Individual plates, assume single Sc for all coupons
        self.d = [DilExp(J=J, alpha0=alpha0, alpha=alpha, alphap=alphap, c=c, q=q, pr_a=pr_a, Sc=Sc, Sc_unit=Sc_unit, betabinom=betabinom)\
                      for k in range(self.K)]
        if Sc == 1:
            self.unit = r"$log_{10}\left( CFU + 1 \right)$"
        else:
            self.unit = r"$log_{10}\left(\frac{CFU + 1}{%s} \right)$" % (Sc_unit,)
        self.logten = log(10.0)


    def Data(self, k, s, y):
        """Enter data for DilExp (repetition) k."""
        self.d[k].Data( s=s, y=y, calc_post=True) #Initialize


    def LogPost( self, a, e, s):
        """LogLikelihood, E=e, A=a and s array of log10((N_0^k + 1)/Sc) values."""
        n0 = self.Sc*10**s - 1 #Obtain the n0's,
        if self.d[0].betabinom: #All are betbinom ll or not
            lp = sum([self.d[k].LogLikelihoodBetaBinom(int(n0[k])) for k in range(self.K)])
        else:
            lp = sum([self.d[k].LogLikelihood(int(n0[k])) for k in range(self.K)])
        ### s = log10(n0*d[0].Sc + 1)
        ### Here the dist. of N^k_0 given E=e, A=a
        ### -self.logten is the jacobian correction
        gashape = a #e**2 * a ### Shape and scale parameters for the Gamma
        gasc = e/a #1/(e*a)     ### Hierarchycal prior of the s's
        lp += self.K*(-gashape*log(gasc) - gammaln(gashape)) +\
            sum([(gashape-1)*log(sk) - sk/gasc  for sk in s]) #+ sk*self.logten
        ### HERE THE HYPER PRIOR FOR E ~ U(0,10) (constant), AND A ~exp(b)
        lp += self.Cte  - a/self.b#- (2-1)*log(a) - (3-1)*log(a)  - e/self.M a/self.b# - a/self.b  #- log(a) 
        return lp


    def SimInitValues(self):
        """Simulate initial values for A, E and S (=log10(N/Sc+1))."""
        
        #First, simuate values from the free post of each n0
        S = array([d.SimPost( size=1, CFU_X_Area=True, log10_p1=True)\
                  for d in self.d]).flatten()
        # Take E as its mean
        E = array(mean(S))
        # And A from its prior
        A = expon.rvs(size=1, scale=self.b)
        return append(append( A, E), S)
    
    
    def Energy(self, x):
        return -self.LogPost( x[0], x[1], x[2:])


    def Supp(self, x):
        ###        A                  E and S's                E and S's
        return (all(x > 0.0) and all(self.MinE < x[1:]) and all(x[1:] <= self.M))
    
    
    def RunTwalk( self, T):
        """Run the twalk MCMC with Tr iterations, simulate the two
           initial vaues for theta from self.SimInit."""
        x0=self.SimInitValues()
        xp0=self.SimInitValues()
        while any(abs(x0 -xp0) <= 0.001): # with 0 is the test in pytwalk
            sleep(0.500)
            xp0=self.SimInitValues()
        self.twalk = pytwalk( n=self.n, U=self.Energy, Supp=self.Supp) #Open the twalk object
        self.twalk.Run( T=T, x0=x0, xp0=xp0 )
        
        
    def CalcProb(self, th=2.0):
        """Calculate P( E <= th | Data)."""
        return sum(self.twalk.Output[:,1] <= th)/float(self.twalk.T)


    def TwalkOutput(self):
        """Results of the twalk."""
        return self.twalk.Output
    
    def TwalkE(self):
        """twalk simulation for parameter E."""
        return self.twalk.Output[:,1]


    def PlotvPar(self, par, position, alpha=1.0, color='k', linewidth=2):
        kde = gaussian_kde(self.twalk.Output[:,1])
        e = linspace(0, self.M, 100)
        post = kde.evaluate(e)
        cdf = cumsum(post[:-1]*diff(e))
        i1 = 0
        i2 = len(cdf)
        for i in range(len(cdf)):
            if cdf[i] <= 0.05:
                i1=i
            if cdf[-i] >= 0.95:
                i2=-i
        e = e[i1:i2]
        post = post[i1:i2]
        plot( position+alpha*post, e, '-', color=color, linewidth=linewidth)
        plot( position-alpha*post, e, '-', color=color, linewidth=linewidth)


    def PlotResults(self, plot_data=True, fignum=0, color="blue"):        
        # Plot marg. post for E
        print("Ploting Results ...""")

        ### Plot marg. post for S's, count data and mean
        self.mean_s = []
        print("Data Mean:\n tube: mean")
        hold_data = array([])
        figure(fignum)
        for k in range(self.K):
            e, pdf = GetKDE(self.twalk.Output[:,2+k])
            plot( e, pdf, '-', color=color)
            sk_data = log10(self.d[k].y*self.d[k].mult_factors/self.d[k].Sc +1)
            if plot_data:
                hold_data = append(hold_data, sk_data)
            tmp = mean(self.d[k].y*self.d[k].mult_factors)
            self.mean_s += [log10(tmp/self.d[k].Sc +1)]
            #plot([self.mean_s[k]], [0], 'r*', markersize=7)
            print((" %4d: %f  " % (k,self.mean_s[k])))

        ### Plot marginal for E, commonly there is no need for burn-in
        e, pdf = GetKDE(self.twalk.Output[:,1])
        if color != "blue":
            plot( e, pdf, 'k--', linewidth=1)
        else:
            plot( e, pdf, 'k-', linewidth=2)

        ### Axis labels
        xlabel(self.unit)
        ylabel("Density")

        ### Create and plot classical estimates
        self.mean_s = array(self.mean_s)
        self.mean_s_mean = mean(self.mean_s)
        self.mean_s_sd = sqrt(var(self.mean_s))
        plot( [self.mean_s_mean - 3*self.mean_s_sd,self.mean_s_mean + 3*self.mean_s_sd], [0,0], 'r-', linewidth=4)
        plot( [self.mean_s_mean - 2*self.mean_s_sd,self.mean_s_mean + 2*self.mean_s_sd], [0,0], 'k|', markersize=6)
        plot( [self.mean_s_mean - 1*self.mean_s_sd,self.mean_s_mean + 1*self.mean_s_sd], [0,0], 'k|', markersize=6)
        plot( [self.mean_s_mean], [0], 'b^', markersize=6)        
        print((r"Classical estimate: %f$\pm$%f" % (self.mean_s_mean,self.mean_s_sd)))
        print(("Interval 1sd: [ %f, %f]" % (self.mean_s_mean - 1*self.mean_s_sd,self.mean_s_mean + 1*self.mean_s_sd)))
        print(("Interval 2sd: [ %f, %f]" % (self.mean_s_mean - 2*self.mean_s_sd,self.mean_s_mean + 2*self.mean_s_sd)))
        print(("Interval 3sd: [ %f, %f]" % (self.mean_s_mean - 3*self.mean_s_sd,self.mean_s_mean + 3*self.mean_s_sd)))
        ym, yM = ylim()
        ylim(( -0.05*(yM-ym), yM))
        vlines( 0, -0.05*(yM-ym), yM)
        plot( hold_data, [-0.025*(yM-ym)]*len(hold_data), 'k*', markersize=3)



class InterLabDilExp:
    
    def __init__( self, mdlist, M=10, b=500):
        """Dilution CFU count experiment analysis on several labs:
            mdlist: a list of MultiDilExp instances, representing each lab.
            We assume the twalk has been run on each MultiDilExp instance.
            
            A hierarchical model will be fitted to all E's, which are already
            normalized for Sc.  All labs shoould used the same unit for Sc, but
            may in fact used different Sc's.
            
            E_i | \mathcal{E} = e, \mathcal{A} = a \sim E=e, A=a \sim Ga( a, e a^{-1} ).
        """
        
        self.mdlist = mdlist
        self.L = len(self.mdlist) #Number of labs
        self.M = M
        ###Number of parameters in each md isntance.  First tow will be assigned
        ### E and A, the rest mdn-2 are the s para emeters, log10((N_0^k + 1)/Sc). 
        self.parsnum = [md.n for md in mdlist] 
        ###The fist two will be \mathcal{A} and \mathcal{E} global interlab pars.
        self.n = 2 + sum(self.parsnum) ## total number of parameters
        self.b = float(b)
        self.Cte = -log(self.M) - log(b)
        self.twalk = pytwalk( n=self.n, U=self.Energy, Supp=self.Supp) #Open the twalk object


    def Energy( self, x):
        """Un pack all parameters from selfn+2 vector x and send then to individual
           md LogPost instances."""
        lp = 0.0 #Log posterior
        ca = x[0]
        ce = x[1]
        indx = 2
        for i, md in enumerate(self.mdlist):
            ################### a         e             s
            lp += md.LogPost( x[indx], x[indx+1], x[indx+2:(indx+md.n)])
            ei = x[indx+1]
            ##Hyper prior for E_i
            lp += (ca-1)*log(ei)-ei*(ca/ce) 
            indx = indx+md.n
        ### The normalization constant in each hyper prior for E_i
        lp += self.L*(ca*log(ca) - ca*log(ce) - gammaln(ca))
        ### HERE THE HYPER-HYPER PRIOR FOR E ~ U(0,10) (constant), AND A ~exp(b)
        lp += self.Cte - ca/self.b        
        return -lp #Return energy
            
    def Supp( self, x):
        """Un pack all parameters from selfn+2 vector x and send then to individual
           md Supp instances."""
        if not(all(x[0:2] > 0)):
            return False ###E or A not positive
        indx = 2
        for i, md in enumerate(self.mdlist):
            ################### a         e             s
            if not(md.Supp( x[indx:(indx+md.n)])):
                return False
            #self.E[i] = x[indx+1]
            indx = indx+md.n
        return True


    def SimInitValues(self):
        """Simulate initial values, from the corresponding function in each
           md instance, adding su¡imulated values for \mathcal{A} and \mathcal{E}.
        """
        
        x = zeros(self.n)
        meane = 0.0
        indx = 2
        for i, md in enumerate(self.mdlist):
            x[indx:(indx+md.n)] = md.SimInitValues()
            meane += x[indx+1]
            indx = indx+md.n
        # Take E as the mean of the simulated E's
        ce = array([meane/self.L])
        # And A from its prior
        ca = expon.rvs(size=1, scale=self.b)
        return append(append( ca, ce), x[2:])
        

    def RunTwalk( self, T):
        """Run the twalk MCMC with Tr iterations, simulate the two
           initial vaues for theta from self.SimInit."""
        x0=self.SimInitValues()
        xp0=self.SimInitValues()
        while any(abs(x0 -xp0) <= 0.001): # with 0 is the test in pytwalk
            sleep(0.005)
            xp0=self.SimInitValues()
        self.twalk.Run( T=T, x0=x0, xp0=xp0 )
    
    def TwalkcE(self):
        return self.twalk.Output[:,1]

        
    def AnaEs( self, col, calc_hierachical_Es=False, alpha=0.001, N=200):
        """col: list of colosr to use to plot individual lab post. E's ."""
        kde_list = []
        e0=15
        e1=0
        for i,md in enumerate(self.mdlist):    
            q = mquantiles( md.TwalkE(), [alpha,1-alpha])
            if (q[0] < 1.0): #cut_at_zero: Positive densities
                q[0]=0.0
            e0 = min(q[0],e0)
            e1 = max(q[1],e1)
        e = linspace( e0, e1, N)

        for i,md in enumerate(self.mdlist):
            print("Calculating kde for free E_%d ..." % (i,))
            kde_instance = gaussian_kde( md.TwalkE(), bw_method=lambda kde: kde.n**(-1.0/((kde.d + 4)*1.6)))
            kde = kde_instance.evaluate(e)
            C = sum(diff(e)*kde[:-1]) ### Keep the integral == 1
            kde /= C
            kde_list += [[kde_instance,e,kde]]
            plot( e, kde, '-', color=col[i], linewidth=1)
        if calc_hierachical_Es:
            for i in range(self.L):
                print("Calculating kde for hierachical E_%d ..." % (i,))
                kde_instance = gaussian_kde( self.twalk.Output[:,2+i*self.mdlist[i].n+1], bw_method=lambda kde: kde.n**(-1.0/((kde.d + 4)*1.6)))
                kde = kde_instance.evaluate(e)
                C = sum(diff(e)*kde[:-1]) ### Keep the integral == 1
                kde /= C
                kde_list += [[kde_instance,e,kde]]
                plot( e, kde, '-', color=col[i], linewidth=1)
        print("Calculating kde of \mathcal{E}_%d ..." % (i,))
        kde_instance = gaussian_kde( self.TwalkcE(), bw_method=lambda kde: kde.n**(-1.0/((kde.d + 4)*1.6)))
        kde = kde_instance.evaluate(e)
        C = sum(diff(e)*kde[:-1]) ### Keep the integral == 1
        kde /= C
        kde_list += [[kde_instance,e,kde]]    
        plot( e, kde, 'k-', linewidth=2)
        xlabel(r"$log_{10}\left( CFU + 1 \right)$")
        ylabel("Density")
        """
        TVs = zeros((len(self.mdlist),len(self.mdlist)))
        for i,md1 in enumerate(self.mdlist):
            for j,md2 in enumerate(self.mdlist):
                print("Calculating ||f_%d - f_%d||_{TV} ..." % (i,j))
                if i==j:
                    TVs[i,j] = 0.0
                else:
                    tv_eval = abs(kde_list[i][2]-kde_list[j][2]) 
                    TVs[i,j] = 0.5*sum(diff(kde_list[i][1])*tv_eval[:-1])
        """
        return kde_list

            