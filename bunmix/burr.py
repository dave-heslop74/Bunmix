import numpy as np
from scipy.optimize import minimize
import pymc3 as pm
from pymc3.math import erf, sqrt, log, minimum, abs_, sgn, exp, clip, maximum, le
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.interpolate import interp1d
import matplotlib.ticker as mtick


##### Approximations of Burr distribution with specific mean, std and skewness for THENO #####
def Ghat(skew):
    
    #estimate G given skewness
    pG = np.array([-0.18519617,-0.36581369,-0.26346167,-0.10462069,-0.34142823,0.8666556])
    
    A = 0
    for i in range(6):
        A += pG[i]*skew**(5-i)
    
    B = -0.848
    
    return (A*B+1)**(1/B)

def Mhat(G):
    #estimate mean given G
    pM = np.array([-2.42043861e-05,7.07715287e-04,-8.72589408e-03,5.78173631e-02,-2.12287390e-01,3.65070722e-01,2.39311397e-02,2.39178197e-01])
    
    tot = 0
    for i in range(8):
        tot += pM[i]*log(G)**(7-i)
        
    return tot

def Vhat(G):
    
    #estimate variance given G
    pV = np.array([-2.33501964e-06,8.10258820e-05,-1.21332152e-03,1.01946967e-02,-5.21336282e-02,1.63404382e-01,-2.95110778e-01,2.50969395e-01,-3.73368931e-02])
    
    tot = 0
    for i in range(9):
        tot += pV[i]*log(G)**(8-i)
        
    return tot

def Bcdf(x0,mu,sigma,skew):
    
    A = 6.5
    G = Ghat(-abs_(skew))
    BM = Mhat(G)
    BV = Vhat(G)
    
    xs=(x0-mu)/sigma
    x = xs*sqrt(BV)+BM
    
    #print(sgn(skew).eval()==-1.0)

    #if sgn(skew)==-1.0:
    if le(skew,0.0):
        xc = clip(x,0,5)
        y = 1 - (1+xc**G)**-A
        #y[x<0] = 0
        #print('Bcdf: step1')
    else:
        #print('Bcdf: step2')
        x = -x+2*BM
        xc = clip(x,0,5)
        y = (1+xc**G)**-A        
    return y #, xc

##### Approximations of Burr distribution with specific mean, std and skewness for NUMPY #####
def npBcdf(x0,mu,sigma,skew):    
    
    A = 6.5
    G = npGhat(-np.abs(skew))
    BM = npMhat(G)
    BV = npVhat(G)
    
    xs=(x0-mu)/sigma
    x = xs*np.sqrt(BV)+BM
    
    if np.sign(skew)==-1:
    #if skew<=0.0:
        xc = np.clip(x,0,5)
        y = 1 - (1+xc**G)**-A
        #print('npBcdf: step1')
    else:
        #print('npBcdf: step2')
        x = -x+2*BM
        xc = np.clip(x,0,5)
        y = (1+xc**G)**-A        
    
    return y #, xc

def npMhat(G):
    #estimate mean given G
    pM = np.array([-2.42043861e-05,7.07715287e-04,-8.72589408e-03,5.78173631e-02,-2.12287390e-01,3.65070722e-01,2.39311397e-02,2.39178197e-01])
    
    tot = 0
    for i in range(8):
        tot += pM[i]*np.log(G)**(7-i)
        
    return tot

def npVhat(G):
    
    #estimate variance given G
    pV = np.array([-2.33501964e-06,8.10258820e-05,-1.21332152e-03,1.01946967e-02,-5.21336282e-02,1.63404382e-01,-2.95110778e-01,2.50969395e-01,-3.73368931e-02])
    
    tot = 0
    for i in range(9):
        tot += pV[i]*np.log(G)**(8-i)
        
    return tot

def npGhat(skew):
    
    #estimate G given skewness
    pG = np.array([-0.18519617,-0.36581369,-0.26346167,-0.10462069,-0.34142823,0.8666556])
    
    A = 0
    for i in range(6):
        A += pG[i]*skew**(5-i)
    
    B = -0.848
    
    return (A*B+1)**(1/B)

def npBpdf(x0,mu,sigma,skew):
    
    A = 6.5
    G = npGhat(-np.abs(skew))
    BM = npMhat(G)
    BV = npVhat(G)
    
    xs=(x0-mu)/sigma
    x = xs*np.sqrt(BV)+BM
    
    if skew>0:
        x = -x+2*BM
    
    xc = np.clip(x,0,5)
    y = A*G*(xc)**(G-1)*(1+(xc)**G)**(-A-1)*np.sqrt(BV)

    
    return y/sigma

## Objective function for MLE
def Bsse(x,x0,y0,Ncomp):

    hat = np.zeros(len(x0))
    for i in range(Ncomp):
        hat += npBcdf(x0,x[i*4],x[i*4+1],x[i*4+2])*x[i*4+3]
  
    return np.sum((hat-y0)**2)

## Find MLE for model complexities and return BIC-based probability
def Burr_AIC(x0,y0,Nmax):

    n = len(x0)
    bnds = ((0.0,0.9*np.max(x0)),(0.08,0.6),(-0.9,0.0),(0.05,1.5))
    
    Ncomp0 = np.arange(Nmax)+1 #number of components to compare
    AIC = np.zeros(Nmax)
    s2 = np.zeros(Nmax)

    xmle = []
    for Ncomp in Ncomp0:
        #Setup bounds
        bounds = bnds
        for i in range(Ncomp-1):
            bounds += bnds 
    
        # setup minimization iterations    
        Niter = 50
        sse = np.zeros(Niter)
        xout = []
        for j in range(Niter):
            
            #generate random initial state
            x = np.zeros(Ncomp*4)
            
            for i in range(Ncomp):
                x[i*4] = np.random.rand(1)*np.diff(bnds[0])+bnds[0][0]
                x[i*4+1] = np.random.rand(1)*np.diff(bnds[1])+bnds[1][0]
                x[i*4+2] = np.random.rand(1)*np.diff(bnds[2])+bnds[2][0]
                x[i*4+3] = np.random.rand(1)*np.diff(bnds[3])+bnds[3][0]    
        
            res = minimize(Bsse,x,args=(x0,y0,Ncomp),bounds=bounds)
            sse[j] = res.fun
            xout.append(res.x)
            

        idx = np.argmin(sse)
        xmle.append(xout[idx])
        s2[Ncomp-1] = np.min(sse)/n
        k = Ncomp*4+1

        AIC[Ncomp-1] = 2*k + n*np.log(np.min(sse)) + 2*k*(k+1)/(n-k-1)
        
    DAIC = AIC-np.min(AIC)
    W = np.exp(-0.5*DAIC)/np.sum(np.exp(-0.5*DAIC))
    
    return W, s2, xmle

#### Functions for PyMC ####

def mix(x,loc,scale,shape,weight,Ncomp):
    
    tot = np.zeros(len(x))
    
    for i in range(Ncomp):
        tot += Bcdf(x,loc[i],scale[i],shape[i])*weight[i]
        
    return tot

def make_model(nB,x,y,sigma_beta):

    mixture_model = pm.Model()

    with mixture_model:
    
        # Define common priors
        sigma = pm.HalfCauchy('sigma', beta=sigma_beta)
        #weight = [pm.Uniform('weight_%d' %i, lower=0.05, upper=1.5) for i in range(nB)]
        weight = [pm.TruncatedNormal('weight_%d' %i,mu=1.0/nB, sigma=0.25,lower=0.05,upper=1.5) for i in range(nB)]
        #loc = [pm.Uniform('loc_%d' %i,lower=0, upper=np.max(x)*1.0) for i in range(nB)]
        loc = [pm.TruncatedNormal('loc_%d' %i,mu=1.7, sigma=0.6,lower=0.8,upper=1.5*np.max(x)) for i in range(nB)]
        #scale = [pm.HalfCauchy('scale_%d' %i,beta=0.15,testval=0.5) for i in range(nB)]
        #scale = [pm.Uniform('scale_%d' %i,lower=0.08,upper=0.6) for i in range(nB)]
        scale = [pm.TruncatedNormal('scale_%d' %i,mu=0.3, sigma=0.15,lower=0.08,upper=0.6) for i in range(nB)]
        shape = [pm.TruncatedNormal('shape_%d' %i,mu=-0.35,sigma=0.25,lower=-1.0,upper=0.1) for i in range(nB)]
        
        likelihood = pm.Normal('y', mu = mix(x,loc,scale,shape,weight,nB), sigma=sigma, observed=y)    
    
    return mixture_model

def unmix(nB,x,y,nsample=2000,tune=5000,sigma_beta=0.01,nchains=1,ncores=1):

    model = make_model(nB,x,y,sigma_beta)    
    trace = pm.sample(nsample, tune=tune, chains=nchains,cores=ncores,model=model,compute_convergence_checks=False)
        
    return trace, model

################# Plotting routines #################

def npCDFmix(x,loc,scale,shape,weight,Ncomp):
    
    tot = np.zeros(len(x))
    
    for i in range(Ncomp):
        tot += npBcdf(x,loc[i],scale[i],shape[i])*weight[i]
        
    return tot

def npPDFmix(x,loc,scale,shape,weight,Ncomp):
    
    tot = np.zeros(len(x))
    
    for i in range(Ncomp):
        tot += npBpdf(x,loc[i],scale[i],shape[i])*weight[i]
        
    return tot

def evaluate_CI(x,trace,nB,plow=2.5,phigh=97.5):

    nS = len(trace['model_logp']) #number of samples
    nX = len(x)
    
    cdf_lower = np.zeros((nX,nB+1))
    cdf_median = np.zeros((nX,nB+1))
    cdf_upper = np.zeros((nX,nB+1))
    pdf_lower = np.zeros((nX,nB+1))
    pdf_median = np.zeros((nX,nB+1))
    pdf_upper = np.zeros((nX,nB+1))
  
    for k in range(nX):
        cdf_tot = np.zeros(nS)
        pdf_tot = np.zeros(nS)
        
        for i in range(nB):
            Bweight = (trace.get_values('weight_%d' %i))
            Bloc = (trace.get_values('loc_%d' %i))
            Bscale = (trace.get_values('scale_%d' %i))
            Bshape = (trace.get_values('shape_%d' %i))
            
            cdf_B = np.zeros(nS)
            pdf_B = np.zeros(nS)
                       
            for j in range(nS):
                cdf_B[j] = npBcdf(x[k],Bloc[j],Bscale[j],Bshape[j])*Bweight[j]
                cdf_tot[j] += cdf_B[j]
                pdf_B[j] = npBpdf(x[k],Bloc[j],Bscale[j],Bshape[j])*Bweight[j]
                pdf_tot[j] += pdf_B[j]
            
            cdf_lower[k,i] = np.percentile(cdf_B,plow)
            cdf_median[k,i] = np.percentile(cdf_B,50)
            cdf_upper[k,i] = np.percentile(cdf_B,phigh)
            pdf_lower[k,i] = np.percentile(pdf_B,plow)
            pdf_median[k,i] = np.percentile(pdf_B,50)
            pdf_upper[k,i] = np.percentile(pdf_B,phigh)
            
        cdf_lower[k,-1] = np.percentile(cdf_tot,plow)
        cdf_median[k,-1] = np.percentile(cdf_tot,50)
        cdf_upper[k,-1] = np.percentile(cdf_tot,phigh)
        pdf_lower[k,-1] = np.percentile(pdf_tot,plow)
        pdf_median[k,-1] = np.percentile(pdf_tot,50)
        pdf_upper[k,-1] = np.percentile(pdf_tot,phigh)
    
        
    return cdf_lower, cdf_median, cdf_upper, pdf_lower, pdf_median, pdf_upper


def plot_mixture(x,y,trace,nB,ounits,Mn=1):
    
    xi = np.linspace(np.min(x),np.max(x),101)
    cdf_lower, cdf_median, cdf_upper, pdf_lower, pdf_median, pdf_upper = evaluate_CI(xi,trace,nB)

    colors = ((68/255,119/255,170/255),(238/255,102/255,119/255),(34/255,136/255,51/255),(102/255,204/255,238/255),
              (204/255,187/255,68/255),(187/255,187/255,187/255))
    
    plt.figure(figsize=(7.5,14))
    
    ########## PLOT IRM
    ax0 = plt.subplot2grid((5, 1), (0, 0), colspan=1,rowspan=2)
    for i in range(nB):
        ax0.fill_between(10**xi,cdf_lower[:,i]*Mn,cdf_upper[:,i]*Mn,alpha=0.5,color=colors[i],edgecolor=colors[i])
        ax0.plot(10**xi,cdf_median[:,i]*Mn,color=colors[i],label=('Comp %d' %(i+1)))
    
    ax0.fill_between(10**xi,cdf_lower[:,-1]*Mn,cdf_upper[:,-1]*Mn,alpha=0.5,color=colors[-1],edgecolor=colors[-1])
    ax0.plot(10**xi,cdf_median[:,-1]*Mn,color=colors[-1],label='Mixture')
        
    ax0.plot(10**x,y*Mn,'ok')
    ax0.set_xscale('log')
    ax0.minorticks_on()
    ax0.set_xlim([np.min(10**x),np.max(10**x)])
    ax0.set_ylim([0,np.max(y*Mn)*1.05])
    ax0.legend(fontsize=14)
    ax0.tick_params(labelsize=14)
    ax0.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax0.tick_params(labelsize=14)
    
    if ounits==0:
        ax0.set_ylabel('Magnetic moment [emu]',fontsize=14)
    elif ounits==1:
        ax0.set_ylabel('Magnetization [emu / cm'+u'\u00B3'+']',fontsize=14)
    elif ounits==2:
        ax0.set_ylabel('Mass magnetization [emu / g]',fontsize=14)
    elif ounits==3:
        ax0.set_ylabel('Magnetic moment [Am'+u'\u00B2'+']',fontsize=14)  
    elif ounits==4:
        ax0.set_ylabel('Magnetization [A / m]',fontsize=14)
    else:
        ax0.set_ylabel('Mass magnetization [Am'+u'\u00B2'+' / kg]',fontsize=14)

    
    #### PLOT CDF Residuals
    ax1 = plt.subplot2grid((5, 1), (2, 0), colspan=1,rowspan=1)
    hat = np.interp(x,xi,cdf_median[:,-1]*Mn)
    #
    hat_lower = np.interp(x,xi,cdf_lower[:,-1]*Mn)
    hat_upper = np.interp(x,xi,cdf_upper[:,-1]*Mn)
    ax1.fill_between(10**x,y*Mn-hat_lower,y*Mn-hat_upper,alpha=0.5,color='grey',edgecolor='grey')
    ax1.plot(10**x,y*Mn-hat,'k')
    ax1.set_xscale('log')
    ax1.minorticks_on()
    ax1.set_xlim([np.min(10**x),np.max(10**x)])
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax1.tick_params(labelsize=14)
    
    if ounits==0:
        ax1.set_ylabel('Residual [emu]',fontsize=14)
    elif ounits==1:
        ax1.set_ylabel('Residual [emu / cm'+u'\u00B3'+']',fontsize=14)
    elif ounits==2:
        ax1.set_ylabel('Residual [emu / g]',fontsize=14)
    elif ounits==3:
        ax1.set_ylabel('Residual [Am'+u'\u00B2'+']',fontsize=14)  
    elif ounits==4:
        ax1.set_ylabel('Residual [A / m]',fontsize=14)
    else:
        ax1.set_ylabel('Residual [Am'+u'\u00B2'+' / kg]',fontsize=14)
    #
    ax1.plot([1,10**np.max(x)],[0,0],'--k')
    
    
    ### PLOT DERIVATIVE #####
    ax2 = plt.subplot2grid((5, 1), (3, 0), colspan=1,rowspan=2)
    for i in range(nB):
        ax2.fill_between(10**xi,pdf_lower[:,i]*Mn,pdf_upper[:,i]*Mn,alpha=0.5,color=colors[i],edgecolor=colors[i])
        ax2.plot(10**xi,pdf_median[:,i]*Mn,color=colors[i],label=('Comp %d' %(i+1)))
    
    ax2.fill_between(10**xi,pdf_lower[:,-1]*Mn,pdf_upper[:,-1]*Mn,alpha=0.5,color=colors[-1],edgecolor=colors[-1])
    ax2.plot(10**xi,pdf_median[:,-1]*Mn,color=colors[-1],label='Mixture')
        
    Dy = np.diff(y*Mn)/np.diff(x)
    Dx = x[:-1]+np.diff(x)/2
    ax2.plot(10**Dx,Dy,'ok')
    ax2.set_xscale('log')
    ax2.minorticks_on()
    ax2.set_xlim([np.min(10**x),np.max(10**x)])
    ax2.legend(fontsize=14)
    ax2.set_xlabel('B [mT]',fontsize=14)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax2.tick_params(labelsize=14)

    if ounits==0:
        ax2.set_ylabel('Derivative [emu / $\log_{10}$(B[mT])]',fontsize=14)
    elif ounits==1:
        ax2.set_ylabel('Derivative [emu / cm'+u'\u00B3'+' / $\log_{10}$(B[mT])]',fontsize=14)
    elif ounits==2:
        ax2.set_ylabel('Derivative [emu / g / $\log_{10}$(B[mT])]',fontsize=14)
    elif ounits==3:
        ax2.set_ylabel('Derivative [Am'+u'\u00B2'+' / $\log_{10}$(B[mT])]',fontsize=14)  
    elif ounits==4:
        ax2.set_ylabel('Derivative [A / m / $\log_{10}$(B[mT])]',fontsize=14)
    else:
        ax2.set_ylabel('Derivative [Am'+u'\u00B2'+' / kg / $\log_{10}$(B[mT])]',fontsize=14)
