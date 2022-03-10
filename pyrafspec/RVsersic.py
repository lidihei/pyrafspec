import matplotlib.pylab as plt
# import PyAstronomy
# from PyAstronomy.pyTiming import pyPDM  #need to install PyAstronomy

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d, interp2d
import matplotlib
from matplotlib import cm
from matplotlib.colors import Normalize
matplotlib.rc('xtick',labelsize=12)
matplotlib.rc('ytick',labelsize=12)
from astropy.table import Table
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'inline')
from glob import glob
from pathlib import Path

from astropy import time
import astropy.io.fits as pyfits, numpy as np
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.utils import iers
import astropy.units as u
import os, sys
from astropy.io import fits
from astropy.time import Time
from scipy.optimize import curve_fit


def gaussian(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))/(np.sqrt(2*np.pi)*sigma)

def func(x,mu, a1, sig1,a2, sig2,a3, sig3,a4, sig4):
    f1 = gaussian(x,a1,mu,sig1)
    f2 = gaussian(x,a2,mu,sig2)
    f3 = gaussian(x,a3,mu,sig3)
    f4 = gaussian(x,a4,mu,sig4)
    flux = f1+ f2+f3+f4
    return flux

def sersic(x,a,b,c,x0):
    return 1.-a*np.exp(-(np.abs(x-x0)/np.abs(b))**c)


def fit_sersic(ww,ff,p0=[0.8,4.,1.,4101.0]):
    popt,pcov = curve_fit(sersic,ww,ff,p0)
    ffm = sersic(ww,popt[0],popt[1],popt[2],popt[3])
    #plt.plot(ww,ff,'k')
    #plt.plot(ww,ffm,'b--')
    dy = (ffm-ff)
    ind = dy<np.percentile(dy,95)
    popt,pcov = curve_fit(sersic,ww[ind],ff[ind],p0)
    ffm = sersic(ww,popt[0],popt[1],popt[2],popt[3])
    #plt.plot(ww,ffm,'g--')
    dy = (ffm-ff)
    ind = dy<np.percentile(dy,95)
    popt,pcov = curve_fit(sersic,ww[ind],ff[ind],p0)
    ffm = sersic(ww,popt[0],popt[1],popt[2],popt[3])
    #plt.plot(ww,ffm,'r--')
    dy = (ffm-ff)
    ind = dy<np.percentile(dy,95)
    popt,pcov = curve_fit(sersic,ww[ind],ff[ind],p0)
    ffm = sersic(ww,popt[0],popt[1],popt[2],popt[3])
    #plt.plot(ww,ffm,'m--')
    #plt.show()
    return popt,pcov
    
def lnprob_sersic(x,w,y):
    ym = sersic(w,x[0],np.abs(x[1]),x[2],x[3])
    dy = np.abs(ym-y)
    ind = dy<np.percentile(dy,99)
    return -np.sum(dy[ind]**2)

def sersic_mcmcfit(w,f,y0):
    # MCMC sampling
    n = len(f)
    
    #start to configure emcee
    nwalkers = 20
    ndim = 4
    p0=np.zeros((nwalkers,ndim))
    p0[:,0] = np.random.rand(nwalkers)*1
    p0[:,1] = np.random.rand(nwalkers)*10
    p0[:,2] = np.random.rand(nwalkers)*1
    p0[:,3] = np.random.rand(nwalkers)*1+y0
    #p0[:,4] = np.random.rand(nwalkers)*0.1+1.
      
    sampler = emcee.EnsembleSampler(nwalkers,             ndim, lnprob_sersic, args=[w,f])
    
    pos, prob, state = sampler.run_mcmc(p0, 100)
    sampler.reset()
    
    sampler.run_mcmc(pos, 10000)
    
    samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
    corner.corner(samples)
    popt = np.median(samples, axis=0)
    pcov = np.zeros((ndim,ndim))
    for i in range(ndim):
        for j in range(ndim):
            pcov[i,j] = (np.sum((samples[:,i]-popt[i])*                (samples[:,j]-popt[j])))/len(samples)
    return popt, pcov


def dlambda2rv(lam, lamerr, lam0):
    '''
    parameters:
    -------------
    lam: [float] the wavelength of emission or absorption lines (in angstrom)
    lamerr: [float] the error of the line
    lam0: the rest reference wavelength of the line
    returns:
    ----------
    rv: [float] velocity (km/s)
    rverr: [float] velocity error (km/s)
    '''
    c = 299792.458
    rv = (lam - lam0)/lam0 * c
    rverr = lamerr*c/lam0
    return rv, rverr

def rv_bysersic(wave, flux, lam0=4861, waverange=[4840, 4890], p0=[0.6,4.5,0.8,4862.], show=False, v =0):
    '''Halpha: lam0 = 6562.79, waverange = [6550, 6580], p0 = [0.8,4.,1.,6562.]
       Hbeta:  lam0 = 4861.35, waverange = [4840, 4890], p0 = [0.6,4.5,0.8,4862.]
       Hgamma: lam0 = 4340.46, waverange = [4310, 4360], p0 = [0.8,4.,1.,4340.]
    returns:
    ------------
    rv
    rverr
    '''
    c = 299792.458
    wave = wave + v/c*wave
    ws, we = waverange
    indw = (wave > ws) & (wave < we) & (flux  > 0) & (flux < 2.)
    try:
        popt,pcov = fit_sersic(wave[indw],flux[indw],p0=p0)
        lam, lamerr = popt[-1], np.sqrt(pcov[-1,-1])
        rv, rverr = dlambda2rv(lam, lamerr, lam0)
        if show:
           xx = np.linspace(ws, we, 200)
           ee = np.ones_like(xx)
           vv, _ = dlambda2rv(xx, ee, lam0) 
           yfit = sersic(xx, *popt)
           fig, ax1 = plt.subplots(1,1,figsize=[7,3])
           plt.plot(wave[indw], flux[indw], 'k', lw=0.8)
           plt.plot(xx, yfit, 'r', lw=1.2)
           plt.axvline(x = lam0, color='r', lw=1, ls='--')
           plt.xlabel(r'wavelength $\AA$')
           ax2 = ax1.twiny()
           plt.plot(vv, yfit, 'r', lw=1.2)
           plt.axvline(x=rv, color='b', lw=1)
           plt.xlabel('v (km/s)'+f' rv = {rv:.2f}' + r'$\pm' f'{rverr:.2f}''$')
    except:
        print('can not fit this line by sersic')
        rv, rverr = -99999, -99999
    return rv, rverr

def rv_bysersic2(wave, flux, fluxerr, lam0=4861, waverange=[4840, 4890], p0=[0.6,4.5,0.8,4862.], show=False, v =0):
    '''Halpha: lam0 = 6562.79, waverange = [6550, 6580], p0 = [0.8,4.,1.,6562.]
       Hbeta:  lam0 = 4861.35, waverange = [4840, 4890], p0 = [0.6,4.5,0.8,4862.]
       Hgamma: lam0 = 4340.46, waverange = [4310, 4360], p0 = [0.8,4.,1.,4340.]
    returns:
    ------------
    rv
    rverr
    chi2nv: [float] the reduced chi2
    n: [int] the number of points used to sersic fit
    '''
    c = 299792.458
    wave = wave + v/c*wave
    ws, we = waverange
    indw = (wave > ws) & (wave < we) & (flux  > 0) & (flux < 2.)
    n = np.sum(indw)
    try:
        popt,pcov = fit_sersic(wave[indw],flux[indw],p0=p0)
        yfit0 = sersic(wave[indw], *popt)
        lam, lamerr = popt[-1], np.sqrt(pcov[-1,-1])
        rv, rverr = dlambda2rv(lam, lamerr, lam0)
        chi2 = np.nansum((flux[indw] - yfit0)**2/fluxerr**2)
        chi2nv = chi2/(n-4)
        if show:
           xx = np.linspace(ws, we, 200)
           ee = np.ones_like(xx)
           vv, _ = dlambda2rv(xx, ee, lam0) 
           yfit = sersic(xx, *popt)
           fig, ax1 = plt.subplots(1,1,figsize=[7,3])
           plt.plot(wave[indw], flux[indw], 'k', lw=0.8)
           plt.plot(xx, yfit, 'r', lw=1.2)
           plt.axvline(x = lam0, color='r', lw=1, ls='--')
           plt.xlabel(r'wavelength $\AA$')
           ax2 = ax1.twiny()
           plt.plot(vv, yfit, 'r', lw=1.2)
           plt.axvline(x=rv, color='b', lw=1)
           plt.xlabel('v (km/s)'+f' rv = {rv:.2f}' + r'$\pm' f'{rverr:.2f}''$')
    except:
        print('can not fit this line by sersic')
        rv, rverr, chi2nv = -99999, -99999, -99999
    return rv, rverr, chi2nv, n
