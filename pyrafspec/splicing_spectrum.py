import collections
import os
import warnings
from glob import glob
import joblib
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy import coordinates as coord
from astropy.io import fits
from laspec import normalization as lanorm
from .fitslist import *
import lightkurve as lk
#from log216tolist import *
from .bfosclog import *

'''
# binning law
$\lambda = 10^x$

$\lambda^\prime = ln(10) 10^x$

$\frac{\Delta \lambda}{\lambda} = ln(10)\Delta x$  & $R =  \frac{\lambda}{\Delta \lambda} $--> $\Delta x = \frac{1}{Rln(10)*n}$
'''

def combine_spectrum_sum(waves, fluxs, flux_errs, wave_dens=None, speclist=None, show=False, log10=True,funcnorm=None):
    '''combine spectra of E9+G10 by summing the overlap zone
    page 28 of specnan2014
    parameters
    -----------
    waves [2D array] e.g. waves.shape = (11, 2048)
    fluxs [2D array] same as the waves.shape
    flux_errs [2D array] same as the waves.shape
    speclist [1D array of list] the indices of spectra; e.g. [0, 1, 2, 3]
    weight [bool] sum with weight
    funcnorm [function] a normalization function, if is None, lanorm.lanorm.normalize_spectrum_spline
    return
    -----------
    wave_dens [1D array] the wavelength of the combine spectrum
    flux_norm [1D array] the normalized flux of combined spectrum
    flux_err [1D array] the flux error of combined spectrum
    '''
    if speclist is None: speclist = np.arange(waves.shape[0])
    if wave_dens is None:
       delta_wv = 0.1#np.min(np.diff(waves, axis=1))
       wv_min, wv_max= np.min(waves), np.max(waves)
       wave_dens = np.arange(wv_min, wv_max, delta_wv)
    

    #delta_wv = np.median(np.diff(np.log10(wave_dens)))
    delta_wv = np.median(np.diff(wave_dens))
    
    speclist = np.arange(waves.shape[0])
    N_wv = len(wave_dens)
    fluxs_dens = np.zeros((fluxs.shape[0], N_wv))
    fluxnorms_dens = np.zeros((fluxs.shape[0], N_wv))
    fluxerr2s_dens = np.zeros((fluxs.shape[0], N_wv))
    smoothed_dens = np.zeros((fluxs.shape[0], N_wv))
    #waves_dens = np.zeros((fluxs.shape[0], N_wv))
    
    for fluxi in speclist:
        #spec0 = rmcosimic(waves[fluxi], fluxs[fluxi], flux_errs[fluxi], sigma =5, itera=2, percentile_up=94, percentile_low=0.01,window_length=7, polyorder=2, 
        #                  show=False)
        #### norlize spectrum
        #flux_norm, flux_smoothed2 =lanorm.normalize_spectrum_spline(spec0.time, spec0.flux, p=1E-6, q=0.5, lu=(-1, 3), binwidth=30,
        #                              niter=5)
        #flux_norm, flux_smoothed2 =lanorm.normalize_spectrum_spline(waves[fluxi], fluxs[fluxi], p=1E-6, q=0.5, lu=(-1, 3), binwidth=30,
        #                              niter=5)
        ind = fluxs[fluxi] > 0
        _wave, _flux, _fluxerr = waves[fluxi][ind], fluxs[fluxi][ind],  flux_errs[fluxi][ind]
        if log10:
           wavenorm =10**(_wave)
        else:
           wavenorm = _wave
        if funcnorm is None:
           flux_norm, flux_smoothed2 =lanorm.normalize_spectrum_spline(wavenorm, _flux, p=1E-6, q=0.6, lu=(-1, 3), binwidth=40,niter=5)
        else:
           flux_norm, flux_smoothed2 =funcnorm(_wave, _flux, **kwargs)
        fluxs_dens[fluxi] = np.interp(wave_dens, _wave[50:-50], _flux[50:-50], right=0., left=0.)
        smoothed_dens[fluxi] = np.interp(wave_dens,  _wave[50:-50], flux_smoothed2[50:-50], right=0., left=0.)
        #fluxnorms_dens[fluxi] = np.interp(wave_dens, spec0.time[50:-50], flux_norm[50:-50], right=0., left=0.)
        err_frac = np.append(1, np.diff(_wave)/delta_wv)
        #print(delta_wv)
        
        err2_tmp =  _fluxerr**2*err_frac
        #print(np.median(err_frac))
        fluxerr2s_dens[fluxi] = np.interp(wave_dens, _wave[50:-50], err2_tmp[50:-50], right=0., left=0.)
        #waves_dens[fluxi] = wave_dens
   
    fracs = 1./np.sum(fluxs_dens > 0, axis=0)
    fracs[np.isinf(fracs)] = 0.
    _fluxs_sum = np.sum(fluxs_dens*fracs, axis=0)
    _smooth_sum = np.sum(smoothed_dens*fracs, axis=0)
    _fluxerr2s_sum = np.sum(fluxerr2s_dens*fracs**2, axis=0)
    _fluxnorm_sum = _fluxs_sum/_smooth_sum
    _fluxnorm_err = np.sqrt(_fluxerr2s_sum)/_smooth_sum
    if show:
       spec = lk.LightCurve(time=wave_dens,  flux=_fluxnorm_sum, flux_err=None)
       spec.plot()
       plt.xlabel('wavelength')
       plt.ylim(0, 1.3)
    return wave_dens, _fluxnorm_sum, _fluxnorm_err

def binning_spectrum(waves, fluxs, flux_errs, wave_bin=np.arange(3574, 8940, 0.8)):
    '''binning the spectrum
    
    parameters:
    -------------------
    waves [2D array]
    fluxs [2D array]
    flux_errs [2D array]
    returns:
    --------------------
    wave_new [2D array]
    flux_bin [2D array]
    fluxerr_bin [2D array]
    n [2D array]
    '''
    norder = waves.shape[0]
    wave = (wave_bin[1:] + wave_bin[:-1])/2.
    nwvl = len(wave)
    wave_new = np.zeros((norder, nwvl))
    flux_bin = np.zeros((norder, nwvl))
    fluxerr_bin = np.zeros((norder, nwvl))
    n = np.zeros((norder, nwvl))
    for _i, _wave in enumerate(waves):
        _flux = fluxs[_i]
        _fluxerr2 = flux_errs[_i]**2
        _fluxerr2_bin,_ = np.histogram(_wave, bins=wave_bin, weights=_fluxerr2)
        _flux_bin,_ = np.histogram(_wave, bins=wave_bin, weights=_flux)
        _n,_ = np.histogram(_wave, bins=wave_bin)
        flux_bin[_i] = _flux_bin/_n
        fluxerr_bin[_i] =  np.sqrt(_fluxerr2_bin)/_n
        wave_new[_i] = wave
        n[_i] = _n
    _ind = np.isnan(flux_bin)
    flux_bin[_ind] = 0.
    fluxerr_bin[_ind] = 0.
    return  wave_new, flux_bin, fluxerr_bin, n


def combine_wave_flux(fwave, fflux):
    '''combine the wavelength and fflux to a dic
    parameters:
    ------------------
    fwave [str] file name of wavelength e.g. /media/backup/216BFOSC/20201218_bfosc/fear-20201218-0056.fit.dump
    fflux [str] file name of flux, e.g. '/media/backup/216BFOSC/20201218_bfosc/star-20201218-0055.fit.dump'
    '''
    star = {}
    star['wavelength'] = joblib.load(fwave)
    star['flux'] = joblib.load(fflux)
    return star

def get_loglam(R, lam_start, lam_end, N=3):
    ''' get log10(lambda) with a proper sample interval 
    parameters:
    ---------------
    R: [float] spectral resolution (BFOSC: R ~ 1600)
    lam_start: [float]: the start of wavelength
    lam_end: [float]: the end of wavelength
    N: [int] oversampling, (typical oversampling: N=3 --> 3pixel=FWHM)
       R = lambda/FWHM --> (log10(lambda))' =1/(lambda*ln(10)) --> FWHM = 1/(R*ln(10))
    returns:
    log10lam: [array] log10(lambda)
    '''
    deltax = 1/R/N/np.log(10)
    log10lam = np.arange(np.log10(lam_start), np.log10(lam_end), deltax)
    return log10lam

def rmcosimics(waves, fluxs, fluxerrs,**kwargs):
    for _i, wave in enumerate(waves):
        flux = fluxs[_i].copy()
        fluxerr = fluxerrs[_i].copy()
        lc = rmcosimic(wave, flux, fluxerr, **kwargs)
        waves[_i] = lc.time
        fluxs[_i] = lc.flux
        fluxerrs[_i] = lc.flux_err
    return waves, fluxs, fluxerrs


def rmcosimic(wave, flux, fluxerr, sigma =8, itera=2, percentile_up=100, percentile_low=0.01, window_length=7, 
              polyorder=2, 
              show=False):
    '''remove cosimic line of spec
    parameters:
    --------------
    wave [1d array]
    flux [1d array]
    fluxerr [1d array]
    sigma [float] exclude the flux which sigma*std
    itera [int] iteration number to calculate std
    percentile_up, percentile_low [float] calculate std by using the points that flux-trend_flux in the interval of [percentile_low, percentile_up]
    window_length [odd]
    '''
    spec = lk.LightCurve(time=wave,  flux=flux, flux_err=fluxerr)
    lc = spec.copy()
    ind_good = np.ones(len(lc.flux), dtype=np.bool)
    lc_detrend, lc_trend = lc.flatten(window_length=window_length, polyorder=polyorder, return_trend=True,
                break_tolerance=5, niters=3, sigma=3, mask=None)
    lcflux = spec.flux[ind_good] 
    lctrendflux = lc_trend.flux[ind_good]
    for i in np.arange(itera):
        diff_flux = lcflux  - lctrendflux
        std = np.std(diff_flux[(diff_flux < np.percentile(diff_flux, percentile_up)) & (diff_flux > np.percentile(diff_flux, percentile_low))])
        ind_good = diff_flux < sigma*std
        lcflux = lcflux[ind_good]
        lctrendflux= lctrendflux[ind_good]
        
    ind_bad = ((lc.flux - lc_trend.flux) > sigma*std) | (lc.flux < 0)
    #print(std)
    lc.flux[ind_bad] = np.interp(lc.time[ind_bad], lc.time[~ind_bad], lc.flux[~ind_bad])
    lc.flux_err[ind_bad] = np.sqrt(np.interp(lc.time[ind_bad], lc.time[~ind_bad], (lc.flux_err[~ind_bad])**2))
    if show:
       spec.scatter(color='b', )
       plt.plot(lc.time, lc.flux)
       plt.plot(lc_trend.time, lc_trend.flux, lw=0.5, c='r')
    return lc

def splicing_irafspectrum(filename, R=1500, N=3, lam_start=3600, lam_end=8900, pix=[0, -1], funcnorm= None, orders=None, **kwargs):
    ''' splicing the spectral segmentation of E9G10 to an spectrum
    parameters:
    --------------
    filename: [str] the fits file name (1ds spectra extracted by iraf)
    R: [float] spectral resolution
    N: [int] oversampling
    pix: [list] len(pix) =2 and the elements should be int; the [start, end] pixel of each order of spectra
    orders: [list]
    returns:
    logwave: [array] log10(wavelength)
    fluxnorm: [array] normalized flux
    fluxnorm_err: [array] the error of normalized flux
    '''
    spec = load_multispec(filename)
    if orders is None: orders = spec.orders
    N_orders = len(orders)
    waves = np.zeros((N_orders, pix[1]- pix[0]))
    fluxs = np.zeros_like(waves)
    flux_errs = np.zeros_like(waves)
    
    for _i, order in enumerate(orders):
        _wave = spec[order].wv[pix[0]:pix[1]]
        _flux = spec[order].flux[pix[0]:pix[1]]
        _indsort = np.argsort(_wave)
        waves[_i] = _wave[_indsort]
        fluxs[_i] = spec[order].flux[pix[0]:pix[1]][_indsort] 
    waves, fluxs, flux_errs = rmcosimics(waves, fluxs, flux_errs, **kwargs)
    logwaves = np.log10(waves)
    log10lam = get_loglam(R, lam_start, lam_end, N=3)
    waves, fluxs, flux_errs, _n = binning_spectrum(logwaves, fluxs, flux_errs, wave_bin=log10lam)
    logwave, fluxnorm, fluxnorm_err = combine_spectrum_sum(waves, fluxs, flux_errs, wave_dens=log10lam, speclist=None, funcnorm=funcnorm)
    return logwave, fluxnorm, fluxnorm_err

def splicing_spectrum(fwave, fflux, R=1500, N=3, lam_start=3600, lam_end=8900, funcnorm= None):
    ''' splicing the spectral segmentation of E9G10 to an spectrum
    parameters:
    --------------
    fflux: [str] a dump file of star produced by bfosc_pipeline.py (e.g. star-***.fit.dump)
    fwave: [str] a dump file of lamp produced by bfosc_pipeline.py (e.g. fear-***.fit.dump)
    R: [float] spectral resolution
    N: [int] oversampling
    returns:
    logwave: [array] log10(wavelength)
    fluxnorm: [array] normalized flux
    fluxnorm_err: [array] the error of normalized flux
    '''
    star = combine_wave_flux(fwave, fflux)
    waves = star['wavelength']['wave_solu']
    fluxs = star['flux']['spec_extr']
    flux_errs = star['flux']['err_extr']
    waves, fluxs, flux_errs = rmcosimics(waves, fluxs, flux_errs, **kwargs)
    logwaves = np.log10(waves)
    log10lam = get_loglam(R, lam_start, lam_end, N=3)
    waves, fluxs, flux_errs, _n = binning_spectrum(logwaves, fluxs, flux_errs, wave_bin=log10lam)
    logwave, fluxnorm, fluxnorm_err = combine_spectrum_sum(waves, fluxs, flux_errs, wave_dens=log10lam, speclist=None, funcnorm=funcnorm)
    return logwave, fluxnorm, fluxnorm_err

def write2fits(fstar, flamp, dire, fout=None):
    ''' write the header of fits file of star raw data  and jointed spectrum to fits file
    parameters:
    -------------
    fstar: [str] file name of star spectrum (e.g. 202110220019_SPECSTARGET_BD+25d4655_slit16s_G10_E9)
    flamp: [flamp] file name of lamp (e.g. '202110220021_SPECSLAMP_FeAr_slit16s_G10_E9')
    dire: [str] directory stored raw data and dump file
    fout: [str] the name of output fits file, if None: fout = f'{dire}/{fstar}.fits'
    -------------
    '''
    fwave, fflux, fstarall, flampall = dumpname(fstar, flamp, dire)
    logwave, flux, flux_err = splicing_spectrum(fwave, fflux, R=1500, N=3, lam_start=3600, lam_end=8900)
    table = fits.BinTableHDU.from_columns([
                    fits.Column(name='loglambda',  format='E', array=np.array(logwave, dtype=np.float32)),
                    fits.Column(name='flux',format='E', array=np.array(flux, dtype=np.float32)),
                    fits.Column(name='error',format='E', array=np.array(flux_err, dtype=np.float32))
                                         ])
    starheader = fits.getheader(fstarall)
    data0 = np.zeros(1,dtype=np.float32)
    hdu = fits.PrimaryHDU(data0)
    hdu.header = starheader
    hdu.header['fstar'] = (fstar, 'file name of raw image of star')
    hdu.header['flamp'] = (flamp, 'file name of raw image of lamp')
    hdul = fits.HDUList([hdu, table])
    if fout is None:
       fout = os.path.join(dire, f'{fstar}.fits')
    hdul.writeto(fout, overwrite=True)


def pyrafspc1d2fits(filename, loglambda, flux, fluxerr, dire=None, fout=None):
    ''' write the header of fits file of star raw data  and jointed spectrum to fits file
    parameters:
    -------------
    filename: [str] file name of 1d spectrum (e.g. 202110220019_SPECSTARGET_BD+25d4655_slit16s_G10_E9)
    dire: [str] directory stored raw data and dump file
    fout: [str] the name of output fits file, if None: fout = f'{dire}/{fstar}.fits'
    -------------
    '''
    if dire is None:
       dire = os.path.dirname(filename)

    hdu = fits.open(filename)[0]
    header= hdu.header
    dateobs = header['DATE-OBS'].split('T')[0]
    utstat = dateobs+'T' +header['UTSTART']
    utstop = dateobs +'T' +header['UTSTOP']
    times = Time([utstat, utstop],scale='utc', format='isot' )
    jd = np.mean(times.jd)
    header['JD'] = jd
    ##########################################################################################
    #--------------light_travel_time----------------------------------------------------------
    radec = header['RA']+' '+header['DEC']
    ip_peg = coord.SkyCoord(radec, unit=(u.hourangle, u.deg), frame='icrs')
    site = coord.EarthLocation.from_geodetic(lat=26.6951*u.deg, lon=100.03*u.deg, height=3200*u.m)
    times_ltt = Time(jd, format='jd', scale='utc', location=site)
    ltt = times_ltt.light_travel_time(ip_peg,'barycentric')
    barycorr = ip_peg.radial_velocity_correction(obstime=Time(times_ltt.iso), location=site)
    #----------------------------------------------------------------------------------------
    header['BJD'] = jd+ltt.jd
    header['barycorr'] = (barycorr.to('km/s').value, 'km/s')
    ###########################################################################################
    
    table = fits.BinTableHDU.from_columns([
                        fits.Column(name='loglambda',  format='E', array=np.array(loglambda, dtype=np.float32)),
                        fits.Column(name='flux',format='E', array=np.array(flux, dtype=np.float32)),
                        fits.Column(name='error',format='E', array=np.array(fluxerr, dtype=np.float32))
                                         ])  
    hdul = fits.HDUList([hdu, table])

    if fout is None:
       basename = os.path.basename(filename)[:-5]+'_splicing.fits'
       fout = os.path.join(dire, basename)
    hdul.writeto(fout, overwrite=True)


def dumpname(fstar, flamp, dire):
    ''' get dump file from lamp and star file
    parameters:
    -------------
    fstar: [str] file name of star spectrum (e.g. 202110220019_SPECSTARGET_BD+25d4655_slit16s_G10_E9)
    flamp: [flamp] file name of lamp (e.g. '202110220021_SPECSLAMP_FeAr_slit16s_G10_E9')
    dire: [str] directory stored raw data and dump file
    returns:
    -------------
    fwave: [str] dump file name of lamp
    fflux: [str] dump file name of star
    fstarall: [str]
    flampall: [str]
    '''
    stardump = f'star-{fstar}.fit.dump'
    lampdump = f'star-{flamp}.fit.dump'
    fwave = os.path.join(dire, lampdump)
    fflux = os.path.join(dire, stardump)
    fstarall = os.path.join(dire, f'{fstar}.fit')
    flampall = os.path.join(dire, f'{fstar}.fit')
    return fwave, fflux, fstarall, flampall


if __name__ == '__main__':
   dire = '/home/lcq/media/backup/216BFOSC/20211022_bfosc'
   logfile = 'liuchao_bfosc.log' 
   fname = os.path.join(dire, logfile)
   star_list, lamp_list = match_star2lamp(fname, fout=None)
   for fstar, flamp in zip(star_list, lamp_list):
       write2fits(fstar, flamp, dire, fout=None)
