from laspec.extern.interpolate import SmoothSpline
from laspec import convolution
from PyAstronomy import pyasl
import numpy as np
from .broadsed import vgconv, rotconv
from scipy.fftpack import fft
from scipy.fftpack import ifft


def reinterp_wave(ipar, Tpars, Tlusty_sps, wave0, wave, linearinterp=False):
    '''interpolate spectrum in a new wavelength
    parameters:
    ------------------
    ipar [int]
    Tpars [ 2D array] parameters array of TLUSTY grid e.g. np.array([[teff, logg, z, vt]])
    Tlusty_sps [2D array] spectra array of TLUSTY grid
    wave0 [1D array] the original wavelength of TLUSTY spectra
                     e.g. BTLUSTY optical: wave = np.arange(3200.01, 9998.7, 0.01); 
                          OTLUSTY: wave = np.arange(3000, 7500, 0.01)
    wave [1D array] the new wavelength
    returns:
    -------------
    pars [1D array] which has the same parameters with Tpars
    spec [1D array] has same length as wave
    '''
    #print(len(Tlusty_sps[ipar]), len(wave0))
    if linearinterp:
        spec = np.interp(wave, wave0, Tlusty_sps[ipar])
    else:
        func = SmoothSpline(wave0, Tlusty_sps[ipar], p=1)
        spec = func(wave)
    return Tpars[ipar], spec


def smooth_by_bin(time, flux, bins):
    '''smooth light curve by bins
    parameters:
    -------------------
    time: [array] could be phase
    flux: [flux]
    bins: [array]
    returns:
    -------------------
    time_bins: [array] time_bins = (bins[1:] + bins[:-1])/2
    flux_smooth: [array]
    '''
    flux_bins, _bins = np.histogram(time, weights=flux, bins=bins)
    n_bins, _bins = np.histogram(time, bins=bins)
    time_bins = (bins[1:] + bins[:-1])/2.
    flux_smooth = flux_bins/n_bins
    return time_bins, flux_smooth

def broad_tlusty(ipar, pix=0.01, mwv= 4861, R = 2000, pars=None, sps=None, wave=None, laspec_conv=False):
    '''broad tlusty spectrum with a resolution
    parameters:
    ------------------
    ipar [int]
    pix [float] the sample interval of spectra
    mwv [float] the median (or mean) wavelength of spectrum which is used to caculate the gaussian window by resoution
    R [float] spectrum resolution which we want to get
    Tpars [ 2D array] parameters array of TLUSTY grid e.g. np.array([[teff, logg, z, vt]])
    Tlusty_sps [2D array] spectra array of TLUSTY grid
    wave0 [1D array] the original wavelength of TLUSTY spectra
                     e.g. BTLUSTY optical: wave = np.arange(3200.01, 9998.7, 0.01); 
                          OTLUSTY: wave = np.arange(3000, 7500, 0.01)
    wave [1D array] the new wavelength
    returns:
    -------------
    pars [1D array] which has the same parameters with Tpars
    spec [1D array] has same length as wave
    '''
    flux = sps[ipar]
    if laspec_conv:
       if wave is None: wave  =  nnp.arange(3201, 7499, 0.01)
       wv_new, convsp = convolution.conv_spec(wave, flux, R_hi=300000., R_lo=R, over_sample_additional=1,
              gaussian_kernel_sigma_num=5., wave_new=wave,
              wave_new_oversample=1, verbose=False, return_type='array')
    else:
        fwhm = mwv/R
        sigma =  fwhm / (2.0 * np.sqrt(2. * np.log(2.)))
        nn2 = sigma/pix
        x = np.arange(-nn2*2,nn2*2+1,1)
        s2 = stats.norm.pdf(x,0,nn2)
        s = s2/np.sum(s2)
        convsp = np.convolve(flux,s,mode='same')
    return Tpars[ipar], convsp


def normalize_template(ipar, wvl, pars, sps):
    '''Apply rotational broadening to a spectrum. The formulae given in Gray's "The Observation
       and Analysis of Stellar Photospheres". 
       
    parameters:
    ------------------
    ipar [int]
    vsini : [float]
        Projected rotational velocity [km/s].
    epsilon : [float]
        Linear limb-darkening coefficient (0-1).
    wvl : array
        The wavelength array [A]. Note that a
        regularly spaced array is required.
    edgeHandling : string, {"firstlast", "None"}
        The method used to handle edge effects.
    Tpars [ 2D array] parameters array of TLUSTY grid e.g. np.array([[teff, logg, z, vt]])
    Tlusty_sps [2D array] spectra array of TLUSTY grid
    
    returns:
    -------------
    pars [1D array] which has the same parameters with Tpars
    spec [1D array] has same length as wave
    
    '''
    flux = sps[ipar]
    nflux = normalization.normalize_spectrum_spline(wvl, flux)
    return pars[ipar], nflux[0]

def rotbroad_template(ipar, epsilon, vsini, wvl, edgeHandling='firstlast', pars=None, sps=None):
    '''Apply rotational broadening to a spectrum. The formulae given in Gray's "The Observation
       and Analysis of Stellar Photospheres". 
       
    parameters:
    ------------------
    ipar [int]
    vsini : [float]
        Projected rotational velocity [km/s].
    epsilon : [float]
        Linear limb-darkening coefficient (0-1).
    wvl : array
        The wavelength array [A]. Note that a
        regularly spaced array is required.
    edgeHandling : string, {"firstlast", "None"}
        The method used to handle edge effects.
    Tpars [ 2D array] parameters array of TLUSTY grid e.g. np.array([[teff, logg, z, vt]])
    Tlusty_sps [2D array] spectra array of TLUSTY grid
    
    returns:
    -------------
    pars [1D array] which has the same parameters with Tpars
    rflux [1D array] rotational broadening spectra
    
    '''
    flux = sps[ipar]
    if vsini > 0:
       rflux = pyasl.rotBroad(wvl, flux, epsilon, vsini, edgeHandling=edgeHandling)
    else:
       rflux = flux
    return np.append(Tpars[ipar], vsini), rflux

def interpolate_fluxrv(times, phase, flux_phase, tc0, period=None):
    '''interpolate flux or rv by the fluxs at one period
    parameters:
    --------------
    times: [1D array] the times which one want to inteplate
    phase: [1D array] the phase which is given by phoebe
    flux_phase: [1D array] the flux of phase
    tc0: [float] Zeropoint date at superior conjunction (periastron passage) of the primary component 
    period: [period]
    phase0: if phase0 = 0.5 phase: -0.5, 0.5
    returns: 
    -------------
    fluxes: [1 D array]
    '''
    if period is not None:
       tminust0 =times - tc0
       tdividep =  tminust0/period
       _phase = np.mod(tminust0, period)/period
    else:
       _phase = times
    phase = np.append(phase, phase+1)
    flux_phase = np.append(flux_phase, flux_phase)
    fluxes = np.interp(_phase, phase, flux_phase)
    return fluxes


def rvcorr_spec(wave, flux, fluxerr, rv, wave_new=None, left=np.nan, right=np.nan, interp1d=None):
    ''' correct spectrum with radial velocity
    parameters:
    ------------
    wave [1d array]
    flux [1d array]
    fluxerr [1d array]
    barycorr [float] barycentric radial velocity in units of km/s
    
    returns:
    ----------
    flux_bc [1d array]
    fluxerr_bc [1d array]
    '''
    wvl = wave
    flux = flux
    ## Page 71 of An Introduction to Close Binary Stars
    c = 299792.458
    beta = rv/c
    lgwvl = np.log(wvl)
    gamma =(1+beta)/(1-beta)
    _lgwvl = lgwvl + 0.5*np.log(gamma)
    
    if wave_new is None:
       return _lgwvl, flux, fluxerr
    else: lgwvl = np.log(wave_new)
       if interp1d is None:
          flux_bc = np.interp(lgwvl, _lgwvl, flux, left=left, right=right)
          err2 = np.interp(lgwvl, _lgwvl, fluxerr**2, left=left, right=right)
       else:
          flux_bc = interp1d(_lgwvl, flux, kind='linear',fill_value='extrapolate')(lgwvl)
          err2 = interp1d(_lgwvl, fluxerr**2, kind='linear',fill_value='extrapolate')(lgwvl)
       fluxerr_bc = np.sqrt(err2)
    return flux_bc, fluxerr_bc

def lambda2vel(lambda0, wvl, flux, dwvl= 50, vel_dens=None):
    '''conver wavelength to velocity for a specific line (lambda0)
    parameters:
    lambda0 [float] the center wavelength of an absoption (or emission) line, which is in the units of angstrom
    wlv [1d array] the wavelength of a spectrum
    flux [1d array] the flux of a spectrum
    dwvl [float] the width of a absorption (emssion) line
    vel_dens [array] 
    returns:
    vel [1d array] the corresponding velocity of wavelength
    flux1 [1d array] the flux betweent [lambda0-dwvl, lambda0+dwvl]
    '''
    ind = (wvl < (lambda0 + dwvl)) & (wvl > (lambda0 - dwvl))
    flux1 = flux[ind]
    wvl1 = wvl[ind]
    c = 299792.46
    vel = (wvl1 - lambda0)/lambda0 * c
    if vel_dens is not None:
       flux_dens = np.interp(vel_dens, vel, flux1, right=1., left=1.)
       return vel_dens, flux_dens
    return vel, flux1

def composite_spectra(wave1, flux1, R1, wave2, flux2, R2,\
                     wave=None, rv1=1, vsini1=0, rv2=0, vsini2=0, epsilon1=0.2, epsilon2=0.2):
    '''composite spectra of two stars
    parameters:
    -------------------
    wave1: [array] the wavelength of star1
    flux1: [array] the flux of star1
    R1: [array] the radius of star1
    wave2: [array] the wavelength of star2
    flux2: [array] the flux of star2
    R2: [array] the radius of star2
    wave: [array] the new wavelength
    rv1: [float] the radial velocity of star1 km/s
    vsini1: [float] the rotational velocity of star1 km/s
    rv2: [float] the radial velocity of star2 km/s
    vsini2: [float] the rotational velocity of star2 km/s
    epsilon1: [float] the limb-darkening coefficient of star1
    epsilon2: [float] the limb-darkening coefficient of star2 
    returns:
    ------------
    flux: [array]  flux = flux1*R1**2 + flux2*R2**2
    '''
    # shift with rv
    flux2, _fluxerr2 = rvcorr_spec(wave2, flux2, np.zeros(len(wave2)),\
                                    rv2, wave_new=wave, left=np.nan, right=np.nan, interp1d=None)
    flux1, _fluxerr1 = rvcorr_spec(wave1, flux1, np.zeros(len(wave1)),\
                                    rv1, wave_new=wave, left=np.nan, right=np.nan, interp1d=None)
    #-----------------rotating broad
    if vsini1 != 0:
        wave1_vsini, flux1_vsini = rotconv(wave, flux1, epsilon1, vsini1)
    else:
        wave1_vsini, flux1_vsini = wave, flux1
    if vsini2 != 0:
        wave2_vsini, flux2_vsini = rotconv(wave, flux2, epsilon2, vsini2)
    else:
        wave2_vsini, flux2_vsini = wave, flux2
    #--------------------reinterp wave
    flux2, _fluxerr2 = rvcorr_spec(wave2_vsini, flux2_vsini, np.zeros(len(wave2_vsini)),\
                                0, wave_new=wave, left=np.nan, right=np.nan, interp1d=None)
    flux1, _fluxerr1 = rvcorr_spec(wave1_vsini, flux1_vsini, np.zeros(len(wave1_vsini)),\
                                    0, wave_new=wave, left=np.nan, right=np.nan, interp1d=None)
    flux = R1**2*flux1 + R2**2*flux2
    #if show:
    #   fig, ax = plt.subplots(1,1, figsize=[10,5])
    #   nflux = normalization.normalize_spectrum_spline(wave, flux, niter=5)
    #   #nflux1 = normalization.normalize_spectrum_spline(wave, flux1, niter=5) #(R1**2+R2**2)*flux1/nflux[1]
    #   #nflux2 = normalization.normalize_spectrum_spline(wave, flux2, niter=5) #(R1**2+R2**2)*flux2/nflux[1]
    #   flux1m = np.median(flux1[~np.isnan(flux1)])
    #   flux2m = np.median(flux2[~np.isnan(flux2)])
    #   print(flux1m)
    #   nflux1 = R1**2*flux1/nflux[1] + R2**2*flux2m/(R1**2*flux1m + R2**2*flux2m)
    #   nflux2 = (R2**2)*flux2/nflux[1] + R1**2*flux1m/(R1**2*flux1m + R2**2*flux2m)
    #   plt.plot(wave, nflux[0], 'k')
    #   plt.plot(wave, nflux1, '--r', label=f'primary rv={rv1} km/s')
    #   plt.plot(wave, nflux2, '--b', label=f'secondary rv = {rv2} km/s')
    #   #plt.xlim(4900, 5075)
    #   #plt.ylim(0.8, 1.05)
    #   plt.legend()
    #   return flux, nflux1, nflux2
    return flux

def vacuum2air(w):
    return w / (1.0 + 2.735182e-4 + 131.4182 / w**2 + 2.76249e8 / w**4)

def sampling_uniform_in_velocity(wave_base, wave_top, velocity_step):
    """
    Create a uniformly spaced grid in terms of velocity:

    - An increment in position (i => i+1) supposes a constant velocity increment (velocity_step).
    - An increment in position (i => i+1) does not implies a constant wavelength increment.
    - It is uniform in log(wave) since:
          Wobs = Wrest * (1 + Vr/c)^[1,2,3..]
          log10(Wobs) = log10(Wrest) + [1,2,3..] * log10(1 + Vr/c)
      The last term is constant when dealing with wavelenght in log10.
    - Useful for building the cross correlate function used for determining the radial velocity of a star.
    """
    # Speed of light
    c = 299792.4580 # km/s
    #c = 299792458.0 # m/s

    ### Numpy optimized:
    # number of elements to go from wave_base to wave_top in increments of velocity_step
    i = int(np.ceil( (c * (wave_top - wave_base)) / (wave_base*velocity_step)))
    grid = wave_base * np.power((1 + (velocity_step / c)), np.arange(i)+1)

    # Ensure wavelength limits since the "number of elements i" tends to be overestimated
    wfilter = grid <= wave_top
    grid = grid[wfilter]

    ### Non optimized:
    #grid = []
    #next_wave = wave_base
    #while next_wave <= wave_top:
        #grid.append(next_wave)
        ### Newtonian version:
        #next_wave = next_wave + next_wave * ((velocity_step) / c) # nm
        ### Relativistic version:
        ##next_wave = next_wave + next_wave * (1.-np.sqrt((1.-(velocity_step*1000.)/c)/(1.+(velocity_step*1000.)/c)))                                                     

    return np.asarray(grid)


def cross_correlation_function_uniform_in_velocity(wave_obs, flux_obs, flux_err, wave_template, flux_template,\
                                                   lower_velocity_limit, upper_velocity_limit, velocity_step, fourier=True):
    """
    Calculates the cross correlation value between the spectrum and the specified mask
    by shifting the mask from lower to upper velocity.

    - The spectrum and the mask should be uniformly spaced in terms of velocity (which
      implies non-uniformly distributed in terms of wavelength).
    - The velocity step used for the construction of the mask should be the same
      as the one specified in this function.
    - The lower/upper/step velocity is only used to determine how many shifts
      should be done (in array positions) and return a velocity grid.

    If fourier is set, the calculation is done in the fourier space. More info:

        VELOCITIES FROM CROSS-CORRELATION: A GUIDE FOR SELF-IMPROVEMENT
        CARLOS ALLENDE PRIETO
        http://iopscience.iop.org/1538-3881/134/5/1843/fulltext/205881.text.html
        http://iopscience.iop.org/1538-3881/134/5/1843/fulltext/sourcecode.tar.gz
    returns:
    -----------------
    velocities, ccf, ccf_err, len(flux), xfilter
    """
    last_reported_progress = -1
    # Speed of light in m/s
    c = 299792458.0

    # 1 shift = 1.0 km/s (or the specified value)
    shifts = np.arange(np.int32(np.floor(lower_velocity_limit)/velocity_step), np.int32(np.ceil(upper_velocity_limit)/velocity_step)+1)
    velocity = shifts * velocity_step

    waveobs = sampling_uniform_in_velocity(np.min(wave_obs), np.max(wave_obs), velocity_step)
    flux = np.interp(waveobs, wave_obs, flux_obs, left=0.0, right=0.0)
    flux = flux - np.median(flux)
    err = np.interp(waveobs, wave_obs, flux_err, left=0.0, right=0.0)


    depth = flux_template #np.abs(np.max(mask['flux']) - mask['flux'])
    resampled_mask = np.interp(waveobs, wave_template, depth, left=0.0, right=0.0)
    
    resampled_mask = resampled_mask - np.median(resampled_mask)
    if fourier:
        # Transformed flux and mask
        tflux = fft(flux)
        tresampled_mask = fft(resampled_mask)
        #conj_tflux = np.conj(tflux)
        conj_tresampled_mask = np.conj(tresampled_mask)
        num = np.int(len(resampled_mask)/2+1)
        #F_flux = abs(ifft(tflux*conj_tflux))
        #F_temp = abs(ifft(tresampled_mask*conj_tresampled_mask))
        tmp = abs(ifft(tflux*conj_tresampled_mask))
        #tmp = ifft(tflux*conj_tresampled_mask)
        ccf = np.hstack((tmp[num:], tmp[:num]))

        # Transformed flux and mask powered by 2 (second)
        #ccf_err = np.zeros(len(ccf))
        # Conservative error propagation
        terr = fft(err)
        tmp = abs(ifft(terr*conj_tresampled_mask))
        ccf_err = np.hstack((tmp[num:], tmp[:num]))
        ## Error propagation
        #tflux_s = fft(np.power(flux, 2))
        #tresampled_mask_s = fft(np.power(resampled_mask, 2))
        #tflux_err_s = fft(np.power(err, 2))
        #tresampled_mask_err_s = fft(np.ones(len(err))*0.05) # Errors of 5% for masks

        #tmp = abs(ifft(tflux_s*np.conj(tresampled_mask_err_s)))
        #tmp += abs(ifft(tflux_err_s*np.conj(tresampled_mask_s)))
        #ccf_err = np.hstack((tmp[num:], tmp[:num]))
        #ccf_err = np.sqrt(ccf_err)
        # Velocities
        
        sigma_s = np.sqrt(np.sum(flux*flux))
        sigma_t = np.sqrt(np.sum(resampled_mask*resampled_mask))
        velocities = velocity_step * (np.arange(len(resampled_mask), dtype=float)+1 - num)
        ccf = ccf/sigma_s/sigma_t
        ccf_err = ccf_err/sigma_s/sigma_t

    else:
        sigma_s = np.sqrt(np.sum(flux**2))        #Zucker, S. (2003) mnrs
        sigma_t = np.sqrt(np.sum(resampled_mask**2))        #Zucker, S. (2003) mnrs
        num_shifts = len(shifts)
        # Cross-correlation function
        ccf = np.zeros(num_shifts)
        ccf_err = np.zeros(num_shifts)
        velocities = velocity
        for shift, i in zip(shifts, np.arange(num_shifts)):
            #shifted_mask = resampled_mask
            if shift == 0:
                shifted_mask = resampled_mask
            elif shift > 0:
                #shifted_mask = np.hstack((shift*[0], resampled_mask[:-1*shift]))
                shifted_mask = np.hstack((resampled_mask[-1*shift:], resampled_mask[:-1*shift]))
            else:
                #shifted_mask = np.hstack((resampled_mask[-1*shift:], -1*shift*[0]))                                                                                      
                shifted_mask = np.hstack((resampled_mask[-1*shift:], resampled_mask[:-1*shift]))
            #ccf[i] = np.correlate(flux, shifted_mask)[0]
            #ccf_err[i] = np.correlate(err, shifted_mask)[0] # Propagate errors
            ccf[i] = np.sum(flux*shifted_mask)
            ccf_err[i] = np.sum(err*shifted_mask) # Propagate errors
            #ccf[i] = np.average(np.tanh(flux*shifted_mask))
            #ccf_err[i] = np.average(np.tanh(err*shifted_mask)) # Propagate errors

            current_work_progress = ((i*1.0)/num_shifts) * 100
            if report_progress(current_work_progress, last_reported_progress):
                last_reported_progress = current_work_progress
                #logging.info("%.2f%%" % current_work_progress)
                if frame is not None:
                    frame.update_progress(current_work_progress)
        ccf = ccf/sigma_s/sigma_t#resampled_mask_mod # Normalize
        ccf_err = ccf_err/sigma_s/sigma_t # Propagate errors

    
    
    

    # Filter to area of interest
    xfilter = np.logical_and(velocities >= lower_velocity_limit, velocities <= upper_velocity_limit)
    ccf = ccf[xfilter]
    ccf_err = ccf_err[xfilter]
    velocities = velocities[xfilter]

    return velocities, ccf, ccf_err, len(flux), xfilter
