from astropy.io import fits
import numpy as np
from pyrafspec.fitslist import *
import matplotlib.pyplot as plt
from pyrafspec import splicing_spectrum
from pyrafspec.bfosclog import *
import pyrafspec.bfosce9g10 as bfosce9g10
import os
from astropy.time import Time
from astropy import constants
from PyAstronomy import pyasl
from astropy import coordinates as coord
from pyrafspec.spec_tools import rvcorr_spec
from laspec import normalization as lanorm
import joblib
import collections
from astropy.table import Table
plt.style.use('lijiaostyle')

# raw data directory
dire = '/home/lcq/media/backup/216BFOSC/sd0B_216_202109/20211114_bfosc/'

# log file
logfile = os.path.join(dire, '20211114_New.log')

###-------- the data observed after 20211023
starlist, lamplist = match_star2lamp(logfile, equipment='G10_E9', lamp_expt= 300, fout=None)

###-------- the data observed before 20211023
#starlist, lamplist = match_star2lamp_2020(logfile, equipment='G10_E9', lamp_expt= 300, fout=None)

def funcnorm(wave, flux, show=False):
    flux_norm, flux_smoothed =lanorm.normalize_spectrum_spline(wave, flux, p=1E-6, q=0.5, lu=(-1, 3), binwidth=150, niter=5)
    if show:
       plt.plot(wave, flux)
       plt.plot(wave, flux_smoothed)
    return flux_norm, flux_smoothed


show= True
fig, ax = plt.subplots(1,1, figsize=[13,7])
for _i, flamp in enumerate(lamplist[2:3]):
    fstar = starlist[2:3][_i]
    print(flamp, fstar)
    fwave, fflux, starfit, lampfit = bfosce9g10.dumpname(flamp, fstar, dire)
    star = splicing_spectrum.combine_wave_flux(fwave, fflux)
    waves = star['wavelength']['wave_solu']
    fluxs = star['flux']['spec_extr']
    flux_errs = star['flux']['err_extr']
    logwave, flux, flux_err = splicing_spectrum.splicing_spectrum(fwave, fflux, R=3000, N=3, lam_start=3700, lam_end=8900, 
                                                                  pix=[300, 2030], orders=np.arange(0, 11), funcnorm=funcnorm, 
                                                                   divide_blaze=False, threshold_blaze=700,
                                                                  show=show, ax=ax)
    for _i, wave in enumerate(waves):
        flux_norm, flux_smoothed = funcnorm(wave, fluxs[_i])
        plt.plot(wave[250:1950], flux_norm[250:1950], 'k', lw=0.5)

    
plt.plot(10**logwave, flux, color='r', label='splice')
plt.legend()
plt.xlim(3800, 9000)



'''
print('#-----------------splice the spectra--------------------------------')
for _i, flamp in enumerate(lamplist[2:3]):
    fstar = starlist[2:3][_i]
    print(flamp, fstar)
    fwave, fflux, starfit, lampfit = bfosce9g10.dumpname(flamp, fstar, dire)
    star = splicing_spectrum.combine_wave_flux(fwave, fflux)
    logwave, fluxnorm, fluxnorm_err = splicing_spectrum.splicing_spectrum(fwave, fflux, R=3000, N=3, lam_start=3700, 
                                                                          lam_end=8900, pix=[200, 2030], orders=np.arange(0, 11), 
                                                                          funcnorm=funcnorm, divide_blaze=True, threshold_blaze=100, ax=ax,
                                                                          sigma =8, percentile_up=99, window_length=11, show=False)
                                                                          
    table = fits.BinTableHDU.from_columns([
                    fits.Column(name='loglambda',  format='E', array=np.array(logwave, dtype=np.float32)),
                    fits.Column(name='flux',format='E', array=np.array(flux, dtype=np.float32)),
                    fits.Column(name='error',format='E', array=np.array(flux_err, dtype=np.float32))
                                         ])
    header = fits.getheader(starfit)
    data0 = np.zeros(1,dtype=np.float32)
    hdu = fits.PrimaryHDU(data0)
    hdu.header = header
    hdu.header['fstar'] = (fstar, 'fname of star raw image')
    hdu.header['flamp'] = (flamp, 'fname of lamp raw image')
    ##########################################################################################
    jd = header['JD'] + header['EXPOSURE']/2./3600./24.
    ##--------------light_travel_time----------------------------------------------------------
    ##radec is the coordinates of the objects e.g. '04:11:10.73 +50:42:29.60'
    radec = f"{header['RA']} {header['DEC']}"
    ip_peg = coord.SkyCoord(radec, unit=(units.hourangle, units.deg), frame='icrs')
    site = coord.EarthLocation.of_site('Beijing Xinglong Observatory')
    times_ltt = Time(jd, format='jd', scale='utc', location=site)
    ltt = times_ltt.light_travel_time(ip_peg,'barycentric')
    barycorr = ip_peg.radial_velocity_correction(obstime=Time(times_ltt.iso), location=site)
    #----------------------------------------------------------------------------------------
    header['BJD'] = jd+ltt.jd
    header['barycorr'] = (barycorr.to('km/s').value, 'km/s (rv=rv+barycorr+rv*barycorr/c)')
    ###########################################################################################
    hdul = fits.HDUList([hdu, table])
    fout = starfit[:-4] +'_splicing.fits'
    hdul.writeto(fout, overwrite=True)
print('#-----------------------splice end--------------------------------')
'''


## check lamp
fig, ax = plt.subplots(1,1, figsize=(13,6),sharex=True)
plt.xlabel('wavelength')
for _lampname in np.unique(lamplist):
    _lampname = f'{dire}/fear-{_lampname}.fit.dump'
    lampdata = joblib.load(_lampname)
    plt.plot(lampdata['wave_solu'].T, lampdata['fear1d'].T)
    
plt.show()