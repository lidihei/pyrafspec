from astropy.io import fits
import numpy as np
from pyrafspec.fitslist import *
import matplotlib.pyplot as plt
from pyrafspec import splicing_spectrum
from pyrafspec.bfosclog import *
import pyrafspec.bfosce9g10 as bfosce9g10
import os
from astropy.time import Time
from astropy import constants, units
from PyAstronomy import pyasl
from astropy import coordinates as coord
from pyrafspec.spec_tools import rvcorr_spec
from pyrafspec.lighttraveltime import rv2baryrv
from laspec import normalization as lanorm
import joblib
import collections
from astropy.table import Table
from tqdm import tqdm
from pyrafspec.RVsersic import *


# laod and check radial velocity model
rvm = joblib.load('/media/share/lijiao/TLUSTY_grid/tlusty_for_bfosc_sdB_rv.dump')
from laspec.ccf import RVM
pmod_rv = rvm.pmod
wavemod_rv = rvm.wave_mod
fluxmod_rv = rvm.flux_mod
rvm = RVM(pmod_rv, wavemod_rv, fluxmod_rv, npix_lv=5)

#### check rvm model wave is in vacuum or air
#haVac = 6564.66464 # vacuum wavelength of Halpha
#haAir = 6562.85175
#heIAir = 6678.151
#heIVac = 6679.995 # vacuum wavelength of HI
#fig, axs = plt.subplots(2,1,figsize=[10,8])
#plt.sca(axs[0])
#plt.plot(rvm.wave_mod, rvm.flux_mod[0])
#plt.axvline(x=4862.6) #vacum
#plt.axvline(x=4923.3)#vacuum
#plt.axvline(x=5017.07)#vacuum
#plt.xlabel('wavelength')
#plt.xlim(4861, 5020)
#
#plt.sca(axs[1])
#plt.plot(rvm.wave_mod, rvm.flux_mod[0])
#plt.axvline(x=haVac)
#plt.axvline(x=heIVac)
#plt.xlabel('wavelength')
#plt.xlim(6500, 6700)

def getrvind(wave, flux):
    ind = ~( (flux > 1.05) |
        (flux < 0.) |
        (wave < 3940) |
        #((wave >5881) & (wave < 6300))|
        ((wave > 6800) & (wave < 7000)) |
        ((wave > 6600) & (wave < 6640))|
        (wave > 7100) |
        np.isnan(flux)
       )
    return ind

def get_rvs(fname, verbose=True, show=False):
    hdul = fits.open(fname)
    header = hdul[0].header
    data = hdul[1].data
    data0 = np.zeros(1,dtype=np.float32)
    flamp = header['FLAMP']
    rms_lamp = header['RMS_LAMP']
    #wave = 10**data['loglambda']
    wave = 10**data['loglambda_bary']
    flux = data['flux']
    fluxerr = data['error']
    barycorr = header['BARYCORR']
    ra = header['RA']
    dec = header['DEC']
    bjd = header['bjd']
    radec = f"{header['RA']} {header['DEC']}"
    ip_peg = coord.SkyCoord(radec, unit=(units.hourangle, units.deg), frame='icrs')
    ra_deg = np.float32(ip_peg.ra.deg)
    dec_deg = np.float32(ip_peg.dec.deg)
    ########------------CCF rv-------------------------------------------########################
    ind = getrvind(wave, flux)
    wave_obs, flux_obs, fluxerr_obs = wave[ind],  flux[ind], fluxerr[ind]
    rvr = rvm.measure(wave_obs, flux_obs, flux_err=fluxerr_obs, nmc=100, rv_grid=np.arange(-500, 500, 5),
                        )
    rv_ccf = np.float32(rvr['rv_opt'])
    err_ccf = np.float32(np.sqrt(np.sum(np.diff(rvr['rv_pct'])**2)/2))
    basename = os.path.basename(fname)
    if verbose: print(f'{basename}:\nccf baryrv = {rv_ccf}\pm {err_ccf}')
    hdul[0].header['rv_ccf'] = (rv_ccf, '(km/s) barycentric ccf rv')
    hdul[0].header['rverr_ccf'] = (err_ccf, 'radial velocity error(km/s)' )
    ########------------SERSIC rv-------------------------------------------########################
    rv_halpha,rverr_halpha = rv_bysersic(wave, flux, lam0 = 6562.79, waverange = [6550, 6580], p0 = [0.8,4.,1.,6562.], show=show)
    rv_halpha,rverr_halpha = np.float32(rv_halpha), np.float32(rverr_halpha)
    rv_hbeta, rverr_hbeta = rv_bysersic(wave, flux, lam0=4861.35, waverange=[4840, 4890], p0=[0.6,4.5,0.8,4862.], show=show)
    rv_hbeta,rverr_hbeta = np.float32(rv_hbeta), np.float32(rverr_hbeta)
    rv_hgamma,rverr_hgamma = rv_bysersic(wave, flux, lam0=4340.46, waverange = [4310, 4360], p0 = [0.8,4.,1.,4340.], show=show)
    rv_hgamma,rverr_hgamma = np.float32(rv_hgamma), np.float32(rverr_hgamma)
    if verbose: print(f'sersic fitting baryrv: Halpha = {rv_halpha}\pm{rverr_halpha}')
    if verbose: print(f'sersic fitting baryrv: Hbeta={rv_hbeta}\pm{rverr_hgamma}')
    if verbose: print(f'sersic fitting baryrv: Hgamma={rv_hgamma}\pm{rverr_hgamma}')
    hdul[0].header['rv_ha'] = (rv_halpha, '(km/s) barycentric Halpha rv')
    hdul[0].header['rverr_ha'] = (rverr_halpha, '(km/s) barycentric Halpha rv err')
    hdul[0].header['rv_hb'] = (rv_hbeta, '(km/s) barycentric Hbeta rv')
    hdul[0].header['rverr_hb'] = (rverr_hbeta, '(km/s) barycentric Hbeta rv err')
    hdul[0].header['rv_hg'] = (rv_hgamma, '(km/s) barycentric Hgamma rv')
    hdul[0].header['rverr_hg'] = (rverr_hgamma, '(km/s) barycentric Hgamma rv err')
    ###########################################################################################
    hdul.writeto(_fname, overwrite=True)
    _line = basename,bjd,ra,dec,ra_deg,dec_deg,flamp,rms_lamp,rv_ccf,err_ccf,rv_halpha,rverr_halpha,rv_hbeta,rverr_hbeta,rv_hgamma,rverr_hgamma
    _line = ','.join(list(np.array(_line)))
    _line = _line+'\n'
    return _line


from glob import glob

dires = ['20201215', '20201218','20201219','20210225','20210304','20210307','20210316','20210321','20210323','20210424','20210425',
         '20210504','20210505', '20210508', '20210601', '20210602', '20210603', '20211022', '20211023', '20211024', '20211113', '20211114']
for _dire in dires:
    print(_dire)
    dire = f'/home/lcq/media/backup/216BFOSC/sdOB_216/{_dire}_bfosc/'
    fitslists = glob(f'{dire}*_splicing.fits')

    lines = 'basename,bjd,ra,dec,ra_deg,dec_deg,flamp,rms_lamp,rv_ccf,err_ccf,rv_halpha,rverr_halpha,rv_hbeta,rverr_hbeta,rv_hgamma,rverr_hgamma\n'
    for _i, _fname in enumerate(tqdm(fitslists)):
        lines += get_rvs(_fname, verbose=False, show=False)
    tabname =  os.path.join(dire, 'rv_tab.csv')
    with open(tabname, 'w') as ff:
        ff.writelines(lines)
    ff.close()
