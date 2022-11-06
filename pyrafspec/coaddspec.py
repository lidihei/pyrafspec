import numpy as np
from astropy.io import fits
from laspec import normalization
from laspec.ccf import RVM
from .spec_tools import rvcorr_spec


class coaddspec():
    
    def coadd_spec(self, waves, fluxs, fluxerrs, rvs, wave_new, sigma=5, niter=3):
        '''
        waves: [2D list or  2D array] contians the wave of each spectrum, waves.shape = (the number of spectra, the ponit number of each spectrum)
        fluxs: [2D list or  2D array] contians the flux of each spectrum
        fluxerrs: [2D list or  2D array]
        rvs: [1D arrray] the radial velocities of each spectrum
        wave_new: [1D array] the wavelength of the coadd spectra
        returns:
        ---------------
        flux_w: [1D array] flux coadded with weight $f = \frac{\sum{w_i*f_i(\lambda)}}{\sum{w_i}}$; $w_i = median(snr_i)$
        fluxerr_w: [1D array] error of flux which is coadded with weight
        flux_sum: [1D array] flux which is directly added
        fluxerr_sum: [1D array] error of flux which is directly added
        '''
        nflux = len(fluxs)
        fluxs_new = np.ma.zeros((nflux, len(wave_new)))
        fluxerrs_new = np.zeros_like(fluxs_new)
        fluxs_stdc = np.ones_like(fluxs_new)
        fluxs_medianc = np.ones_like(fluxs_new)
        weighti = np.ones_like(fluxs_new)
        mediani = np.ones_like(fluxs_new)
        mask = np.zeros_like(fluxs_new, dtype=np.bool)
        for _i in np.arange(nflux):
            _flux, _fluxerr = rvcorr_spec(waves[_i], fluxs[_i], fluxerrs[_i], -rvs[_i], wave_new=wave_new, left=np.nan, right=np.nan, interp1d=None, returnwvl=False)
            fluxs_new[_i] = _flux
            fluxerrs_new[_i] = _fluxerr
            weighti[_i] = np.median(_flux/_fluxerr)
            mediani[_i] = np.median(_flux)
        fluxs_median = fluxs_new/mediani
        for _ in np.arange(niter):
            median = np.median(fluxs_median, axis=0)
            std = np.std(fluxs_median, axis=0)
            for _i in np.arange(nflux):
                fluxs_medianc[_i] = median
                fluxs_stdc[_i] = std
            dflux = np.abs(fluxs_median -fluxs_medianc)
            mask = dflux > (sigma*fluxs_stdc)
            fluxs_median.mask = mask
        weighti.mask = mask
        weight = weighti.copy()
        fluxs_new.mask = mask
        fluxerrs_new.mask = mask
        weightc_sum =np.sum(weighti, axis=0)
        for _i in np.arange(nflux):
            weight[_i] = weighti[_i]/weightc_sum
        weighti.mask = mask
        weight2 = weight**2
        #flux_m = np.sum(fluxs_median*weight, axis=0)
        flux_w = np.sum(fluxs_new*weight, axis=0)
        fluxerr2_w = np.sum(fluxerrs_new**2*weight2, axis=0)
        fluxerr2_sum = np.sum(fluxerrs_new.data**2, axis=0)
        fluxerr_w = np.sqrt(fluxerr2_w)
        flux_sum = np.sum(fluxs_new.data, axis=0)
        fluxerr_sum = np.sqrt(fluxerr2_sum)
        self.waves = waves
        self.fluxs = fluxs
        self.fluxerrs = fluxerrs
        self.flux_w = flux_w
        self.fluxerr_w = fluxerr_w
        self.flux_sum = flux_sum
        self.fluxerr_sum = fluxerr_sum
        self.wave = wave_new
        return flux_w, fluxerr_w, flux_sum, fluxerr_sum

    def measurerv(self, waves=None, fluxs=None, fluxerrs=None, rvm=None, **keywords):
        waves = self.waves if waves is None else waves
        fluxs = self.fluxs if fluxs is None else fluxs
        fluxerrs = self.fluxerrs if fluxerrs is None else fluxerrs
        nflux = len(fluxs)
        if rvm is None:
           pmod_rv = ['coadd']
           _norm, _cont = normalization.normalize_spectrum_spline(self.wave, self.flux_sum, niter=2)
           _wave = np.arange(self.wave[0], self.wave[-1]+0.01, 0.01)
           #_flux = np.interp(_wave, self.wave, self.flux_sum/np.median(self.flux_sum))
           _flux = np.interp(_wave, self.wave, _norm)
           rvm = RVM(np.array([pmod_rv]), _wave, np.array([_flux]), npix_lv=0)
        rvs,rverrs = np.zeros((2,nflux))
        for _i, wave_obs in enumerate(waves):
            flux_obs = fluxs[_i]
            fluxerr_obs = fluxerrs[_i]
            rvr = rvm.measure(wave_obs, flux_obs, flux_err=fluxerr_obs, nmc=100, rv_grid=(-600, 600, 1),
                        )
            #rvs[_i] = np.float32(rvr['rv_opt'])
            rvs[_i] = np.float32(rvr['rv_pct'][1])
            rverrs[_i] = np.float32(np.sqrt(np.sum(np.diff(rvr['rv_pct'])**2)/2))
        self.rvm = rvm
        return rvs, rverrs
    
    def coadd_spec_iter(self, waves_rv, fluxs_rv, fluxerrs_rv, rvs, tolerance = 1, rvm=None):
       ''' iterately coadd the spectra by using radial velocity measured by the coadded spectra
       rvs [1D array] the initall radial velocities
       tolerance [float] stop interation while np.median(rvs_i - rvs_i-1) < tolerance
       rvm: [laspec.ccf.RVM]
       '''
        drv = 100
        rvs0 = rvs.copy()
        while drv > tolerance:
            rvs, rverrs = self.measurerv(waves=waves_rv, fluxs=fluxs_rv, fluxerrs=fluxerrs_rv, rvm=rvm)
            flux_w, fluxerr_w, flux_sum, fluxerr_sum = self.coadd_spec(self.waves, self.fluxs, self.fluxerrs, rvs, self.wave, sigma=5, niter=3)
            drv = np.median(np.abs(rvs - rvs0))
            print(drv)
            rvs0 = rvs.copy()
        self.rvs = rvs
        self.rverrs = rverrs
        self.flux_w = flux_w
        self.fluxerr_w = fluxerr_w
        self.flux_sum = flux_sum
        self.fluxerr_sum = fluxerr_sum
        
    def writecoaddspec(self, fout, wave=None, flux_w=None, fluxerr_w=None, flux_sum=None, fluxerr_sum=None):
        ''' write the coadd spectrun into fits file
        '''
        wave = self.wave if wave is None else wave
        flux_w = self.flux_w if flux_w is None else flux_w
        fluxerr_w = self.fluxerr_w if fluxerr_w is None else fluxerr_w
        flux_sum = self.flux_sum if flux_sum is None else flux_sum
        fluxerr_sum = self.fluxerr_sum if fluxerr_sum is None else fluxerr_sum
        data = [wave, flux_w, fluxerr_w, flux_sum, fluxerr_sum]
        names = ['wavelength', 'flux_w', 'fluxerr_w', 'flux_sum', 'fluxerr_sum']
        coaddspec = Table(data=data, names=names)
        coaddspec.write(fout, overwrite=True)
        hdu = fits.open(fout)
        header = hdu[1].header
        header['TCOMM2'] = 'flux is coadded with weight [median(snr_i)]'
        header['TCOMMe'] = 'error of flux coadded with weight'
        header['TCOMM4'] = 'flux is coadded directly '
        header['TCOMM5'] = 'error of flux coadded directly'
        hdu[1].header = header
        hdu.writeto(fout, overwrite=True)
        
