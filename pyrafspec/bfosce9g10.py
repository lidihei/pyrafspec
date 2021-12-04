from astropy.io import fits
import os

def dumpname(flamp, fstar, dire, lampprefix='fear', objecprefix = 'star'):
    ''' get dump file from lamp and star file
    parameters:
    -------------
    flamp: [flamp] file name of lamp (e.g. '202110220021_SPECSLAMP_FeAr_slit16s_G10_E9')
    fstar: [str] file name of star spectrum (e.g. 202110220019_SPECSTARGET_BD+25d4655_slit16s_G10_E9)
    dire: [str] directory stored raw data and dump file
    returns:
    -------------
    fwave: [str] dump file name of lamp
    fflux: [str] dump file name of star
    fstarall: [str]
    flampall: [str]
    '''
    stardump = f'{objecprefix}-{fstar}.fit.dump'
    lampdump = f'{lampprefix}-{flamp}.fit.dump'
    fwave = os.path.join(dire, lampdump)
    fflux = os.path.join(dire, stardump)
    fstarall = os.path.join(dire, f'{fstar}.fit')
    flampall = os.path.join(dire, f'{fstar}.fit')
    return fwave, fflux, fstarall, flampall

def _write2fits(fstar, flamp, dire, fout=None):
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



