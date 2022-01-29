import numpy as np
#import scipy.linalg as splin
import numpy.linalg as nl
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('xtick',labelsize=12)
matplotlib.rc('ytick',labelsize=12)

from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.optimize import curve_fit
import astropy.io.fits as fits

import scipy.special as special

from matplotlib.ticker import MultipleLocator
from scipy.interpolate import UnivariateSpline

import re
from mechanize import Browser
from PIL import Image as PILImage
from PIL import ImageOps
try:
   import StringIO
except:
   from io import StringIO, BytesIO
import base64
from astropy.coordinates import SkyCoord
from astropy import units as u
from IPython.display import display, Image, HTML
from astropy.table import Table,vstack,hstack

import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy import units

from astropy import coordinates as coord

def phase2obstime(phase, localdate, period=3.3893, t0 = 2458819.8, ra=62.794724868, dec=50.7082235439, dtime=0, site=None, kind='barycentric', prt=True):
    '''calculate obsertime time by a given phase
    phase [float] [0, 1]e.g. phase= 0.5
    localdate [utc] e.g. '2021-01-02T16:00:00' the local observation date
    period [float] orbital period in day
    t0 [float] hjd or bjd, is the time of zero phase
    ra [float] the right  Ascension of observed object
    dec [float] the declination of observed objec
    dtime [float] is the intervel hours between local time and utc time, e.g. dtime=-8 for Beijing time
    site #defaut site is Xinglong
    mp [str] '-/+' if '-' minus the difference time 
    prt [bool] if True print phase and observation time
    '''
    _localtime = Time(localdate,  format='isot', scale='utc')
    _utctime =  _localtime + dtime*units.hour
    jd_llt, ltt = eval_ltt(ra=ra, dec=dec, jd=_utctime.jd, site=site, kind=kind)
    _phase = np.mod(jd_llt-t0, period)/period
    print(f'phase of the {localdate}: {_phase}')
    phase = np.asarray(phase)
    dphase = phase - _phase
    #if mp is '-':
    #   obstime = Time(jd_llt - dphase*period, format='jd', scale='utc')
    #else:
    obstime = Time(jd_llt+dphase*period, format='jd', scale='utc') 
    utctime = obstime - ltt
    localtime = utctime -  dtime*units.hour
    if prt:
       for _i, ph in enumerate(phase):
           print(f'phase {ph}, localtime {localtime.iso[_i]}')
    return localtime

def eval_ltt(ra=62.794724868, dec=50.7082235439, jd=2456326.4583333, site=None, kind='barycentric'):                                                                                                            
    """ evaluate the jd """
    # conf: https://docs.astropy.org/en/stable/time/
    # defaut site is Xinglong
    if site is None:
        site = coord.EarthLocation.of_site('Beijing Xinglong Observatory')
    # sky position
    ip_peg = coord.SkyCoord(ra, dec, unit=(units.deg, units.deg), frame='icrs')
    # time
    times = Time(jd, format='jd', scale='utc', location=site)
    # evaluate ltt
    ltt = times.light_travel_time(ip_peg,kind)
    jd_llt = times.utc + ltt
    return jd_llt.jd, ltt

class observation():
     def flip_rotate_image(self, im=None, FLIP_LEFT_RIGHT=False,degrees_to_rotate=0, inverse=False):
         '''rotate image and mirror image by y-axis
         '''
         if im is None:
            im = self.im
         fig = plt.figure(figsize=[16,16])
         ax = fig.add_subplot(111)
         if FLIP_LEFT_RIGHT:
            im = im.transpose(PILImage.FLIP_LEFT_RIGHT)
         im = im.rotate(degrees_to_rotate)
         if inverse:
             ax.imshow(ImageOps.invert(im))
         else:
             ax.imshow((im))
         #plt.xlabel('RA',fontsize=18)
         #plt.ylabel('DEC', fontsize=18)
         plt.show()
     
     def downloadChart(self, ra,dec, mag, num=1, fw='6', fh='6', inverse=False,FLIP_LEFT_RIGHT=False,degrees_to_rotate=0, ):
         '''download chart
         parameters
         ----------
         ra [float], degree
         dec [float], degree
         fw [str], field_width
         fh [str], field_high
         FLIP_LEFT_RIGHT [bool], flip left to right of image (mirror image by y-axis)
         degrees_to_rotate [float], in units of degree
         '''
         cc = SkyCoord(ra=ra*u.degree, dec=dec*u.degree,frame='fk5')
         br=Browser()
         br.open("http://astro.swarthmore.edu/finding_charts.cgi")
         br.select_form(nr=0)
         rahms = cc.ra.hms
         decdms = cc.dec.dms
         rastr = "%(h)02d %(m)02d %(s).2f" % {'h':rahms.h,'m':np.int32(np.abs(rahms.m)),'s':np.abs(rahms.s)}
         destr = "%(d)+03d %(m)02d %(s).2f" % {'d':decdms.d,'m':np.int32(np.abs(decdms.m)),'s':np.abs(decdms.s)}
         #rastr0 = "{0:02d}{1:02d}{2:.2f}".format(np.int32(rahms.h),np.int32(np.abs(rahms.m)), np.abs(rahms.s))
         #destr0 = "{0:+03d}{1:02d}{2:.2f}".format(np.int(decdms.d),np.int32(np.abs(decdms.m)), np.abs(decdms.s))
         #print(rastr,destr)
         br["ra"]= rastr
         br["dec"]=destr
         br["field_width"]=fw
         br["field_height"]=fh
         res = br.submit()
         desig = f"J{rastr.replace(' ', '')}{destr.replace(' ', '')}"
         
         name = '{0:02d}: desig={1}; ({2} {3}); Gmag={4:.3f}'.format(ii, desig,ra, dec,Gmag)
         name = 'desig={1}; ({2} {3})'.format(ii, desig,ra, dec,Gmag)
         fig = plt.figure(figsize=[16,16])
         ax = fig.add_subplot(111)
         im = PILImage.open(BytesIO(res.get_data()))
         if FLIP_LEFT_RIGHT:
            im = im.transpose(PILImage.FLIP_LEFT_RIGHT)
         im = im.rotate(degrees_to_rotate)
         if inverse:
             ax.imshow(ImageOps.invert(im))
         else:
             ax.imshow((im))
         plt.title(name, fontsize=22)
         plt.xlabel('RA',fontsize=18)
         plt.ylabel('DEC', fontsize=18)
         plt.show()
         self.im = im
         return im

     def mag2obstime(self, mag, mag_ref, exposure_ref):
         '''calculate observation time by a reference star
         parameters
         ------------
         mag [float] the magnitude of observation star
         mag_ref [float] the magnitude of reference star
         exposure_ref [float] the exposure time of reference star
         returns
         ---------------
         ratio [float] the time ratio between observation star and reference star
         exposure_obs [float] the exposure time of object
         '''
         mag = np.asarray(mag)
         ratio = 10**((mag-mag_ref)/2.5)
         exposure_obs = exposure_ref*ratio
         return ratio, exposure_obs

     def snr2obstime(self, snr, snr_ref, exposure_ref):
         '''calculate observation time of SNR by given a reference SNR
         parameters
         ------------
         snr [float] the expected SNR of observation star
         snr_ref [float] the reference SNR 
         exposure_ref [float] the exposure time of reference star
         returns
         ---------------
         exposure_obs [float] the exposure time of object
         '''
         flux = snr**2
         flux_ref = snr_ref**2
         exposure_obs = exposure_ref *flux/flux_ref
         return exposure_obs
