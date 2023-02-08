from astropy.time import Time, units
from PyAstronomy import pyasl
from astropy import coordinates as coord
def eval_ltt(ra=62.794724868, dec=50.7082235439, jd=2456326.4583333, site=None, kind='barycentric', barycorr=True):                                                                                                            
    """ evaluate the jd 
    parameters
    ----------
    ra, dec: the coordinate of object
    jd: the julian date of observation (UTC time)
    site: the site of observatory
    return
    ------
    jd_llt: the bjd of hjd
    ltt: light_travel_time
    if barycorr=True return the barycentric correction
    # conf: https://docs.astropy.org/en/stable/time/
    # defaut site is Xinglong
    # coord.EarthLocation.from_geodetic(lat=26.6951*units.deg, lon=100.03*units.deg, height=3200*units.m) lijiang
    """
    if site is None:
        site = coord.EarthLocation.of_site('Beijing Xinglong Observatory')
    # sky position
    ip_peg = coord.SkyCoord(ra, dec, unit=(units.deg, units.deg), frame='icrs')
    # time
    times = Time(jd, format='jd', scale='utc', location=site)
    # evaluate ltt
    ltt = times.light_travel_time(ip_peg,kind)
    jd_llt = times.utc + ltt
    if barycorr is True:
       _barycorr = ip_peg.radial_velocity_correction(obstime=Time(times.iso), location=site)
       return jd_llt.jd, ltt, _barycorr
    else:
       return jd_llt.jd, ltt


def rv2baryrv(rv, barycorr):
    '''baryrv =rv+barycorr(1+rv/c) (km/s)
    parameters:
    --------------
    rv: [float] measured radial velocity (km/s)
    barycorr: [float] barycentric radial velocity correction
    returns:
    baryrv: [float] barycentric radial velocity
    '''
    c = 299792.458
    baryrv = rv+barycorr*(1+rv/c)
    return baryrv
