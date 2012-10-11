"""
SkyPos

This 'class' holds a collection of functions to calculate various quantities relevant to position in the sky.
"""

import numpy

# Conversion constants (radians - degrees)
_deg2rad = numpy.pi/180.0
_rad2deg = 180.0/numpy.pi
# MJD of date Jan 1, 2000.00
_MJD2000 = 51544.5


class SkyPos():

    def __init__(self):
        """Instantiate object."""
        return

    def setSite(self, lat=None, lon=None, height=None, pressure=None,
                 temperature=None, relativeHumidity=None, lapseRate=None):
        """Set parameters for the observatory site.
        If 'None' values are given for any parameter, the LSST default values are used.
        Lat/lon are in DEGREES."""
        if lat == None:
            self.siteLat = -29.6667
        else:
            self.siteLat = lat
        if lon ==None:
            self.siteLon = -70.59
        else:
            self.siteLon = lon
        if height==None:
            self.siteHeight= 2650.
        else:
            self.siteHeight = height
        if pressure==None:
            self.sitePressure = 749.3
        else:
            self.sitePressure = pressure
        if temperature==None:
            self.siteTemperature = 285.655
        else:
            self.siteTemperature = temperature
        if relativeHumidity==None:
            self.siteRelativeHumidity=0.40
        else:
            self.siteRelativeHumidity = relativeHumidity
        return

    def wrapRADec(self, ra, dec):
        """Set the RA/Dec values of the pointings. RA and Dec should be in DEGREES."""
        # Check for Dec values beyond -90 / +90.
        # Values larger than 90 are interpreted as being 'on the other side of the sphere'
        condition = (dec > 360)
        dec[condition] = dec[condition] % 360.
        condition ( dec < -360)
        dec[condition] = -1 * (-1*dec[condition] % 360.)
        condition = (dec > 90.)
        dec[condition] = 180. - dec[condition]
        ra[condition] = ra[condition] + 180.
        condition = (dec < -90.)
        dec[condition] = -180. - dec[condition]
        ra[condition] = ra[condition] + 180.
        # Check for RA values larger than 360 degrees or less than 0.
        ra = ra % 360.
        return ra, dec
        
    def calcLMST(self, mjd):
        # Calculate LMST. Assume we're doing this fairly simply, ignoring precession, etc.
        # Also we're assuming all times are UT1 and MJD units (and not TAI .. differences are small though).
        # Note Longitude here is increasing to the east.        
        jc = (mjd - _MJD2000) / 36525.0
        # Calculate gmst, using jc = the centuries since 2000.0
        gmst = (mjd % 1.0) * 360.0 + \
            (24110.54841 + jc*(8640184.812866 + jc*(0.093104 - jc*6.2e-6)))*0.0041666666666666667
        # Wrap into 360 degrees
        gmst = gmst % 360
        # And offset to the location of the site. (siteLon = increasing to the east). 
        lmst = gmst + self.siteLon
        # Wrap again into 360 degrees.
        lmst = lmst % 360
        return lmst

    def radec2altaz(self, ra, dec, mjd):
        """Given the time info to apply to the RA/Dec & site values stored in the class,
        calculate Alt/Az values.
        Uses formula from Pat Wallace (Starlink) for alt/az, and formula from Russell Owen for gmst.
        The time value can be a single MJD value, or a numpy array of the same length as ra/dec. """
        if isinstance(mjd, numpy.ndarray):
            if len(mjd) != len(ra):
                raise Exception('MJD can either be a single value or a numpy array of the same length as ra/dec')
        lmst = self.calcLMST(mjd)
        # Calculate Hour angle.
        hourangle = (lmst - ra) * _deg2rad
        # Calculate alt/az.
        sin_ha = numpy.sin(hourangle)
        cos_ha = numpy.cos(hourangle)
        sin_dec = numpy.sin(dec*_deg2rad)
        cos_dec = numpy.cos(dec*_deg2rad)
        sin_lat = numpy.sin(self.siteLat*_deg2rad)
        cos_lat = numpy.cos(self.siteLat*_deg2rad)
        x = -1. * cos_ha * cos_dec * sin_lat + sin_dec*cos_lat
        y = -1. * sin_ha * cos_dec
        z = cos_ha * cos_dec * cos_lat + sin_dec * sin_lat
        r = numpy.sqrt(x*x+y*y)
        a = numpy.where(r!=0, numpy.arctan2(y, x), 0.0)
        a = numpy.where(a<0, a+numpy.pi*2.0, a)
        az = a * _rad2deg
        alt = numpy.arctan2(z, r) * _rad2deg
        return alt, az
    
    def alt2airmass(self, alt):
        """Given altitude of field, calculate airmass.
        Uses Krisciunas and Schaefer formula for airmass. """
        # Calculate airmass.
        sz = (90.0 - alt) * _deg2rad
        airmass = numpy.power((1 - 0.96*(numpy.sin(sz))**2), -0.5)
        airmass = numpy.where(alt<0, 5.0, airmass)
        return airmass

    def calcDist_vincenty(self, ra1, dec1, ra2, dec2):
        """Calculates distance on a sphere using the Vincenty formula.
        Give this function RA/Dec values in degrees. Returns angular distance(s) in degrees."""
        r1 = ra1*_deg2rad
        r2 = ra2*_deg2rad
        d1 = dec1*_deg2rad
        d2 = dec2*_deg2rad
        v1 = (numpy.cos(d2)*numpy.sin(r2-r1))**2 + \
            (numpy.cos(d1)*numpy.sin(d2) - \
             numpy.sin(d1)*numpy.cos(d2)*numpy.cos(r2-r1))**2
        v1 = numpy.sqrt(v1)
        v2 = (numpy.sin(d1)*numpy.sin(d2) + \
              numpy.cos(d1)*numpy.cos(d2)*numpy.cos(r2-r1))
        dist = numpy.arctan2(v1,v2) * _rad2deg
        return dist

    def calcNight(self, mjd, midnight=0.16):
        """Given a time (in MJD/UTC), and the time at midnight at the site (in lsstsite),
        calculate and return the truncated 'day'. """
        # Note that this could be applied to an array of mjd's if needed.
        night = numpy.floor(mjd - midnight - 0.5)
        return night


