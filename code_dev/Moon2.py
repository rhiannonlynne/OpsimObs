"""
Moon

Class to calculate the approximate position of the Moon at a various times.
"""

import numpy
from Sun import Sun

_deg2rad = numpy.pi/180.0
_rad2deg = 180.0/numpy.pi
_MJD2000 = 51544.5


class Moon():

    def __init__(self):
        """Instantiate the Sun object."""
        return

    def calcPos(self, mjd):
        """Calculate the moon's ecliptic lon/lat and geocentric RA/Dec, for the given MJD(s)."""
        # Using the procedure outlined http://www.stjarnhimlen.se/comp/ppcomp.html 
        self.mjd = numpy.copy(mjd)
        # Calculate the days since J2000.
        days = self.mjd - _MJD2000
        # Calculate the obliquity of the ecliptic at this/these times.
        ecl = (23.4393 - 3.563E-7 * days) * _deg2rad
        # Calculate the moon's orbital elements. (all angles in radians)
        N = ((125.1228 - 0.0529538083 * days) % 360.0) * _deg2rad
        i = (5.1454) * _deg2rad
        w = ((318.0634 + 0.1643573223 * days) % 360.0) * _deg2rad
        a = 60.2666  # earth radii
        e = 0.054900
        M = ((115.3654 + 13.0649929509 * days) % 360.0) * _deg2rad
        # Calculate the eccentric anomaly.
        E = M + e * numpy.sin(M) * ( 1.0 + e * numpy.cos(M) )
        # Calculate the moon's distance and true anomaly.
        xv = a * (numpy.cos(E) - e )
        yv = a * (numpy.sqrt(1.0 - e*e) * numpy.sin(E))
        self.dist = numpy.sqrt(xv*xv + yv*yv)
        v = numpy.arctan2(yv, xv)
        # Calculate the position in 3-d space.
        xh = self.dist * (numpy.cos(N) * numpy.cos(v+w) - numpy.sin(N) * numpy.sin(v+w) * numpy.cos(i))
        yh = self.dist * (numpy.sin(N) * numpy.cos(v+w) + numpy.cos(N) * numpy.sin(v+w) * numpy.cos(i))
        zh = self.dist * (numpy.sin(v+w) * numpy.sin(i) )
        # Calculate the position in the ecliptic plane.
        self.ec_lon = numpy.arctan2(yh, xh) * _rad2deg
        self.ec_lat = numpy.arctan2(zh, numpy.sqrt(xh*xh+yh*yh)) * _rad2deg
        # Add perturbations to the moon's position. (needed for accuracy < 2 deg).
        # Calculate the Sun's Mean Anomaly. (radians)
        Msun = ((356.0470 + 0.9856002585 * days) % 360.0) * _deg2rad
        # Calculate the Sun's argument of perihelion. (radians)
        wsun = ((282.9404 + 4.70935E-5 * days) % 360.0) * _deg2rad
        # Calculate the mean longitude of the Sun (Nsun = 0).  (.. radians...)
        Lsun = Msun + wsun
        # Calculate the mean longitude of the Moon.
        Lmoon = M + w + N
        # Calculate the mean elongation of the Moon.
        D = Lmoon - Lsun
        # Calculate the argument of latitude for the Moon.
        F = Lmoon - N
        # Calculate terms to add to longitude of the moon. (degrees)
        LonCorr = (-1.274 * numpy.sin(M - 2*D) + 0.658 * numpy.sin(2*D) - 0.186 * numpy.sin(Msun) \
                   -0.059 * numpy.sin(2*M - 2*D) - 0.057 * numpy.sin(M-2*D+Msun) + 0.053 *numpy.sin(M + 2*D) \
                   +0.046 * numpy.sin(2*D - Msun) + 0.041 * numpy.sin(M - Msun) - 0.035 * numpy.sin(D) \
                   -0.031 * numpy.sin(M + Msun) - 0.015 * numpy.sin(2*F - 2*D) + 0.011 * numpy.sin(M - 4*D))
        self.ec_lon = self.ec_lon + LonCorr
        print self.ec_lon.min(), self.ec_lon.max()
        # Calculate terms to add to latitude of the Moon. (degrees)
        LatCorr = (-0.173 * numpy.sin(F - 2*D) - 0.055 * numpy.sin(M - F - 2*D) - 0.046 * numpy.sin(M + F - 2*D) \
                   + 0.033 * numpy.sin(F + 2*D) + 0.017 * numpy.sin(2*M + F))
        self.ec_lat = self.ec_lat + LatCorr
        # Calculate terms to add to distance to the Moon.
        DistCorr = (-0.58 * numpy.cos(M - 2*D) - 0.46 * numpy.cos(2*D))
        self.dist = self.dist + DistCorr
        from pyslalib import slalib
        ra = numpy.zeros(len(self.mjd), 'float')
        dec = numpy.zeros(len(self.mjd), 'float')
        for i in range(len(self.mjd)):
            ra[i], dec[i] = slalib.sla_ecleq(self.ec_lon[i]*_deg2rad, self.ec_lat[i]*_deg2rad, self.mjd[i])
        self.ra = ra * _rad2deg
        self.dec = dec * _rad2deg
        return
        # Compute geocentric position in x/y/z.
        xh = self.dist * numpy.cos(self.ec_lon) * numpy.cos(self.ec_lat)
        yh = self.dist * numpy.sin(self.ec_lon) * numpy.cos(self.ec_lat)
        zh = self.dist * numpy.sin(self.ec_lat)
        # Convert to equatorial coordinates (rotate by ecliptic obliquity)
        xe = xh
        ye = yh * numpy.cos(ecl) - zh * numpy.sin(ecl)
        ze = yh * numpy.sin(ecl) + zh * numpy.cos(ecl)
        # Convert to ra/dec. (note that these are geocentric RA/Dec .. not topocentric).
        #  The difference between topocentric and geocentric can be as much as a degree due to parallax. 
        self.ra = numpy.arctan2(ye, xe) * _rad2deg
        self.dec = numpy.arctan2(ze, numpy.sqrt(xe*xe+ye*ye)) * _rad2deg
        self.dist = numpy.sqrt(xh*xh+yh*yh+zh*zh)
        # Calculate the lunar phase. 
        # Calculate the solar longitude. 
        sun = Sun()
        lon_sun = sun.getLon(self.mjd)
        # Calculate the solar elongation of the Moon.
        solarelong = numpy.arccos(numpy.cos((lon_sun - self.ec_lon)*_deg2rad) * numpy.cos(self.ec_lat*_deg2rad))
        # Calculate the phase of the moon.
        self.phase = 180.0 - solarelong*_rad2deg
        # Calculate the illumination of the Moon. 
        self.illum = ( 1 + numpy.cos(self.phase*_deg2rad) )/2.0
        return
    
    def getAltAz(self, skypos):
        """Get Alt/Az/Airmass information for each observation. Pass an already instantiated SkyPos object."""
        # Calculate the alt/az for these times, 
        self.alt, self.az = skypos.radec2altaz(self.ra, self.dec, self.mjd)
        # Add a correction for the geocentric -> topocentric position. 
        self.parallax = numpy.arcsin(1/self.dist) * _rad2deg
        self.alt = self.alt - self.parallax * numpy.cos(self.alt)
        print self.alt
        return

    def test_slalib(self, mjd, config):
        from pyslalib import slalib
        ra = numpy.zeros(len(mjd), 'float')
        dec = numpy.zeros(len(mjd), 'float')
        for i in range(len(mjd)):
            # Calculate moon's position. 
            #(ra_RAD, dec_RAD, diam) = slalib.sla_rdplan(mjd[i], 3,
            #                                            config['longitude']*_deg2rad, config['latitude']*_deg2rad)
            moon_posvel = slalib.sla_dmoon(mjd[i])
            moon_sph = slalib.sla_dc62s(moon_posvel)
            ra_RAD = moon_sph[0]
            dec_RAD = moon_sph[1]
            ra[i] = ra_RAD * _rad2deg
            dec[i] = dec_RAD * _rad2deg
        print "Moon Test_slalib results :"
        print " .. ra(slalib)"
        print ra
        print " .. ra(self)"
        print self.ra
        print self.ra.min(), self.ra.max()
        print " .. dec(slalib)"    
        print dec
        print " .. dec(self)"
        print self.dec
        
