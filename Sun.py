"""
Sun

Class to calculate the approximate position of the sun at a particular time.

"""

import numpy

_deg2rad = numpy.pi/180.0
_rad2deg = 180.0/numpy.pi
_MJD2000 = 51544.5

class Sun():

    def __init__(self):
        """Instantiate the Sun object."""
        return


    def getLon(self, mjd):
        """Calculate the Sun's true longitude only, for the given MJD(s). """
        # Using the procedure outlined http://www.stjarnhimlen.se/comp/ppcomp.html
        # Calculate the days since J2000.
        days = mjd - _MJD2000
        # Calculate the obliquity of the ecliptic at this/these times.
        ecl = 23.4393 - 3.563E-7 * days
        # Calculate the 'orbital elements' of the Sun at this/these times. 
        N = 0.0
        i = 0.0
        w = (282.9404 + 4.70935E-5 * days) % 360.0
        a = 1.000000 
        e = 0.016709 - 1.151E-9 * days
        M = (356.0470 + 0.9856002585 * days) % 360.0
        # Calculate the eccentric anomaly for the sun.
        E = M + e*(_rad2deg)*numpy.sin(M*_deg2rad)*(1.0 + e * numpy.cos(M*_deg2rad))
        # Calculate the distance to the Sun and the true anomaly. 
        xv = numpy.cos(E*_deg2rad) - e
        yv = numpy.sqrt(1.0 - e*e) * numpy.sin(E*_deg2rad)
        v = numpy.arctan2( yv, xv ) * _rad2deg
        # Calculate the true longitude of the Sun.
        lon = v + w
        lon = lon % 360.0    
        return lon
    
    def calcPos(self, mjd):
        """Calculate the Sun's true longitude, distance, and RA/Dec, for the given MJD(s). """
        # Using the procedure outlined http://www.stjarnhimlen.se/comp/ppcomp.html
        self.mjd = numpy.copy(mjd)
        # Calculate the days since J2000.
        days = mjd - _MJD2000
        # Calculate the obliquity of the ecliptic at this/these times.
        ecl = 23.4393 - 3.563E-7 * days
        # Calculate the 'orbital elements' of the Sun at this/these times. 
        N = 0.0
        i = 0.0
        w = (282.9404 + 4.70935E-5 * days) % 360.0
        a = 1.000000 
        e = 0.016709 - 1.151E-9 * days
        M = (356.0470 + 0.9856002585 * days) % 360.0
        # Calculate the eccentric anomaly for the sun.
        E = M + e*(_rad2deg)*numpy.sin(M*_deg2rad)*(1.0 + e * numpy.cos(M*_deg2rad))
        # Calculate the distance to the Sun and the true anomaly. 
        xv = numpy.cos(E*_deg2rad) - e
        yv = numpy.sqrt(1.0 - e*e) * numpy.sin(E*_deg2rad)
        self.dist = numpy.sqrt(xv*xv + yv*yv)
        v = numpy.arctan2(yv, xv) * _rad2deg
        # Calculate the true longitude of the Sun.
        self.lon = v + w
        self.lon = self.lon % 360.0    
        # Now translate back to RA/Dec.
        # Calculate the sun's position in ecliptic x/y/z coordinates. 
        xs = self.dist * numpy.cos(self.lon*_deg2rad)
        ys = self.dist * numpy.sin(self.lon*_deg2rad)
        # zs = distance out of ecliptic plane, but this is always 0 for the Sun.
        # Calculate the sun's position in equatorial x/y/z geocentric coordinates. 
        xe = xs
        ye = ys * numpy.cos(ecl * _deg2rad)
        ze = ys * numpy.sin(ecl * _deg2rad)
        # And convert to RA/Dec coordinates (this is geocentric position). 
        self.ra  = numpy.arctan2(ye, xe) * _rad2deg
        self.dec = numpy.arctan2(ze, numpy.sqrt(xe*xe+ye*ye)) * _rad2deg
        return
    
    def getAltAz(self, skypos):
        """Get Alt/Az/Airmass information for each observation. Pass an already instantiated SkyPos object."""
        self.alt, self.az = skypos.radec2altaz(self.ra, self.dec, self.mjd)
        return

    
    def test_slalib(self, mjd, config):
        from pyslalib import slalib
        ra = numpy.zeros(len(mjd), 'float')
        dec = numpy.zeros(len(mjd), 'float')
        for i in range(len(mjd)):
            # Calculate sun's (topocentric) position. 
            #ra_RAD, dec_RAD, diam = slalib.sla_rdplan(mjd[i], 0,
            #                                          config['longitude']*_deg2rad, config['latitude']*_deg2rad)
            # Sun's (geocentric?) position.
            bary_vel, bary_pos, helio_vel, helio_pos = slalib.sla_evp(mjd, 2000)
            sun_pos = -1 * (bary_pos)
            ra_RAD, dec_RAD = slalib.sla_dcc2s(sun_pos)
            ra[i] = ra_RAD * _rad2deg
            dec[i] = dec_RAD * _rad2deg
        print "Sun: Test_slalib results"
        print " ... ra(slalib)"
        print ra
        print "  ... ra(self)"
        print self.ra
        print self.ra.min(), self.ra.max()
        print " ... dec(slalib)"    
        print dec
        print " ... dec(self)"
        print self.dec
        
