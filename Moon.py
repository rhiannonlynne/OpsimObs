"""
Moon

Class to calculate the approximate position and phase of the Moon. 
"""

import numpy
from pyslalib import slalib
from Sun import Sun

_deg2rad = numpy.pi/180.0
_rad2deg = 180.0/numpy.pi
_MJD2000 = 51544.5


class Moon():

    def __init__(self):
        """Instantiate the Sun object."""
        return

    def calcPos(self, mjd, config):
        """Calculate the moon's ecliptic lon/lat and geocentric RA/Dec, for the given MJD(s)."""
        self.mjd = numpy.copy(mjd)
        self.ra = numpy.zeros(len(mjd), 'float')
        self.dec = numpy.zeros(len(mjd), 'float')
        for i in range(len(mjd)):
            # Calculate moon's position. 
            ra_RAD, dec_RAD, diam = slalib.sla_rdplan(mjd[i], 3,
                                                      config['longitude']*_deg2rad, config['latitude']*_deg2rad)             
            self.ra[i] = ra_RAD * _rad2deg
            self.dec[i] = dec_RAD * _rad2deg
        # Calculate the lunar phase. 
        eclon = numpy.zeros(len(mjd), 'float')
        eclat = numpy.zeros(len(mjd), 'float')
        for i in range(len(mjd)):
            eclon[i], eclat[i] = slalib.sla_eqecl(self.ra[i]*_deg2rad, self.dec[i]*_deg2rad, self.mjd[i])            
        # Calculate the solar longitude. 
        sun = Sun()
        lon_sun = sun.getLon(self.mjd)*_deg2rad
        # Calculate the solar elongation of the Moon.
        solarelong = numpy.arccos(numpy.cos((lon_sun - eclon)) * numpy.cos(eclat))
        # Calculate the phase of the moon.
        self.phase = 180.0 - solarelong*_rad2deg
        # Calculate the illumination of the Moon. 
        self.illum = ( 1 + numpy.cos(self.phase*_deg2rad) )/2.0
        return
    
    def getAltAz(self, skypos):
        """Get Alt/Az/Airmass information for each observation. Pass an already instantiated SkyPos object."""
        # Calculate the alt/az for these times, 
        self.alt, self.az = skypos.radec2altaz(self.ra, self.dec, self.mjd)
        return

