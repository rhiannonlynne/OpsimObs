"""
SkyBright - 

Original code written by Perry Gee. Reformulated for multiple field pointings by Lynne Jones.
Questions: ljones@astro.washington.edu

Class data: 
model                  Use Perry or K&S model for sky brightness
alt_field, az_field    position of the field pointing (degrees)
alt_moon, az_moon      position of the moon (degrees)
moon_phase             phase angle of the moon, from -180 (new) to 0 (full) to 180 (new).
                       Note: phase angle = arccos(moon illumination / 50 - 1)
alt_sun, az_sun        position of the sun  (not currently used)
kAtm                   extinction coefficients (0-1)
bandpass               bandpass of  ('u', 'g', 'r', 'i', 'z', 'y3', 'y4', 'y')

All angles are in DEGREES.
The methods could all handle arrays instead of single values, except that only a 
single bandpass at a time is supported. 

Note that the color terms (extinction and coefficients for lunar color change with phase) for
  'y' are equal for both 'y3' & 'y4' at the moment, since I don't have the values for 'y4'. 

Typical usage for this class would be as follows:
 s = SkyBright(model='Perry', solar_phase='ave')
 s.setSkyBright(field_alt, field_az, moon_alt, moon_az, moon_phase, bandpass='g')
 skybrightness = s.getSkyBright()

However, you could repeat runs with the same field position and change the bandpass (s.setBandpassAndK(bandpass, k))
or retrieve the moon and dark sky  skybrightness values separately (s.getMoonSky(), s.getDarkSky()).
Or swap between models (KS vs 'Perry') or solar phases. 

"""

import warnings
import numpy

rad2deg = 180.0/numpy.pi
deg2rad = numpy.pi/180.0

##
# Parameters for Perry's sky model (results of empirical fits)
#
bandpass_types = ['u', 'g', 'r', 'i', 'z', 'y', 'y3', 'y4']
# Extinction coefficient dictionary, used if no other value given
# Based on average values taken from CTIO from 2000-2005 (from Perry Gee) 
k_Extinction = {'u': 0.4837, 'g': 0.1756, 'r': 0.09398, 'i': 0.05958, 'z': 0.05217, 'y': 0.03, 'y3':0.03, 'y4':0.04}
#
# Coefficients for change in lunar brightness as a function of phase
# Based on data from Southern Standard Stars survey from CTIO 0.9m
lunarp0 = {'u':1.8819e+01, 'g':1.8301e+01, 'r':1.8938e+01, 
           'i':1.9032e+01, 'z':1.8547e+01, 'y':1.8108e+01, 'y3':1.8108e+01, 'y4':1.8108e+01}
lunarp1 = {'u':3.5802e-02, 'g':3.6370e-02, 'r':2.6097e-02, 
           'i':1.3086e-02, 'z':4.4932e-03, 'y':1.2908e-04, 'y3':1.2908e-04, 'y4':1.2908e-04}
lunarp2 = {'u':-7.6183e-05,'g':-9.0022e-05,'r':-8.0058e-05,
           'i':-4.1280e-05,'z':-1.3532e-05, 'y':2.9816e-16, 'y3':2.9816e-16, 'y4':2.9816e-16}
# Coefficients to convert lunar magnitudes to LSST bandpasses
lunarToLSSTMag = {'u':-0.0621,'g':-0.0039, 'r':0.0004, 'i':0.0006, 'z':-0.0096, 'y':-0.0189, 'y3':-0.0189, 'y4':-0.0189}
#
# Values of dark sky backgrounds from SDSS/DLS/Standard Southern Stars/UKIDDS bandpasses
# for solar minimum
dark_min = {'u':22.95, 'g':22.32, 'r':21.43, 'i':20.40, 'z':19.22, 'y':18.11, 'y3':18.11, 'y4':18.11}
# for solar average
dark_ave = {'u':22.80, 'g':22.20, 'r':21.30, 'i':20.30, 'z':19.10, 'y':17.96, 'y3':17.96, 'y4':17.96}
# for solar maximum
dark_max = {'u':22.65, 'g':22.10, 'r':21.13, 'i':20.12, 'z':18.92, 'y':17.81, 'y3':17.81, 'y4':17.81}
# Coefficients to convert dark sky magnitudes to LSST bandpasses
darkToLSSTMag = {'u':0.1757, 'g':-0.3135, 'r':-0.0081,'i':-0.0445, 'z':0.5726, 'y':0.2735, 'y3':0.2735, 'y4':0.2735}
#
# Coefficients to calculate change in lunar brightness as field/moon separation changes
scatp0 = {'u':2.5573e-01, 'g':3.9323e-01, 'r':2.4854e-01, 
          'i':2.3688e-01, 'z':3.1613e-01, 'y':7.5894e-17, 'y3':7.5894e-17, 'y4':7.5894e-17}
scatp1 = {'u':3.1554e-13, 'g':1.7097e+00, 'r':2.9850e+00, 
          'i':4.1227e+00, 'z':5.5089e+00, 'y':2.7105e-18, 'y3':2.7105e-18, 'y4':2.7105e-18}


class SkyBright:

    def __init__(self, model='Perry', solar_phase='ave'):
        """ Initialize sky brightness class object.
        Set type of model and solar phase (although these can be reset later)."""
        self.setModelType(model)
        self.setSolarPhase(solar_phase)
        return

    def setModelType(self, model):
        # Set model type.
        model_types = ['KS', 'Perry']
        if model not in model_types:
            raise ValueError("Model must be one of %s" %(model_types))
        self.model = model        
        return

    def setSolarPhase(self, solar_phase):
        # Set solar phase.
        solar_phase_types = ['min', 'ave', 'max']
        if solar_phase not in solar_phase_types:
            raise ValueError("Solar_Phase must be one of %s" %(solar_phase_types))
        self.solar_phase = solar_phase
        return

    def setBandpassAndK(self, bandpass, k=None):
        # Set bandpass.
        if bandpass not in bandpass_types:
            raise ValueError("Bandpass must be one of %s" %(bandpass_types))
        self.bandpass = bandpass
        if self.model == "KS":
            self.bandpass = "V"  # not one of the types, but dictionaries not called.
        # Set extinction coefficient.
        if k == None:
            self.k=self._setExtinction()
        else:
            self.k = k
        return

    def _setExtinction(self):
        # Set default values for extinction coefficients. 
        if self.model == "KS":
            k = 0.113
        elif self.model == "Perry":
            k = k_Extinction[self.bandpass]
        return k
        
    
    def setSkyBright(self, field_alt, field_az, moon_alt, moon_az, moon_phase,
                     bandpass='g', k=None, model=None, solar_phase=None):
        """Set parameters (field pointings, moon and sun position) that need
        before calculating sky brightness.
        All angle values are given in degrees."""
        # Have an opportunity to reset model and solar_phase, if desired. 
        if model != None:
            self.setModelType(model)
        if solar_phase != None:
            self.setSolarPhase(solar_phase)
        # Set bandpass and extinction coefficient. 
        self.setBandpassAndK(bandpass, k)
        # Set other values. Remember all angles in DEGREES.
        self.moon_phase = moon_phase
        self.moon_alt = moon_alt
        self.moon_az = moon_az
        if numpy.shape(self.moon_alt) != numpy.shape(self.moon_az):
            raise ValueError("Moon altitude and azimuth must have the same length.")
        self.field_alt = field_alt
        self.field_az = field_az
        if numpy.shape(self.field_alt) != numpy.shape(self.field_az):
            raise ValueError("Field altitude and azimuth must have the same length.")        
        # Check that incoming values are compatible .. 
        #  that is, if field values are array - then moon alt/az has to be array of same length, or single value
        if isinstance(self.moon_alt, numpy.ndarray):
            if isinstance(self.field_alt, numpy.ndarray):
                if (numpy.shape(self.moon_alt) != numpy.shape(self.field_alt)):
                    raise ValueError('If both field and moon positions are arrays, they must be the same length.')
        if isinstance(self.moon_phase, numpy.ndarray):
            if (numpy.shape(self.moon_alt) != numpy.shape(self.moon_phase)):
                raise ValueError('If the moon phase is given as an array, it must be the same length as the moon positions.')
        return

    def printVals(self):
        """Print values currently being used in sky brightness calculations / class."""
        print "Moon_phase  moon_alt  moon_az   field_alt  field_az  model  bandpass k  solar_phase"
        print self.moon_phase, self.moon_alt, self.moon_az, self.field_alt, self.field_az, \
            self.model, self.bandpass, self.k, self.solar_phase
        print "Total Sky Brightness"
        print self.getSkyBright()
        return

    def _calcAngle(self, altitude1, azimuth1, altitude2, azimuth2):
        """Calculates distance on a sphere using the Vincenty formula. 
        Give this function values in DEGREES. Returns angular distance(s), in DEGREES."""
        # This one is supposed to be accurate for all distances, but is more complex. 
        alt1 = altitude1 * deg2rad
        az1 = azimuth1 * deg2rad
        alt2 = altitude2 * deg2rad
        az2 = azimuth2 * deg2rad
        D1 = (numpy.cos(alt2)*numpy.sin(az2-az1))**2 + \
            (numpy.cos(alt1)*numpy.sin(alt2) - numpy.sin(alt1)*numpy.cos(alt2)*numpy.cos(az2-az1))**2
        D1 = numpy.sqrt(D1)
        D2 = (numpy.sin(alt1)*numpy.sin(alt2) + numpy.cos(alt1)*numpy.cos(alt2)*numpy.cos(az2-az1))
        D = numpy.arctan2(D1,D2) 
        D = D * rad2deg
        return D

    def _calcAngle2(self, altitude1, azimuth1, altitude2, azimuth2):
        """Calculate the angle between two alt/az coordinates. Input values in DEGREES."""
        # Tested this with comparison to above version, and give similar results.
        # This is what Perry was using in his ETC.py code.
        alt1 = altitude1 * deg2rad;
        az1 = azimuth1 * deg2rad;
        alt2 = altitude2 *  deg2rad;
        az2 = azimuth2 * deg2rad;
        dx = numpy.cos(az1)*numpy.cos(alt1) - numpy.cos(az2)*numpy.cos(alt2);
        dy = numpy.sin(az1)*numpy.cos(alt1) - numpy.sin(az2)*numpy.cos(alt2);
        dz = numpy.sin(alt1) - numpy.sin(alt2);
        d = numpy.sqrt(dx*dx + dy*dy + dz*dz)/2.0;
        return 2 * (numpy.arctan2(d, numpy.sqrt(1-(d**2)))*rad2deg)

    def _calcAirmass(self, altitude):
        """Calculate the airmass for fields with low-airmass or 
        for when estimating the scattered light from the moon.
        Zenith_distance in DEGREES.  Good up to ~2 airmasss."""
        z = (90 - altitude) * deg2rad
        X = numpy.power((1 - 0.96*(numpy.sin(z))**2), -0.5)
        return X
    
    def _calcAirmassLargeDirect(self, altitude):
        """Calculate the airmass appropriate for direct observed 
        when the field is close to the horizon.
        Zenith_distance in DEGREES.
        Note that K&S show this is *not correct* for use in estimating
        scattered light from the moon, even when the moon is low on horizon."""
        z = (90.0 - altitude) * deg2rad
        X = (numpy.cos(z) + 0.025*numpy.exp(-11*numpy.cos(z)))**-1
        return X
    
    def _calcLunarScattering(self, angsep):
        """Calculate the change in lunar brightness as the distance between the field
        and the moon changes.  angsep in DEGREES."""
        if self.model == "KS":       
            # Formula from K&S. 
            f_rho = (10**(5.36))*(1.06 + numpy.cos(angsep*deg2rad)**2) + \
                10**(6.15-angsep/40.0)
        elif self.model == "Perry":   
            # Perry's model is based on the K&S formula, but uses coefficients fit 
            # from data in each bandpass.
            f_rho = scatp0[self.bandpass] * (1.06 + numpy.cos(angsep*deg2rad)**2) + \
                scatp1[self.bandpass] * (10**(-angsep/40.0)) 
        return f_rho

    def _calcLunarScatteringIntoFOV(self, airmass, k=None):
        """Calculate change in sky brightness due to change in field airmass. 
        Higher airmass / higher k means more air volume for scattering, thus higher brightness."""
        # fraction of scattered moonlight that reaches the observer
        if k == None:
            k = self.k
        # this airmass should be the airmass of the field.
        scattering = 1 - 10**(-0.4*k*airmass) 
        return scattering

    def _calcAirmassExtinction(self, airmass, k=None):
        """Calculate extinction along line of sight.
        Higher airmass / higher k means objects above atmosphere are reduced in brightness more."""
        if k==None:
            k = self.k
        # this airmass should be the airmass of the moon, calculated using normal airmass func.
        extinction = 10**(-0.4 * k * airmass)
        return extinction
        

    def _calcMoonIllumAboveAtmo(self):
        """Calculate the moon illuminance above the atmosphere."""
        # Notice that the units for K&S and for Perry are NOT the same.
        # KS = nanoLamberts, Perry = ...
        if self.model == "KS":
            # Basic K&S value (no bandpass correction)
            alpha = numpy.abs(self.moon_phase) 
            mm = 3.84 + 0.026*alpha + 4e-9*alpha**4
            I_o = 10**(-0.4*mm)
            if isinstance(self.moon_phase, float):
                if alpha  < 7:
                    I_o = I_o + (I_o * 0.35 * (7-alpha)/7.0)
                #raise warnings.warn("Moon phase < 7 degrees; adding estimate for opposition effect")
            elif isinstance(self.moon_phase, numpy.ndarray):
                condition = (alpha < 7.0)
                I_o[condition] = I_o[condition] + I_o[condition]*0.35*(7.0-alpha[condition])/7.0
        elif self.model == "Perry":
            # Perry Gee model (from ETC.py)
            alpha = numpy.abs(self.moon_phase)
            # Calculate lunar brightness at the current phase angle.
            mmp = (lunarp0[self.bandpass] + alpha*lunarp1[self.bandpass]
                   + (alpha**2)*lunarp2[self.bandpass])
            # Calculate the lunar brightness at new moon. 
            mmo = (lunarp0[self.bandpass] + 180.0*lunarp1[self.bandpass] 
                   + (180.0**2)*lunarp2[self.bandpass])
            # Calculate the brightness of the moon as the difference between the above values.
            I_s = 10**(-0.4*mmp) - 10**(-0.4*mmo)
            if isinstance(I_s, float):
                if I_s <= 0:
                    I_s = 1e-100
            else:
                condition = (I_s <= 0)
                I_s[condition] = 1e-100
            # Add lunar 'color' term
            I_s = I_s* 10**(-0.4*lunarToLSSTMag[self.bandpass])
            # Perry's model above is actually for a set of 'STANDARD' values
            # moon-field separation = 90 deg, field airmass = 1, 
            # moon zenith angle = 50  (altitude = 40)
            # So if we want to convert this intensity to a 'top of the atmosphere' value we must
            # 'undo' these effects
            # First:  moon-field separation
            I_s = I_s / self._calcLunarScattering(90)
            # Next: undo the extinction of the moonlight through the atmosphere that Perry 
            #  originally used! (note, this is not what should be used .. different airmass value)
            x = self._calcAirmassLargeDirect(40.0)
            tmp = 10**(-0.4*k_Extinction[self.bandpass]*x)
            I_s = I_s / tmp
            # Finally: undo the addition of scattering into the volume of air along the line of sight
            I_o = I_s / self._calcLunarScatteringIntoFOV(1, k_Extinction[self.bandpass])
        return I_o
            
    def getMoonSky(self, returnInMag=False):
        """Calculate the moon brightness."""        
        # Calculate the moon brightness above the atmosphere.
        I_o = self._calcMoonIllumAboveAtmo()
        # Calculate the separation between the field and the moon.
        rho = self._calcAngle2(self.moon_alt, self.moon_az, self.field_alt, self.field_az)
        # Include lunar scattering into field of view dependent on angle between field/moon.
        I_s = I_o * self._calcLunarScattering(rho)
        # Include extinction of lunar light.
        I_s = I_s * self._calcAirmassExtinction(self._calcAirmass(self.moon_alt), self.k)
        # Include lunar scattering into line of sight.
        I_s = I_s * self._calcLunarScatteringIntoFOV(self._calcAirmass(self.field_alt), self.k)
        if returnInMag:
            if self.model=="KS":
                return ((20.7233 - numpy.log(I_s/34.08))/0.92104)
            if self.model=="Perry":
                return (-2.5*numpy.log10(I_s))
        return I_s        
                   
    def getDarkSky(self, returnInMag=False):
        """Calculate dark sky brightness."""
        if self.model=="KS":
            darkmag = 21.8
            I_o = 34.08 * numpy.exp(20.7233 - 0.92104*darkmag)
            # note that Perry and KS use different units.
        elif self.model=="Perry":
            if self.solar_phase == "min":
                darkmag = dark_min[self.bandpass]
            elif self.solar_phase == "ave":
                darkmag = dark_ave[self.bandpass]
            elif self.solar_phase == "max":
                darkmag = dark_max[self.bandpass]
            darkmag = darkmag + darkToLSSTMag[self.bandpass]
            I_o = 10**(-0.4*darkmag)
        # Add extinction for dark sky proportional to 10^(-0.4k(X-1)),
        #  and additional emission proportional to X
        # Note this does not use airmass directly, but rather X-1
        airmass = self._calcAirmass(self.field_alt)
        I_s = I_o * airmass * self._calcAirmassExtinction((airmass-1), self.k)
        if returnInMag:
            if self.model=="KS":
                return ((20.7233 - numpy.log(I_s/34.08))/0.92104)
            if self.model=="Perry":
                return (-2.5*numpy.log10(I_s))
        else:
            return I_s

    def getSkyBright(self):
        """Calculate the total sky brightness."""
        # Find the dark sky value.
        darkI = self.getDarkSky()
        # Find the moon value.
        # First check if the moon is above the horizon.
        if (isinstance(self.moon_alt, float)) | (isinstance(self.moon_alt, int)):
            if self.moon_alt <= 0:
                moonI = 0
            else:
                moonI = self.getMoonSky()
        else: # assume otherwise dealing with numpy array
            moonI = self.getMoonSky()
            condition = (self.moon_alt <= 0 )
            moonI[condition] = 0
        # Add and convert to magnitudes. 
        if self.model == "KS":
            skymag = (20.7233 - numpy.log((moonI + darkI)/34.08))/0.92104
        elif self.model=="Perry":
            skymag = -2.5 * numpy.log10(darkI + moonI)
        return skymag

    
