"""
ObsFields

A python class to hold the information about each observation for opsimSky.py.
This class reads the input file and does some basic checks on whether the timing of each observation is appropriate.
It can update each observation with information about 'where' in the sky (alt/az/airmass) it would be taken,
and retrieve the cloud and seeing information for each observation.
"""

import numpy

_sec2day = 1/60.0/60.0/24.0
_day2sec = 24.0*60.0*60.0

class ObsFields():
    """A class to hold information about each observation for opsimSky."""
    def __init__(self):
        return

    def clearVals(self):
        self.ra = None
        self.dec =  None
        self.mjd = None
        self.filter = None
        self.exptime = None
        self.alt = None
        self.az = None

    def readInputObs(self,inputfilename):
        """Read input values - RA Dec Date(UTC-MJD) Filter Exptime
        Sets RA/Dec/MJD/Filter/exptime for class. """
        # Read data from file.
        infile = open(inputfilename, 'r')
        self.ra = []
        self.dec = []
        self.mjd = []
        self.filter=  []
        self.exptime = []
        print '# Reading observation data from file %s' %(inputfilename)
        for line in infile:
            if line.startswith('#') | line.startswith('!'):
                continue
            try:
                (ra, dec, mjd, filter, exptime) = line.split()
                self.ra.append(ra)
                self.dec.append(dec)
                self.mjd.append(mjd)
                self.filter.append(filter)
                self.exptime.append(exptime)
            except:
                continue
        infile.close()
        # Convert to numpy arrays. 
        self.ra = numpy.array(self.ra, float)
        self.dec = numpy.array(self.dec, float)
        self.mjd = numpy.array(self.mjd, float)
        self.filter = numpy.array(self.filter, str)
        self.exptime = numpy.array(self.exptime, float)
        print '# Read %d input observations.' %(len(self.ra))
        return 

    def _calcNight(self, midnight=0.16):
        """Given a time (in MJD/UTC), and the time at midnight at the site (in lsstsite),
        calculate and return the truncated 'day'. """
        # Note that this could be applied to an array of mjd's if needed.
        self.night = numpy.floor(self.mjd - midnight - 0.5)
        return     

    def checkInputObs(self, config):
        """Do some simple consistency checks to see if observations could be scheduled as given.
        Will raise exception if checks fail. """
        # Check data.
        # Do the filters match the filters we know about?
        print "# Assuming %d exposures per visit/line in input file" %(config['nexp_visit'])
        filters = numpy.unique(self.filter)
        for f in filters:
            if f not in config['filter_names']:
                raise Exception('Unknown filter in input observation file: %s' %(f))
        # Do the observations need more than config['n_simultaneous_filters'] filters in a night?
        self._calcNight(config['midnight'])
        for n in numpy.unique(self.night):
            condition = (self.night == n)
            n_filters_night = len(numpy.unique(self.filter[condition]))
            if n_filters_night > config['n_simultaneous_filters']:
                raise Exception('Cannot request more than %d filters within a night (violated in night %.0f requesting %d filters)' % (config['n_simultaneous_filters'], n, n_filters_night))
        # Is the spacing between observations sufficient?
        self.opentime = (self.exptime[0:len(self.exptime)-1]) \
            + config['nexp_visit']* config['add_shutter']*config['shutter_time'] \
            + (config['nexp_visit'] - 1 * ((1+config['add_slew'])%2)* config['readout_time']) \
            + config['add_slew'] * config['slew_time']
        deltatime = numpy.diff(self.mjd) * _day2sec
        diff = deltatime - self.opentime
        thresh = -2
        if diff.min()<thresh:
            # This thresh should be 0, but there are rounding errors which seem associated with
            #   getting information from opsim DB & using MJD that mean the seconds can be off by a bit. 
            condition = (diff<thresh)
            raise Exception('Not enough time between observations near %s to include shutter and readout time' %(self.mjd[condition]))
        # Is the filter time change sufficient?
        f_prev = self.filter[0]
        for i in range(1, len(self.filter)):
            if self.filter[i] != f_prev:
                f_prev = self.filter[i]
                if ((self.mjd[i] - self.mjd[i-1])*_day2sec) < config['filter_time']:
                    raise Exception('Cannot request observations with filter changes spaced less than %f seconds apart. Fails at %f'  %(config['filter_time'], self.mjd[i-1]))
        # Finished basic checks.    
        print '# Passed basic checks on observations.'
        return


    def getAltAzAirmass(self, skypos):
        """Get Alt/Az/Airmass information for each observation. Pass an already instantiated SkyPos object."""
        self.alt, self.az = skypos.radec2altaz(self.ra, self.dec, self.mjd)
        self.airmass = skypos.alt2airmass(self.alt)
        return

        
    def getObsWeather(self, weather, config):
        """Get weather information (clouds/seeing) at observation time for each field.
        Pass an already-initiated 'weather' object (i.e. has data already in memory). """
        # Get cloud information. 
        self.cloud = weather.getCloud(self.mjd, config)
        # Get seeing information. 
        self.seeing = weather.getSeeing(self.mjd, self.airmass, self.filter, config)
        return
    
    def getDowntime(self, downtime, config):
        """Get downtime status at observation time for each field.
        Pass an already-initiated 'downtime' object (i.e. with data already in memory). """
        # Get downtime status.
        self.downstatus = downtime.checkDownstatus(self.mjd, config)
        return


    def printObs(self, sun, moon, sky, maglimit, config):
        survey_day = numpy.floor(self.mjd - config['sim_start'] - config['midnight'] + 0.5)
        print '# RA  Dec  MJD  filter  Alt  Az  Airmass  Cloud   Seeing  SunAlt SunAz   MoonAlt  MoonAz  MoonPhase MoonPhase(opsim) SkyBrightness Maglimit SurveyNight Downtime?(0=No,1=Yes)'
        for i in range(len(self.mjd)):
            print '%.5f %.5f %.7f %s %.5f %.4f %.3f %.2f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %d' \
                % (self.ra[i], self.dec[i], self.mjd[i], self.filter[i],
                   self.alt[i], self.az[i], self.airmass[i],
                   self.cloud[i], self.seeing[i],
                    sun.alt[i], sun.az[i], moon.alt[i], moon.az[i], moon.phase[i], moon.illum[i], sky[i], maglimit[i], survey_day[i], self.downstatus[i])
        return
