"""
Weather

A python class to hold site weather data for opsimObs.py.
This class reads the seeing and cloud data from text files.

The cloud and seeing files should have two pieces of information for each record:
  date (in seconds)   seeing or cloud value
  seeing == raw atmospheric seeing in arcseconds at 500nm at zenith.
  cloud == value between 0 to 1 representing amount of cloud extinction
In the default files (prepared from original text files SeeingPachon.txt and CloudTololo.txt),
  the seeing values are sampled every 300s for 2 years (data from Pachon)
  the cloud values are sampled four times per night and are reported in 1/8 increments by observers at CTIO

"""

import numpy
    
_sec2day = 1/60.0/60.0/24.0
_day2sec = 24.0*60.0*60.0

class Weather():
    """Weather: a class to setup and query in-memory sqlite databases of cloud and seeing values."""
    def __init__(self):
        return

    def readWeather(self, config):
        # Read text files from disk, store in memory. 
        # config is needed as it stores the names of the text files on disk.
        self.maxtime = {}
        self.dates = {}
        self.weather = {}
        for w in ('cloud', 'seeing'):
            self._setupWeather(w, config)
        return

    def getCloud(self, mjd, config):
        """Query the cloud data to return an interpolated cloud value. """
        w= 'cloud'
        cloud = self._checkWeather(mjd, w, config)
        return cloud
        
    def getSeeing(self, mjd, airmass, filter, config):
        """Query the seeing data to return an interpolated seeing value.
        Must provide the airmass and filter of observation to provide adjustments for the
        effect of the atmosphere and filter (and telescope components). """
        w = 'seeing'
        rawseeing = self._checkWeather(mjd, w, config)
        # If the raw seeing is 'too good to be true', swap to 'good_seeing_limit' instead.
        rawseeing = numpy.where(rawseeing<config['good_seeing_limit'], config['good_seeing_limit'], rawseeing)
        # Adjust seeing for airmass of observation.
        seeing = rawseeing * numpy.power(airmass, 0.6)
        # Adjust seeing for the filter.
        #  do this if the filters are a numpy array (and possibly different filters).
        if isinstance(filter, numpy.ndarray) | isinstance(filter, list):
            filterwavelength = numpy.empty(len(filter), float)
            for i in range(len(filter)):
                filterwavelength[i] = config['filter_wavelen'][filter[i]]
        else:
            # do this if there is only one filter given.
            filterwavelength = config['filter_wavelen'][filter] 
        basewavelength = config['seeing_wavelen']
        seeing = seeing * numpy.power(basewavelength / filterwavelength, 0.2)
        # Adjust seeing for the telescope systematic floor.        
        seeing = numpy.sqrt(seeing**2 + config['seeing_Telescope']**2)
        return seeing
    
    def _setupWeather(self, w, config):
        """Read a data file into  numpy arrays. (w == 'seeing' or 'cloud'). """
        wnames = ('cloud', 'seeing')
        if w not in wnames:
            raise Exception('w should be one of %s' %(wnames))
        filename = config['%s_datafile' %(w)]
        file = open(filename, 'r')
        # Also assume flat file contains only date / value in a space or tab separated file. 
        self.dates[w] = []
        self.weather[w] = []
        # Read the data file.
        print '# Reading weather data file %s' %(filename)
        for line in file:
            if line.startswith('#') | line.startswith('!'):
                continue
            self.dates[w].append(line.split()[0])
            self.weather[w].append(line.split()[1])
        file.close()
        self.dates[w] = numpy.array(self.dates[w], float)
        self.weather[w] = numpy.array(self.weather[w], float)
        # Check the total amount of data (mostly for user awareness):
        print '# Read %d weather values from %s file. ' %(len(self.weather[w]), filename)
        # Get the total length of time included in this (seeing/cloud) file,
        #  so that we can determine a wrap-around date if we need that.
        self.maxtime[w] = self.dates[w].max()
        return 

    def _checkWeather(self, mjd, w, config):
        """Check the weather (clouds or seeing) for a given mjd, either single value or multiples."""
        if isinstance(mjd, float):
            return self._do_checkWeather(mjd, db, config)
        if (isinstance(mjd, float) == False):
            # So let's assume mjd was either a list or a numpy array. Just iterate through it. 
            values = numpy.empty(len(mjd), float)
            for i in range(len(mjd)):
                values[i] = self._do_checkWeather(mjd[i], w, config)
            return values


    def _do_checkWeather(self, mjd, w, config):
        """Check in-memory cloud or seeing information for a single MJD.
        Returns an interpolated or extrapolated value (from nearest neighbors in time).
        Since weather databases may cover a shorter length of time than the simulated survey,
        the queries are 'wrapped' in time into the weather database information
        (i.e. if there are 2 years of seeing values, a survey observation at 2.5 years would
        be wrapped to return the seeing at 0.5 years). """
        # Convert mjd to the relevant time units of the weather dates.
        time = (mjd - config['sim_start'] + config['%s_start' %(w)]) * _day2sec
        # And wrap the time, if we need to. 
        time = time % self.maxtime[w]
        # Find the observations which are closest in time to our requested time.
        time_order = (abs(self.dates[w] - time)).argsort()
        date1 = self.dates[w][time_order[0]]
        date2 = self.dates[w][time_order[1]]
        weather1 = self.weather[w][time_order[0]]
        weather2 = self.weather[w][time_order[1]]
        # Do interpolation for weather at this particular time.
        weather = (weather2 - weather1) / (date2 - date1) * (time - date1) + weather1
        # Check within 0-1 limits.
        if weather<0:
            weather = 0
        if weather>1:
            weather = 1
        return weather




        
