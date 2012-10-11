"""
Weather

### don't use this .. it's much slower than using numpy arrays to store and return this information.

A python class to hold site weather data for opsimObs.py.
This class reads the seeing and cloud data from on-disk sqlite or text files, and places this
information into an in-memory sqlite database. This database is queried to generate the
cloud and seeing information required for each pointing. 

The cloud and seeing database files should have two pieces of information for each record:
  date (in seconds)   seeing or cloud value
  seeing == raw atmospheric seeing in arcseconds at 500nm at zenith.
  cloud == value between 0 to 1 representing amount of cloud extinction
In the default files (prepared from original text files SeeingPachon.txt and CloudTololo.txt),
  the seeing values are sampled every 300s for 2 years (data from Pachon)
  the cloud values are sampled four times per night and are reported in 1/8 increments by observers at CTIO

"""

import numpy
import sqlite3

_sec2day = 1/60.0/60.0/24.0
_day2sec = 24.0*60.0*60.0

class Weather():
    """Weather: a class to setup and query in-memory sqlite databases of cloud and seeing values."""
    def __init__(self):
        return

    def readWeather(self, config):
        """Instantiate object and set up databases in memory. """
        # config is needed as it stores the names of the database dump files on disk.
        self.maxtime = {}
        self.conn = {}
        self.cursor = {}
        # Set up the cloud database, then set up the seeing database.
        # Get maxtime values (required for wraparound of queries in time if db length < survey length)
        #  and get conn and cursors for database.
        for db in ('cloud', 'seeing'):
            self.maxtime[db], self.conn[db], self.cursor[db] = self._setupWeatherDB(db, config)
        return

    def getCloud(self, mjd, config):
        """Query the sqlite in-memory database to return an interpolated cloud value. """
        db= 'cloud'
        cloud = self._queryDB(mjd, db, config)
        return cloud
        
    def getSeeing(self, mjd, airmass, filter, config):
        """Query the sqlite in-memory database to return an interpolated seeing value.
        Must provide the airmass and filter of observation to provide adjustments for the
        effect of the atmosphere and filter (and telescope components). """
        db = 'seeing'
        rawseeing = self._queryDB(mjd, db, config)
        # Adjust seeing for airmass of observation.
        seeing = rawseeing * numpy.power(airmass, 0.6)
        # Adjust seeing for the filter.
        if isinstance(filter, numpy.ndarray) | isinstance(filter, list):
            filterwavelength = numpy.empty(len(filter), float)
            for i in range(len(filter)):
                filterwavelength[i] = config['filter_wavelen'][filter[i]]
        else:
            filterwavelength = config['filter_wavelen'][filter] 
        basewavelength = config['seeing_wavelen']
        seeing = seeing * numpy.power(basewavelength / filterwavelength, 0.2)
        # Adjust seeing for the telescope systematic floor.
        seeing = numpy.sqrt(seeing**2 + config['seeing_Telescope']**2)
        return seeing
    
    def _setupWeatherDB(self, db, config):
        """Read a data file into a sqlite in-memory database. (dbname == 'seeing' or 'cloud'). """
        dbnames = ('seeing', 'cloud')
        if db not in dbnames:
            raise Exception('dbname should be one of %s' %(dbnames))
        filename = config['%s_datafile' %(db)]
        print '# Will read %s data from %s' %(db, filename)
        # If filename ends with .sqlite, assume it holds data from a sqlite file.
        if filename.endswith('.sqlite'):
            print '# Assuming this file is full sqlite dump. '
            file = open(filename, 'r')
            # Read the sqlite file. 
            with file:
                data = file.read()
            # Set up the sqlite in memory database.
            conn = sqlite3.connect(':memory:')
            cursor = conn.cursor()
            # Do the data insert. 
            cursor.executescript(data)
            file.close()
        # If the filename ends with .txt, assume it is a flat data file and must be inserted manually.
        # Also assume flat file contains only date / value in a CSV file. 
        else:
            print '# Assuming this file is text or data file only, and should insert each individual line of data.'
            file = open(filename, 'r')
            # Set up the in memory database
            conn = sqlite3.connect(':memory:')
            cursor = conn.cursor()
            cursor.execute('create table %s (date float, %s float)' %(db, db))
            # Read the data file and insert into database. 
            for line in file:
                if line.startswith('#') | line.startswith('!'):
                    continue
                try:
                    (date, val) = line.strip().split(',')                    
                    insertquery  = 'insert into %s(date, %s) values(?,?)' %(db, db)
                    cursor.execute(insertquery, (date, val))
                except:
                    continue
            file.close()
        # Check the total amount of data (mostly for user awareness):
        cursor.execute('select count(date) from %s' %(db))
        print '# Added %d reference observations to %s table.' %(cursor.fetchall()[0][0], db)
        # Add an index on time. 
        cursor.execute('create index date_idx on %s(date)' %(db))
        print '# Added an index on date' 
        # Get the total length of time included in this (seeing/cloud) database,
        #  so that we can determine a wrap-around date for the table if we need that.
        query = 'select max(date) from %s' %(db)
        cursor.execute(query)
        res = cursor.fetchall()
        print res
        date_max = float(res[0][0])
        return date_max, conn, cursor

    def _queryDB(self, mjd, db, config):
        """Query in-memory sqlite database for cloud or seeing information.
        Returns an interpolated or extrapolated value (from nearest neighbors in time).
        Since weather databases may cover a shorter length of time than the simulated survey,
        the queries are 'wrapped' in time into the weather database information
        (i.e. if there are 2 years of seeing values, a survey observation at 2.5 years would
        be wrapped to return the seeing at 0.5 years). """
        # MJD should be single value.
        if isinstance(mjd, float):
            return self.do_queryDB(mjd, db, config)
        if (isinstance(mjd, float) == False):
            # So let's assume mjd was either a list or a numpy array. Just iterate through it. 
            values = numpy.empty(len(mjd), float)
            for i in range(len(mjd)):
                values[i] = self.do_queryDB(mjd[i], db, config)
            return values

    def do_queryDB(self, mjd, db, config):
        """Do the actual (single MJD) query for queryDB."""
        # Adjust time of db query for start of simulated survey and then for potential offset into weather db.
        qtime = (mjd - config['sim_start'] + config['%s_start' %(db)])* _day2sec
        # Then adjust time of db query for length of database (wrap the time). 
        qtime = qtime % self.maxtime[db]
        # Return the two observations which are closest in time to our requested time of db query.
        query = 'select date, %s from %s order by abs(date-%f) limit 2' %(db, db, qtime)
        self.cursor[db].execute(query)
        results = self.cursor[db].fetchall()
        # Assign these results to numpy arrays. 
        dates = numpy.zeros(len(results))
        dbvals = numpy.zeros(len(results))
        for i in range(len(results)):
            dates[i] = float(results[i][0])
            dbvals[i] = float(results[i][1])
        # Do the interpolation.
        value = (dbvals[1] - dbvals[0]) / (dates[1] - dates[0]) * (qtime - dates[0]) + dbvals[0]
        if value > 1.0:
            value = 1.0
        if value < 0.0:
            value = 0.0
        return value

        
