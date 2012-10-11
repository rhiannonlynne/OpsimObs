"""
Downtime

A class to read the opsim downtime configuration files and return yes/no answer to whether telescope is in
downtime at time of observations.

"""

import numpy


class Downtime():
    """Class to read downtime information and evaluate whether a particular time is/is not during Downtime."""
    def __init__(self):
        return

    def readDowntime(self, config):
        """Read the scheduled downtime file. """
        self.downdates = {}
        for s in ('schedDowntime', 'unschedDowntime'):
            self._setupDowntime(s, config)
        self.alldowndates = numpy.concatenate([self.downdates['schedDowntime'],
                                               self.downdates['unschedDowntime']])
        return
    
    def _setupDowntime(self, s, config):
        """Read the downtime data files and build an array of ."""
        snames = ('schedDowntime', 'unschedDowntime')
        if s not in snames:
            raise Exception('t should be one of %s' %(snames))
        filename = config['%s_datafile' %(s)]
        file = open(filename, 'r')
        # Read the data file.
        # Assume that downtime conf files are formatted as a series of
        #  entries, giving startNight & duration for each downtime period.
        print '# Reading downtime data file %s' %(filename)
        startdates = []
        durations = []
        for line in file:
            if line.startswith('#') | line.startswith('!'):
                continue
            values = line.split()
            if len(values)>0:
                if values[0] == 'startNight':
                    startdates.append(int(values[2]))
                if values[0] == 'duration':
                    durations.append(int(values[2]))                
        file.close()
        # Translate startNight & duration into a list of dates covering all downtime in the survey.
        self.downdates[s] = []
        for start, dur in zip(startdates, durations):
            for i in range(0, dur):
                self.downdates[s].append(start+i)
        self.downdates[s] = numpy.array(self.downdates[s], int)
        #print self.downdates[s]
        # Check the total amount of data (mostly for user awareness):
        print '# Read %d downtime nights from %s file. ' %(len(self.downdates[s]), filename)

    def checkDownstatus(self, mjd, config):
        """Check whether telescope is 'open' or not. """    
        # Adjust time in mjd to time contained in downtime records (which counts from 'start' only).
        night = mjd - config['sim_start']
        # Calculate the 'night' of the survey to match against integer values in downtime.
        night = numpy.floor(night - config['midnight'] - 0.5)
        if isinstance(night, float):
            night = numpy.array([night,], float)
        downstate = numpy.zeros(len(night), int)
        for i in range(len(night)):
            if night[i] in self.alldowndates:
                downstate[i] = 1
        if len(night) == 1:
            downstate = downstate[0]
        return downstate
            
                    
        
            
