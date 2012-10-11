"""
opsimObs.py :

A python script intended to compute 'as would be (approximately) seen in opSim' observation conditions
for a given list of RA (degrees)/Dec (degrees)/Time (UTC MJD)/ExpTime values.
The script will read weather data (cloud and seeing histories from Cerro Tololo),
compute seeing appropriate for the bandpass & airmass of each observation,
compute sky brightness (including telescope dome scattered light contributions) for each observation,
compute the 5-sigma limiting magnitude for each observation,
and report these values to the user.

Requirements:
python
numpy

"""

import sys
import warnings
import numpy
from Weather import Weather
from Downtime import Downtime
from SkyPos import SkyPos
from ObsFields import ObsFields
from Sun import Sun
from Moon import Moon
from SkyBright import SkyBright 
from m5calculations import m5calculations 

_deg2rad = numpy.pi/180.0
_rad2deg = 180.0/numpy.pi
_sec2day = 1/60.0/60.0/24.0
_day2sec = 24.0*60.0*60.0

import time
def dtime(time_prev):
    return (time.time() - time_prev), time.time()
    


def setDefaultConfigs(override_file = None):
    """ Set some default configuration parameters for the LSST. The override_file can be used to override
    these parameters - file format should be 'keyword' 'value' (keyword should match one of the dictionary
    keys below and value should be in quotes or have no spaces).
    Returns dictionaries 'config', 'obsreq', 'lsstsite' containing needed parameters. """
    # These next parameters are used to define configuration parameters for opsimObs. 
    config = {}
    config_keys = ['sim_start', 'seeing_start', 'cloud_start', 'seeing_datafile', 'cloud_datafile',
                   'filter_names', 'seeing_wavelen', 'filter_wavelen',
                   'schedDowntime_datafile', 'unschedDowntime_datafile']
    config_keytypes = [float, float, float, str, str, tuple, float, dict, str, str]
    # Simulation 'start' day (UTC MJD) for placing location of Earth / visibility of fields.
    # OpSim starts this on 49353 (to match with start of the weather data).  (sim_start == seeing_start ..)
    config['sim_start'] = 49353.0
    # Seeing data 'start' day (can offset into the seeing data by this many days - will loop if run past end). 
    config['seeing_start'] = 0.0 
    # Cloud data 'start' day (can offset into the cloud data by this many days - will loop if run past end).
    config['cloud_start'] = 0.0 
    # Seeing data file. Default file contains 2 years of seeing data from Pachon, recorded every 300 seconds.
    config['seeing_datafile'] = 'SeeingPachon.txt'
    # Cloud data file. Default file contains 10 years of cloud data from CTIO, sampled 4 times/night.
    config['cloud_datafile'] = 'CloudTololo.txt'
    # List of the LSST filters. Default 'y' = 'y4' (as values for 'y' == 'y4').
    config['filter_names'] = ('u', 'g', 'r', 'i', 'z', 'y', 'y3', 'y4')
    # Effective wavelengths of each filter (for seeing corrections).
    config['seeing_wavelen'] = 500.0
    config['filter_wavelen'] = {'u':359, 'g':479, 'r':620, 'i':753, 'z':869, 
                                     'y':967, 'y3':996, 'y4':967}
    # Downtime files.
    config['schedDowntime_datafile'] = 'schedDown.conf'
    config['unschedDowntime_datafile'] = 'unschedDown.conf'
    # These next parameters are used to do a rudimentary check on if requested observations are possible
    #  based on some simple time requirements (filter changes, number of filters within a night..).
    config_keys += ['n_simultaneous_filters', 'readout_time', 'shutter_time', 'filter_time', 'slew_time',
                    'nexp_visit', 'add_shutter', 'add_slew']
    config_keytypes += [int, float, float, float, float, int, int, int]
    # Note that this check does NOT include any instrument model or slew time needed. 
    # Number of filters that can be used in the telescope in one night.
    config['n_simultaneous_filters'] = 5
    # Readout time (seconds)
    config['readout_time'] = 2.0
    # Shutter time (seconds)
    config['shutter_time'] =  1.0
    # Add shutter time to given expTime (=1) or implicitly include (=0) [opsim implicitly includes in visit expTime]
    config['add_shutter'] = 0.0    
    # Slew/settle time (seconds)
    config['slew_time'] = 5.0
    # Add slew time after each visit (for dither, for example?) [=0 or 1]. 
    config['add_slew'] = 0
    # Filter change time (seconds)
    config['filter_time'] = 120.0
    # Number of exposures per visit (1==one exp per visit) .. telescope slew/settle happens after visit, if 'add_slew'.
    #  NOTE : the input 'observation' file should only record VISITS 
    config['nexp_visit'] = 2.0
    # These next parameters are used to determine some parameters for the LSST telescope.
    config_keys += ['latitude', 'longitude', 'height', 'pressure', 'temperature', 'relativeHumidity', 'lapseRate']
    config_keytypes += [float, float, float, float, float, float, float]
    config['longitude'] = -1.2320792*_rad2deg   #degrees
    config['latitude'] = -0.517781017*_rad2deg  #degrees
    config['height'] = 2650.0  #km 
    config['pressure'] = 749.3
    config['temperature'] = 285.655 #C
    config['relativeHumidity'] = 0.40
    config['lapseRate'] =  0.0065
    config['midnight'] = 0.16  # Midnight at LSST is this much past the start of a day (UTC)
    # Contribution toward the final seeing that comes from the telescope itself (arcseconds)
    #   (the seeing data file is raw atmospheric seeing)
    config['seeing_Telescope'] = 0.323
    # And limit on 'good end' of seeing.
    config['good_seeing_limit'] = 0.25
    # Minimum distance to moon (radians)
    config['min_moonDist'] = 15.0 * _deg2rad
    # Read override file, if provided.
    if override_file != None:
        file = open(override_file, 'r')
        for line in file:
            # Skip comments or incomplete keyword / value pairs. 
            if (line.startswith('#') or line.startswith('!') or (len(line.split())<3)):
                continue
            # Get the keyword / "=" /  value pairs from this line.
            keyword = line.split()[0]
            value = line.lstrip(keyword).lstrip().lstrip('=').lstrip()
            value = value.rstrip()
            # Evaluate and override default dictionary settings.
            if keyword in config_keys:
                ix = config_keys.index(keyword)
                print '# Overriding ', keyword,' from ', config[keyword], 'to ', value, ' (a %s value)' %(config_keytypes[ix])
                config[keyword] = config_keytypes[ix](value)
            else:
                raise warnings.warn('Unrecognized keyword %s: ignoring this value.' %(keyword))
        file.close()
    return config


if __name__ == '__main__':

    inputobs_file = sys.argv[1]
    # Deal with input configuration information. 
    if len(sys.argv)> 2:
        override_config_file = sys.argv[2]
    else:
        override_config_file = None
    # Read configuration parameters.
    config = setDefaultConfigs(override_config_file)

    # Set up a skypos object to hold site information and provide ra/dec -> alt/az/airmass translations.
    skypos = SkyPos()
    skypos.setSite(lat=config['latitude'], lon=config['longitude'], height=config['height'],
                   pressure=config['pressure'], temperature=config['temperature'],
                   relativeHumidity=config['relativeHumidity'], lapseRate=config['lapseRate'])

    # Set up a Weather object to read the site weather data.
    t = time.time()
    weather = Weather()
    weather.readWeather(config)
    dt, t = dtime(t)
    print '# Reading weather required %.2f seconds' %(dt)

    # Set up a Downtime object to read the downtime data.
    downtime = Downtime()
    downtime.readDowntime(config)
    dt, t = dtime(t)
    print '# Reading downtime required %.2f seconds' %(dt)

    # Read observations.
    obs = ObsFields()
    obs.readInputObs(inputobs_file)
    # Check timing of input observations.
    obs.checkInputObs(config)
    nobs = len(obs.ra)
    dt, t = dtime(t)
    print '# Checking %d input observations required %.2f seconds' %(nobs, dt)

    # Calculate alt/az/airmass for all fields.
    obs.getAltAzAirmass(skypos)
    # Calculate weather (cloud/seeing) for all fields.
    dt, t = dtime(t)
    obs.getObsWeather(weather, config)
    dt, t = dtime(t)
    print '# Getting weather information for %d observations required %.2f seconds' %(nobs, dt)
    # Check downtime status for these observations
    obs.getDowntime(downtime, config)

    # Calculate position of sun at the times of these observations.
    sun = Sun()
    dt, t = dtime(t)
    sun.calcPos(obs.mjd)
    sun.getAltAz(skypos)
    dt, t = dtime(t)
    print '# Calculating sun position at %d times required %.2f seconds' %(nobs, dt)

    # Calculate the position, phase and altitude of the Moon. 
    moon = Moon()
    moon.calcPos(obs.mjd, config)
    moon.getAltAz(skypos)
    dt, t = dtime(t)
    print '# Calculating moon position at %d times required %.2f seconds' %(nobs, dt)

    # Will clean this up and put into classes as time is available. 
    # Calculate the sky brightness. 
    skybright = SkyBright(model='Perry', solar_phase='ave')    
    sky = numpy.zeros(len(obs.mjd), 'float')
    for i in range(len(obs.mjd)):
        # Calculate sky brightness for each observation. 
        skybright.setSkyBright(obs.alt[i], obs.az[i], moon.alt[i], moon.az[i], moon.phase[i], 
                               bandpass=obs.filter[i])
        sky[i] = skybright.getSkyBright()
        # Add modification to match 'skybrightness_modified' (which is brighter in twilight)
        sky = numpy.where(sun.alt > -18, 17, sky)
    dt, t = dtime(t)
    print '# Calculating the sky brightness for %d observations required %.2f seconds' %(nobs, dt)

    # Calculate the 5-sigma limiting magnitudes. 
    maglimit = numpy.zeros(len(obs.mjd), 'float')
    m5 = m5calculations()
    # Read the throughput curves for LSST (needed to determine zeropoints). 
    m5.setup_Throughputs(verbose=False)
    # Swap 'y4' for 'y' (or other default y value). 
    tmp_filters = m5.check_filter(obs.filter)
    # Determine the unique exposure times, as have to set up dark current, etc. based on this for telescope ZP's.
    exptimes = numpy.unique(obs.exptime)
    for expT in exptimes:
        condition = [obs.exptime == expT]
        # Calculate telescope zeropoints.        
        opentime = ((expT - config['readout_time']*config['nexp_visit'] -  config['nexp_visit']* config['add_shutter']*config['shutter_time']) 
                    / config['nexp_visit'])
        print "# Calculating depth for %d exposures of %.2f open shutter time" %(config['nexp_visit'], opentime)
        m5.setup_values(expTime=opentime, nexp=config['nexp_visit'])
        # Calculate 5sigma limiting magnitudes. 
        maglimit[condition] = m5.calc_maglimit(obs.seeing[condition], sky[condition], tmp_filters[condition], obs.airmass[condition], snr=5.0)        
    dt, t = dtime(t)
    print '# Calculating the m5 limit for %d observations required %.2f seconds' %(nobs, dt)
        
    # Print out interesting values. 
    obs.printObs(sun, moon, sky, maglimit,  config)
    


    
