import os
import numpy
from Bandpass import Bandpass
from Sed import Sed

class m5calculations():
    def __init__(self, default_y='y4'):
        """Instantiate the object."""
        self.default_y = default_y
        return

    def setup_Throughputs(self, rootDir=None, verbose=True):
        """Read bandpasses and dark sky sed.
        Call this method first before doing anything else. """
        self.filterlist = ('u', 'g', 'r', 'i', 'z', 'y3', 'y4')
        if rootDir == None:
            rootDir = os.getenv('LSST_THROUGHPUTS_DEFAULT')            
        if rootDir == None:
            raise Exception('Either provide "rootDir" or set $LSST_THROUGHPUTS_DEFAULT')
        # Read in the total transmission curves. 
        self.lsst = {}
        hardwarecomponentlist = ['detector.dat','lens1.dat', 'lens2.dat', 'lens3.dat',
                                 'm1_ProtAl_Aged.dat', 'm2_ProtAl_Aged.dat', 'm3_ProtAl_Aged.dat']
                                 #'m1_ProtAl_Ideal.dat', 'm2_ProtAl_Ideal.dat', 'm3_ProtAl_Ideal.dat']
                                 #'m1_ProtAlAg.dat', 'm2_ProtAlAg.dat', 'm3_ProtAlAg.dat']
        commoncomponentlist = hardwarecomponentlist + ['atmos_10.dat',]
        for f in self.filterlist:
            componentlist = commoncomponentlist + ['filter_' +f +'.dat',]
            self.lsst[f] = Bandpass()
            self.lsst[f].readThroughputList(componentlist, rootDir=rootDir, verbose=verbose)
        # Read in the hardware-only transmission curves (no atmosphere).
        self.hardware = {}
        for f in self.filterlist:
            componentlist = hardwarecomponentlist + ['filter_' +f +'.dat',]
            self.hardware[f] = Bandpass()
            self.hardware[f].readThroughputList(componentlist, rootDir=rootDir, verbose=verbose)
        # Read in the dark sky SED.
        darkskyfile = 'darksky.dat'
        self.darksky = Sed()
        self.darksky.readSED_flambda(os.path.join(rootDir, darkskyfile))
        # Set up a flat SED.
        self.flatsed = Sed()
        self.flatsed.setFlatSED()
        return        
    
    def setup_values(self, expTime=15.0, nexp=2., instrument_noise=None, instrument_noise_visit=None):
        """Set appropriate default values. Allows user to change expTime, nexp, or instrument_noise.
        Call this method second - and potentially many times.
        The exposure time is used to calculate dark current noise, m5 and Cm values, as well as ZP_t/ZP_h.
        These should be separated out when I have more time, to allow variable expTime in calcM5 below. """
        # exposure time (seconds)
        self.expTime = expTime
        # number of exposures per visit (for calculating total m5 per visit rather than per exposure)
        self.nexp = nexp
        # Gain, electrons/adu.
        self.gain = 2.3
        if instrument_noise != None:
            self.instrument_noise = instrument_noise
            self.othernoise = 0.0
            self.rdnoise = self.instrument_noise
        elif instrument_noise_visit != None:
            self.instrument_noise = instrument_noise_visit / numpy.sqrt(self.nexp)
            self.othernoise = 0.0
            self.rdnoise = self.instrument_noise
        else:
            # othernoise == camera electronics noise (electrons/exposure/pixel)
            self.othernoise = 3.9   
            # actual readnoise straight off the camera (electrons/exposure/pixel)
            self.rdnoise = 5.9
            # total instrumental noise (camera readnoise + other). electrons/exposure/pixel.
            #  add 0.5 * gain to (allow for potential undersampling of readnoise. with gain*0.5 term)
            #self.instrument_noise = numpy.sqrt(self.othernoise**2 + self.rdnoise**2 + (self.gain*0.5)**2)
            self.instrument_noise = numpy.sqrt(self.othernoise**2 + self.rdnoise**2)
        # Dark current, electrons/pix/second.
        self.dark_current = 0.2         
        # Plate scale, arcseconds per pixel.
        self.platescale = 0.2
        # Telescope primary mirror diameter, effective collecting area (cm^2). 
        self.effarea = numpy.pi*(6.5*100/2.0)**2
        # Set default y band (for inputs where only 'y' is specified). 
        self.default_y = 'y4'
        # Calculate Tb / Sb (integral of transmission/wavelength) as in LSE-40 (Table 7), to calculate kAtm.
        Tb = {}
        Sb = {}
        self.kAtm = {}
        stepsize = self.lsst[self.filterlist[0]].wavelen[1] - self.lsst[self.filterlist[0]].wavelen[0]
        for f in self.filterlist:
            self.kAtm[f] = [0.0, 0.0]
            Tb[f] = (self.lsst[f].sb / self.lsst[f].wavelen)*stepsize
            Sb[f] = (self.hardware[f].sb / self.hardware[f].wavelen)*stepsize
            if f.startswith('y'):
                condition = (self.lsst[f].wavelen >= 975) | (self.lsst[f].wavelen <= 800)
                self.kAtm[f][0] = -2.5*numpy.log10(Tb[f][condition].sum() / Sb[f][condition].sum())
                condition = (self.lsst[f].wavelen > 800) & (self.lsst[f].wavelen < 975)
                self.kAtm[f][1] = -2.5*numpy.log10(Tb[f][condition].sum() / Sb[f][condition].sum())
            else:
                self.kAtm[f][0] = -2.5*numpy.log10(Tb[f].sum() / Sb[f].sum())
        # Calculate telescope zeropoints. 
        self.zpT = {}
        self.zpH = {}
        for f in self.filterlist:
            self.zpT[f] = self.lsst[f].calcZP_t(expTime=self.nexp*self.expTime, effarea=self.effarea,
                                                gain=self.gain)
            self.zpH[f] = self.hardware[f].calcZP_t(expTime=self.nexp*self.expTime, effarea=self.effarea,
                                                    gain=self.gain)
        return

    def print_values(self):
        """Print a report of the values used in the maglimit calculations."""
        # Calculate some additional stuff which can be printed out. 
        self.C_m = {}
        self.darkskymags = {}
        self.m5 = {}
        for f in self.filterlist:
            self.darkskymags[f] = self.darksky.calcMag(self.hardware[f])
            self.m5[f] = self.lsst[f].calcM5(self.darksky, self.hardware[f], expTime=self.expTime,
                                             nexp=self.nexp, readnoise=self.rdnoise,
                                             othernoise=self.othernoise, darkcurrent=self.dark_current,
                                             gain=self.gain, effarea=self.effarea, 
                                             seeing = 0.7, platescale = self.platescale)
            self.C_m[f] = self.m5[f] - 0.50*(self.darkskymags[f] - 21.0)
        # Start printing stuff to screen.
        print 'Exposure time ', self.expTime
        print 'Number of exposures per visit ', self.nexp
        print 'Gain ', self.gain
        print 'Camera readnoise ', self.rdnoise, ' and other readnoise ', self.othernoise
        print 'Instrumental noise per exposure ', self.instrument_noise
        print 'Dark current per second per pixel ', self.dark_current
        print 'Dark current per exposure ', self.dark_current*self.expTime
        print 'Total camera noise per visit (e/pix/visit) ', \
            numpy.sqrt(self.nexp*self.expTime*self.dark_current + self.nexp*(self.instrument_noise)**2)
        print 'Platescale ', self.platescale
        print 'Telescope effective area ', self.effarea        
        print 'Scaling relation C_m and kAtm values: '
        for f in self.filterlist:
            print '\t \tin filter ', f, ' C_m=', self.C_m[f], ' kAtm= ', self.kAtm[f]
        print 'C_m=', self.C_m
        print 'kAtm=', self.kAtm
        print 'Telescope and Hardware zeropoints (%f sec visit):' %(self.expTime*self.nexp)
        for f in self.filterlist:
            print '\t\tin filter ', f, ' zpT = ', self.zpT[f], ' zpH = ', self.zpH[f]
        print 'Dark sky m5 limits :'
        for f in self.filterlist:
            print '\t\tin filter ', f, ' m5 = ', self.m5[f]
        print 'Dark sky noise (seeing=0.7 arcsec, @ zenith) & camera noise (both in electrons | ADU):'
        for f in self.filterlist:
            skynoise = numpy.sqrt(self.darksky.calcADU(self.hardware[f], expTime=self.nexp*self.expTime,
                                                        gain=self.gain, effarea=self.effarea) \
                                                        * self.platescale**2  * self.gain)  # electrons
            instnoise = numpy.sqrt(self.nexp *(self.dark_current*self.expTime + self.instrument_noise**2)) #electrons
            print '\t\tin filter ', f, ' skynoise = ', skynoise, '|', skynoise/self.gain, ' instnoise = ', instnoise, '|', instnoise/self.gain
        return
      
    def check_filter(self, filter):
        """Check filter array for consistency with internal set. Basically this means replace 'y' with
        the default y band choice. """
        filter = numpy.where(filter == 'y', self.default_y, filter)
        filters_used = numpy.unique(filter)
        for f in filters_used:
            if f not in self.filterlist:
                raise Exception('I do not recognize filter %s' %(f))
        return filter
    
    def calc_maglimit(self, seeing, skybrightness, filter, airmass, snr=5.0):
        """Calculate limiting magnitude at snr, for nexp/expTime/gain/instNoise/zeropoint values of class,
        under conditions of seeing (arcseconds) and skybrightness (mag/arcsecond^2).
        Returns mag limit. """        
        # neff = 'effective' pixel area for a point source. (see LSE-40 - SNR doc - eqn 31).
        neff = 2.436 * (seeing/self.platescale)**2
        # Calculate sky counts (counts/pixel/exp) in this bandpass from the skybrightness (mag/arcsecond^2).
        filters_used = numpy.unique(filter)
        skycounts = numpy.zeros(len(skybrightness), float)
        for f in filters_used:
            condition = (filter == f)
            # Convert skycounts from mags/''sq to counts/''sq for visit.
            skycounts[condition] = (10.**(-0.4*(skybrightness[condition] - self.zpH[f])))
        # Convert to skycounts per pixel. 
        skycounts = skycounts * self.platescale**2   #(mag/pixel)
        # Calculate the sky noise (squared) in ADU.
        skynoise_sq = skycounts / self.gain
        # Calculate noise (squared) from the instrument, converting result to ADU.
        instnoise_sq = self.nexp*(self.dark_current*self.expTime + self.instrument_noise**2) / (self.gain**2.0)
        # see equation 42 from SNR doc
        noise_sq = (skycounts / self.gain + instnoise_sq) * neff
        # Translate this to the required counts for a source using equations 45/46 from SNR doc.
        #  Counts are in ADU using this formula. 
        counts = (snr**2.0)/self.gain/2.0 + numpy.sqrt((snr**4.0)/(self.gain)**2/4.0 + (snr**2.0)*noise_sq)
        # And translate to counts, at this airmass in this filter.
        #    Note we did not have to correct for extinction for the skycounts, because the sky background
        #    is already extinction-corrected (which is part of the reason we must use hardware ZP only).
        mags = numpy.zeros(len(counts), 'float')
        for f in filters_used:
            condition = (filter == f)
            # Convert to magnitudes
            mags[condition] = -2.5*numpy.log10(counts[condition]) + self.zpT[f]
            # Correct for atmospheric extinction (note there is already a factor of X=1 in the zeropoint). 
            mags[condition] = (mags[condition] - self.kAtm[f][0]*(airmass[condition]-1)
                               - self.kAtm[f][1]*numpy.sqrt(airmass[condition]-1))
        return mags

            
if __name__ == '__main__':
    m5s = m5calculations()
    m5s.setup_Throughputs()
    m5s.setup_values(instrument_noise=9)
    m5s.print_values()
    exit()

    """
    skybrightness = numpy.arange(18, 24, 0.2)
    airmass = numpy.ones(len(skybrightness), float)
    seeing = airmass * 0.7
    filterlist = ('u', 'g', 'r', 'i', 'z', 'y4')
    mag = {}
    m5 = {}
    colors = {'u':'b', 'g': 'g', 'r':'r', 'i':'m', 'z':'k', 'y4':'y'}
    for f in filterlist:
        filter = numpy.empty(len(skybrightness), 'a3')
        for i in range(len(filter)):    
            filter[i] = f
        mag[f] = m5s.calc_maglimit(seeing, skybrightness, filter, airmass)
        m5[f] = m5s.calc_maglimit_Cm(seeing,skybrightness, filter, airmass)
    import pylab
    for f in filterlist:
        pylab.plot(skybrightness, mag[f], linestyle='-', color=colors[f], label=f)
        pylab.plot(skybrightness, m5[f], linestyle=':', color=colors[f])
    pylab.legend(numpoints=1, loc='lower right')
    pylab.xlabel('Skybrightness')
    pylab.ylabel('mag @ SNR=5')
    pylab.show()
    """
    
    # Test with opsim values.
    import time
    def dtime(time_prev):
        return (time.time() - time_prev), time.time()
    
    # read test data
    import useful_input as ui
    filename ='test_m5'
    #filename ='t'
    keys = ('seeing', 'airmass', 'skybright', 'filter', '5sigma')
    keytypes = ('float', 'float', 'float', 'a3' ,'float')
    data = ui.readDatafile(filename, keys, keytypes)
    
    data['filter'] = m5s.check_filter(data['filter'])
    t = time.time()
    mags = m5s.calc_maglimit(data['seeing'], data['skybright'], data['filter'], data['airmass'])
    dt, t = dtime(t)
    print 'Calculate full SNR mag limits for %d objects in %f s' %(len(data['filter']), dt)
    m5 = m5s.calc_maglimit_Cm(data['seeing'], data['skybright'], data['filter'], data['airmass'])
    dt, t = dtime(t)
    print 'Calculate Cm style for %d objects in %f s' %(len(data['filter']), dt)
    import pylab
    pylab.figure()
    filters_used = numpy.unique(data['filter'])
    for f in filters_used:
        condition = (data['filter'] == f)
        pylab.plot(data['skybright'][condition], (mags[condition] - m5[condition]),
                   marker='.', linestyle='', label=f)
    pylab.legend(loc=(0.9, 0.7), numpoints=1)
    pylab.xlabel('Sky Brightness')
    pylab.ylabel('Full SNR 5sigma - C_m calculated 5sigma')
    pylab.grid()
    pylab.show()
