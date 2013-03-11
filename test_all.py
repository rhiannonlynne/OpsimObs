""" Pull all observations from an opsim run, and test against opsimObs.py """

import os
import useful_input as ui
import numpy
import pylab

opsim= 'output_opsim4_262'
nobs= 10000

def get_opsimObs_inputs(opsim, offset=0, drows=nobs):
    conn, cursor = ui.sqlConnect()
    query = 'select fieldRA*180/PI(), fieldDec*180/PI(), expmjd, filter, exptime from %s group by expmjd' %(opsim)
    results = ui.sqlQueryLimited(cursor, query, offset=offset, drows=drows)
    keys = ('ra', 'dec', 'expmjd', 'filter', 'exptime')
    keytypes = ('float', 'float', 'float', 'str', 'int')
    data = ui.assignResults(results, keys)
    ui.writeDatafile('test_%s' %(opsim), data, keys)
    ui.sqlEndConnect(conn, cursor)
    return


def get_opsimValues(opsim, offset=0, drows=nobs):
    conn, cursor = ui.sqlConnect()
    query = 'select fieldRA*180/PI(), fieldDec*180/PI(), expmjd, filter, altitude*180/PI(), azimuth*180/PI(), airmass, xparency, seeing, sunAlt*180/PI(), sunAz*180/PI(), MoonAlt*180/PI(), moonPhase, moonIllum, skybrightness_modified, 5sigma_modified from %s group by expmjd' %(opsim)
    results = ui.sqlQueryLimited(cursor, query, offset=offset, drows=drows)
    keys = ('ra', 'dec', 'expmjd', 'filter', 'alt', 'az', 'airmass', 'cloud', 'seeing', 'sunalt', 'sunaz', 'moonalt', 'moonphase', 'moonillum', 
            'skybright', 'maglimit')
    opsimdata = ui.assignResults(results, keys)
    ui.sqlEndConnect(conn, cursor)
    return opsimdata


def get_opsimObsValues(filename):
    keys = ('ra', 'dec', 'expmjd', 'filter', 'alt', 'az', 'airmass', 'cloud', 'seeing', 'sunalt', 'sunaz', 'moonalt', 'moonaz', 'moonphase', 'moonillum', 
            'skybright', 'maglimit', 'surveynight', 'downtime')
    keytypes = ('float', 'float', 'float', 'str', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 
                'float', 'float', 'int', 'int')
    opsimObsdata = ui.readDatafile(filename, keys, keytypes)
    return opsimObsdata


def compare(opsim, opsimObs, nobs=nobs):
    compkeys = ('alt', 'az', 'airmass', 'cloud', 'seeing', 'sunalt', 'sunaz', 'moonalt', 'moonphase', 'skybright', 'maglimit')
    for k in compkeys:
        pylab.figure()
        if k == 'moonphase':
            k2 = 'moonillum'
        else:
            k2 = k
        print k, len(opsim[k]), len(opsimObs[k2])
        diff = opsim[k] - opsimObs[k2]
        pylab.hist(diff, bins=nobs / 10)
        pylab.title('%s' % (k))
    # Look in more depth at airmass differences. 
    pylab.figure()
    pylab.plot(opsim['airmass'], opsimObs['airmass'], 'b.')
    pylab.xlabel('opsim airmass')
    pylab.ylabel('opsimObs airmass')
    # Look more depth at moon phase. 
    pylab.figure()
    pylab.plot(opsim['moonphase'], opsimObs['moonillum'], 'r.')
    pylab.xlabel('opsim moon phase')
    pylab.ylabel('opsimObs moon phase')
    # Look more depth at sky brightness / limiting mag.
    pylab.figure()
    pylab.plot(opsim['skybright'], opsimObs['skybright'], 'k.')
    pylab.xlabel('Opsim skybright')
    pylab.ylabel('OpsimObs skybright')
    pylab.title('Comparison of skybrightness estimates')
    return
