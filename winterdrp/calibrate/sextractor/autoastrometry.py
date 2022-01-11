#!/bin/env python

#  autoastrometry.py - a fast astrometric solver
#
#    author: Daniel Perley (dperley@astro.caltech.edu)
#    last significant modifications 2012-04-23
#  
#  Installation:
#     Save this file anywhere on disk, and call it from the command 
#       line: "python autoastrometry.py"
#     Required python packages:  numpy, pyfits, and optionally ephem.
#     You must also have sextractor installed: if the path is
#       nonstandard, edit the global variable below to specify. 
#     For help, type "python autoastrometry.py -help"

# 4/23: program can actually be overwhelmed by too many good matches (too high maxrad).
# need to fix this.

# Modified by Kishalay De (kde@astro.caltech.edu) for removing dependency on deprecated pyfits
# and making it compatible with astropy headers and python 3.6 (June 11, 2018)

import math
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
from math import sin, cos, tan, asin, sqrt
import numpy as np
from astropy.io import fits as af
from winterdrp.paths import base_output_dir
from winterdrp.calibrate.sextractor.sourceextractor import execute_sextractor, sextractor_cmd
import logging
import ephem
from winterdrp.calibrate.sextractor.settings import write_param_file, write_config_file, default_config_path

logger = logging.getLogger(__name__)

default_tolerance = 0.01  # these defaults should generally not be altered.
defaultpatolerance = 1.4
defaultminfwhm = 1.5
defaultmaxfwhm = 40
default_saturation = 1.e10

fastmatch = 1
showmatches = 0


class AstrometryException(Exception):
    pass


class Obj:
    ra = 0.0
    dec = 0.0
    mag = 0.0

    ra_rad = 0.0
    dec_rad = 0.0

    def __init__(
            self,
            inra: float,
            indec: float,
            inmag: float
    ):
        self.ra = float(inra)
        self.dec = indec
        self.ra_rad = inra * math.pi/180
        self.dec_rad = indec * math.pi/180
        self.mag = inmag

    def rotate(self, dpa_deg, ra0, dec0):
        dpa_rad = dpa_deg * math.pi/180
        sindpa = sin(dpa_rad)
        cosdpa = cos(dpa_rad)
        rascale = cos(dec0*math.pi/180)

        #this is only valid for small fields away from the pole.
        x = (self.ra - ra0) * rascale
        y = (self.dec - dec0)

        xrot = cosdpa * x - sindpa * y
        yrot = sindpa * x + cosdpa * y

        self.ra = (xrot / rascale) + ra0
        self.dec = yrot + dec0
        self.ra_rad = self.ra * math.pi/180
        self.dec_rad = self.dec * math.pi/180


class SexObj(Obj):
    x = 0.
    y = 0.
    mag = 0.0
    magerr = 0.0
    ellip = 0.0
    fwhm = 0.0
    flag = 0

    def __init__(self, inline):
        inlinearg = [x.strip() for x in inline.split(" ") if x not in [""]]

        if len(inlinearg) < 8:
            err = f"Expected 8 values in table, found {len(inlinearg)} ({inlinearg})"
            logger.error(err)
            raise ValueError(err)

        self.x = float(inlinearg[0])
        self.y = float(inlinearg[1])

        super().__init__(*[float(x) for x in inlinearg[2:5]])

        self.magerr = float(inlinearg[5])
        self.ellip = float(inlinearg[6])
        self.fwhm = float(inlinearg[7])

        if len(inlinearg) >= 9:
            self.flag = int(inlinearg[8])
        else:
            self.flag = None


#Pixel distance
def imdistance(obj1, obj2):
    return ((obj1.x - obj2.x)**2 + (obj1.y - obj2.y)**2)**0.5


#Great circle distance between two points.
def distance(obj1, obj2):
    # both must be Obj's.

    ddec = obj2.dec_rad - obj1.dec_rad
    dra  = obj2.ra_rad - obj1.ra_rad
    dist_rad = 2 * asin(sqrt( (sin(ddec/2.))**2 + cos(obj1.dec_rad) * cos(obj2.dec_rad) * (sin(dra/2.))**2))

    dist_deg = dist_rad * 180. / math.pi
    dist_sec = dist_deg * 3600.
    return dist_sec


#Non-great-circle distance is much faster
def quickdistance(obj1, obj2, cosdec):
    ddec = obj2.dec - obj1.dec
    dra = obj2.ra - obj1.ra
    if dra > 180:
        dra = 360 - dra
    return 3600 * sqrt(ddec**2 + (cosdec*dra)**2)


#Calculate the (spherical) position angle between two objects.
def posangle(obj1, obj2):
    dra = obj2.ra_rad - obj1.ra_rad
    pa_rad = np.arctan2(cos(obj1.dec_rad) * tan(obj2.dec_rad) - sin(obj1.dec_rad) * cos(dra), sin(dra));
    pa_deg = pa_rad * 180./math.pi;
    pa_deg = 90. - pa_deg  #defined as degrees east of north
    while pa_deg > 200: pa_deg -= 360.   # make single-valued
    while pa_deg < -160: pa_deg += 360.  # note there is a crossing point at PA=200, images at this exact PA
    return pa_deg                        # will have the number of matches cut by half at each comparison level

#Compare objects using magnitude.
def magcomp(obj): #useful for sorting; Altered by KD for compatibility with python 3
    return obj.mag
    #return (obj1.mag > obj2.mag) - (obj1.mag < obj2.mag)


#Check if two values are the same to within a fraction specified.
def fuzzyequal(v1, v2, tolerance):
    return abs(v1/v2 - 1) < tolerance


def median(l):
    a = np.array(l)
    return np.median(a)


def stdev(l):
    a = np.array(l)
    return np.std(a)

def mode(l):
    if len(l) == 0:
        return
    s = np.array(sorted(l))
    d = s[1:] - s[:-1]
    nd = len(d)
    if nd >= 32:
        g = nd/16
    elif nd >= 6:
        g = 2
    else:
        g = 1

    #g = max(nd / 16,1)  #sensitive to clusters up to a little less than 1/16 of the data set
    minmean = d.sum()
    imean = nd / 2

    for i in range(nd):
        r = [int(max(i-g, 0)), int(min(i+g, nd))]
        m = d[r[0]:r[1]].mean()
        if m < minmean:
            minmean = m
            imean = i

    mode = s[int(imean)] #+ s[imean+1])/2
    return mode


def rasex2deg(rastr):
    rastr = str(rastr).strip()
    ra=rastr.split(':')
    if len(ra) == 1: return float(rastr)
    return 15*(float(ra[0])+float(ra[1])/60.0+float(ra[2])/3600.0)


def decsex2deg(decstr):
    decstr = str(decstr).strip()
    dec=decstr.split(':')
    if len(dec) == 1:
        return float(decstr)
    sign = 1
    if (decstr[0] == '-'):
        sign=-1
    return sign*(abs(float(dec[0]))+float(dec[1])/60.0+float(dec[2])/3600.0)


def unique(inlist):
    lis = inlist[:] #make a copy
    lis.sort()
    llen = len(lis)
    i = 0
    while i < llen-1:
        if lis[i+1] == lis[i]:
            del lis[i+1]
            llen = llen - 1
        else:
            i = i + 1
    return lis


def sextract(
        sexfilename,
        nxpix,
        nypix,
        border: float = 3.,
        corner: float = 12.,
        minfwhm: float = defaultminfwhm,
        maxfwhm: float = defaultmaxfwhm,
        maxellip: float = 0.5,
        saturation: float = default_saturation,
        output_dir: str = base_output_dir,
        config_path: str = default_config_path
):
    cmd = f"{sextractor_cmd} {sexfilename} -c {config_path} -SATUR_LEVEL {saturation}"
    execute_sextractor(cmd, output_dir=output_dir)

    # Read in the sextractor catalog
    with open(os.path.join(output_dir, "temp.cat"), 'rb') as cat:
        catlines = [x.replace(b"\x00", b"").decode() for x in cat.readlines()][1:]

    if len(catlines) == 0:
        logger.error('Sextractor catalog is empty: try a different catalog?')
        sys.exit(1)

    minx = border
    miny = border
    maxx = nxpix - border    # This should be generalized
    maxy = nypix - border

    l = -1
    nsexinit = 0
    nsexpass = 0
    xlist = []
    ylist = []
    sexlist = []
    fwhmlist = []
    elliplist = []
    flaglist = []

    rejects = []

    while l < len(catlines)-1:

        l += 1
        if len(catlines[l]) <= 1 or catlines[l][0] == '#':
            continue

        iobj = SexObj(catlines[l]) #process the line into an object
        nsexinit += 1

        # Initial filtering
        if iobj.ellip > maxellip:
            rejects.append("ellip")
            continue
        if iobj.fwhm < minfwhm:
            rejects.append("min fwhm")
            continue
        if iobj.fwhm > maxfwhm:
            rejects.append("max fwhm")
            continue
        if iobj.x < minx:
            rejects.append("min x")
            continue
        if iobj.y < miny:
            rejects.append("min y")
            continue
        if iobj.x > maxx:
            rejects.append("max x")
            continue
        if iobj.y > maxy:
            rejects.append("max y")
            continue
        if iobj.x + iobj.y < corner:
            rejects.append("corner")
            continue
        if iobj.x + (nypix-iobj.y) < corner:
            rejects.append("corner")
            continue
        if (nxpix-iobj.x) < corner:
            rejects.append("corner")
            continue
        if (nxpix-iobj.x) + (nypix-iobj.y) < corner:
            rejects.append("corner")
            continue
        if saturation is not None:
            if iobj.flag > 0:
                rejects.append("saturation")
                continue  # this will likely overdo it for very deep fields.

        sexlist.append(iobj)
        xlist.append(iobj.x)
        ylist.append(iobj.y)
        fwhmlist.append(iobj.fwhm)
        elliplist.append(iobj.ellip)
        flaglist.append(iobj.flag)
        nsexpass += 1

    # Remove detections along bad columns

    threshprob = 0.0001
    ctbadcol = 0
    for i in range(5):
        txp = 1.0
        xthresh = 1
        while txp > threshprob:
            txp *= min((len(sexlist)*1.0/nxpix), 0.8) # some strange way of estimating the threshold.
            xthresh += 1                          #what I really want is a general analytic expression for
        removelist = []                           #the 99.99% prob. threshold for value of n for >=n out
        modex = mode(xlist)                       #of N total sources to land in the same bin (of NX total bins)
        for j in range(len(sexlist)):
            if (sexlist[j].x > modex-1) and (sexlist[j].x < modex+1):
                removelist.append(j)
        removelist.reverse()
        if len(removelist) > xthresh:
            for k in removelist:
                del xlist[k]
                del ylist[k]
                del sexlist[k]
                del fwhmlist[k]
                del elliplist[k]
                del flaglist[k]
                ctbadcol += 1

        typ = 1.0
        ythresh = 1
        while typ > threshprob:
            typ *= min((len(sexlist)*1.0/nypix),0.8)
            ythresh += 1
        removelist = []
        modey = mode(ylist)
        for j in range(len(sexlist)):
            if (sexlist[j].y > modey-1) and (sexlist[j].y < modey+1):
                removelist.append(j)
        removelist.reverse()
        if len(removelist) > ythresh:
            for k in removelist:
                del xlist[k]
                del ylist[k]
                del sexlist[k]
                del fwhmlist[k]
                del elliplist[k]
                del flaglist[k]
                ctbadcol += 1
    if ctbadcol > 0:
        rejects += ["bad columns" for _ in range(ctbadcol)]


    # Remove galaxies and cosmic rays

    if len(fwhmlist) > 5:
        fwhmlist.sort()
        fwhm20 = fwhmlist[int(len(fwhmlist)/5)]
        fwhmmode = mode(fwhmlist)
    else:
        fwhmmode = minfwhm
        fwhm20 = minfwhm

    # formerly a max, but occasionally a preponderance of long CR's could cause fwhmmode to be bigger than the stars
    refinedminfwhm = median([0.75*fwhmmode, 0.9*fwhm20, minfwhm]) # if CR's are bigger and more common than stars, this is dangerous...
    logger.debug(f'Refined min FWHM: {refinedminfwhm} pix')

    #Might also be good to screen for false detections created by bad columns/rows

    ngood = 0
    goodsexlist = []
    for sex in sexlist:
        if sex.fwhm > refinedminfwhm:
            goodsexlist.append(sex)
            ngood += 1
        else:
            rejects.append("refined min fwhm")

    # Sort by magnitude
    goodsexlist.sort(key=magcomp)

    logger.info(f'{ngood} objects detected in image {sexfilename} (a further {nsexinit - ngood} discarded)')

    reject_stats = [(x, rejects.count(x)) for x in list(set(rejects))]
    logger.info(f"Reject reasons: {reject_stats}")

    return goodsexlist


def getcatalog(catalog, ra, dec, boxsize, minmag=8.0, maxmag=-1, maxpm=60.):
    # Get catalog from USNO

    if maxmag == -1:
        maxmag = 999 #default (custom catalog)
        if catalog == 'ub2':
            maxmag = 21.0#19.5
        if catalog == 'sdss':
            maxmag = 22.0
        if catalog == 'tmc':
            maxmag = 20.0


    if (catalog =='ub2' or catalog=='sdss' or catalog=='tmc'):
        usercat = 0
        racolumn = 1
        deccolumn = 2
        magcolumn = 6
        if catalog=='tmc':
            magcolumn=3
        pmracolumn = 10
        pmdeccolumn = 11
        queryurl = "http://tdc-www.harvard.edu/cgi-bin/scat?catalog=" + catalog +  "&ra=" + str(ra) + "&dec=" + str(dec) + "&system=J2000&rad=" + str(-boxsize) + "&sort=mag&epoch=2000.00000&nstar=6400"
        #print queryurl
        cat = urllib.request.urlopen(queryurl)
        catlines = cat.readlines()
        cat.close()
        if len(catlines) > 6400-20:
            logger.warning('Reached maximum catalog query size. Gaps may be '
                           'present in the catalog, leading to a poor solution '
                           'or no solution. Decrease the search radius.')
    else:
        usercat = 1
        try:
            logger.info(f'Reading user catalog {catalog}')

            with open(catalog, 'r') as cat:
                racolumn = 0
                deccolumn = 1   # defaults
                magcolumn = -1    #  (to override, specify in first line using format #:0,1,2)
                catlines = cat.readlines()

        except:
            logger.error(f'Failed to open user catalog {catalog}. File not '
                         f'found or invalid online catalog.  Specify ub2, sdss, or tmc.')
            return []

    l = -1
    catlist = []

    while l < len(catlines)-1:
        l += 1
        inline = catlines[l].strip()

        if len(inline) <= 2:
            continue

        if inline[0:2] == '#:':
            inlinearg = inline[2:].split(',')
            racolumn = int(inlinearg[0])-1
            deccolumn = int(inlinearg[1])-1
            if len(inlinearg) > 2:
                magcolumn = int(inlinearg[2])-1
            continue

        if (inline[0] < ord('0') or inline[0] > ord('9')) and str(inline[0]) != '.':
            continue #this may be too overzealous about
        if (inline[1] < ord('0') or inline[1] > ord('9')) and str(inline[1]) != '.':
            continue # removing comments...

        inlineargByte = inline.split()
        inlinearg = [str(a, 'utf-8') for a in inlineargByte]
        narg = len(inlinearg)

        if inlinearg[racolumn].find(':') == -1:
            ra = float(inlinearg[racolumn])
        else:
            ra = rasex2deg(inlinearg[racolumn])
        if inlinearg[deccolumn].find(':') == -1:
            dec = float(inlinearg[deccolumn])
        else:
            dec = decsex2deg(inlinearg[deccolumn])

        if magcolumn >= 0 and narg > magcolumn:
            try:
                mag = float(inlinearg[magcolumn])
            except:
                mag = float(inlinearg[magcolumn][0:-2])
        else:
            mag = maxmag

        if usercat == 0 and narg > pmracolumn and narg > pmdeccolumn:
            pmra = float(inlinearg[pmracolumn])
            pmdec = float(inlinearg[pmdeccolumn])
        else:
            pmra = pmdec = 0
        #print
        #print ra, dec, mag,
        #print pmra, pmdec,
        if mag > maxmag:
            continue #don't believe anything this faint
        if mag < minmag:
            continue #ignore anything this bright

        if abs(pmra) > maxpm or abs(pmdec) > maxpm:
            continue

        iobj = Obj(ra, dec, mag) #process the line into an object
        catlist.append(iobj)

    catlist.sort(key=magcomp)

    return catlist

def distmatch(sexlist, catlist, maxrad=180, minrad=10, tolerance=0.010, reqmatch=3, patolerance=1.2,uncpa=-1):

    if tolerance <= 0:
        logger.debug('Tolerance cannot be negative!!!')
        tolerance = abs(tolerance)
    if reqmatch < 2:
        logger.debug('Warning: reqmatch >=3 suggested')
    if patolerance <= 0:
        logger.debug('PA tolerance cannot be negative!!!')
        patolerance = abs(patolerance)
    if uncpa is None:
        uncpa = 720

    declist = []
    for s in sexlist:
        declist.append(s.dec_rad)
    avdec_rad = median(declist)       # faster distance computation
    rascale = cos(avdec_rad)          # will mess up meridian crossings, however

    #Calculate all the distances

    #print 'Calculating distances...'
    #dtime0 = time.clock()

    # In image catalog:
    sexdists = []
    sexmatchids = []
    for i in range(len(sexlist)):
        d = []
        dj = []
        for j in range(len(sexlist)):
            if i == j:
                continue
            if abs(sexlist[i].dec - sexlist[j].dec) > maxrad:
                continue
            if rascale*abs(sexlist[i].ra - sexlist[j].ra) > maxrad:
                continue
            dist = quickdistance(sexlist[i], sexlist[j], rascale)
            if dist > minrad and dist < maxrad :
                #print "%.2f" % dist, '['+str(j).strip()+']  ',
                d.append(dist)
                dj.append(j)
        sexdists.append(d)
        sexmatchids.append(dj)
        #print

    # In reference catalog:
    catdists = []
    catmatchids = []
    for i in range(len(catlist)):
        d = []
        dj = []
        for j in range(len(catlist)):
            if i == j:
                continue
            if abs(catlist[i].dec - catlist[j].dec) > maxrad:
                continue
            if rascale*abs(catlist[i].ra - catlist[j].ra) > maxrad:
                continue
            dist = quickdistance(catlist[i], catlist[j], rascale)
            if dist > minrad and dist < maxrad :
                d.append(dist)
                dj.append(j)
        catdists.append(d)
        catmatchids.append(dj)

    # Now look for matches in the reference catalog to distances in the image catalog.

    #print 'All done (in', time.clock()-dtime0, 's)'
    #print 'Finding matches...'
    #mtime0 = time.clock()

    countgreatmatches = 0

    smatch = []
    cmatch = []
    mpa = []
    nmatch = []

    primarymatchs = []
    primarymatchc = []

    for si in range(len(sexdists)):
        sexdistarr = sexdists[si]
        if len(sexdistarr) < 2:
            continue

        for ci in range(len(catdists)):
            catdistarr = catdists[ci]

            if len(catdistarr) < 2:
                continue

            match = 0
            smatchin = []
            cmatchin = []

            for sj in range(len(sexdistarr)):
                sexdist = sexdistarr[sj]
                newmatch = 1
                for cj in range(len(catdistarr)):
                    catdist = catdistarr[cj]
                    if abs((sexdist/catdist)-1.0) < tolerance:
                        match += newmatch
                        newmatch = 0 #further matches before the next sj loop indicate degeneracies
                        smatchin.append(sexmatchids[si][sj])
                        cmatchin.append(catmatchids[ci][cj])

            if match >= reqmatch:

                dpa = []
                # Here, dpa[n] is the mean rotation of the PA from the primary star of this match
                #  to the stars in its match RELATIVE TO those same angles for those same stars
                #  in the catalog.  Therefore it is a robust measurement of the rotation.

                for i in range(len(smatchin)):
                    ddpa = posangle(sexlist[si],sexlist[smatchin[i]]) - posangle(catlist[ci],catlist[cmatchin[i]])
                    while ddpa > 200: ddpa  -= 360.
                    while ddpa < -160: ddpa += 360.
                    #print smatchin[i], '-', cmatchin[i], ': ', posangle(sexlist[si],sexlist[smatchin[i]]), '-', posangle(catlist[ci],catlist[cmatchin[i]]), '=', ddpa
                    dpa.append(ddpa)

                #If user was confident the initial PA was right, remove bad PA's right away
                for i in range(len(smatchin)-1,-1,-1):
                    if abs(dpa[i]) > uncpa:
                        del smatchin[i]
                        del cmatchin[i]
                        del dpa[i]

                if len(smatchin) < 2:
                    continue

                dpamode = mode(dpa)

                #Remove deviant matches by PA
                for i in range(len(smatchin)-1, -1, -1):
                    if abs(dpa[i] - dpamode) > patolerance:
                        del smatchin[i]
                        del cmatchin[i]
                        del dpa[i]

                if len(smatchin) < 2:
                    continue

                #print si, 'matches', ci, ' ('+ str(len(smatchin))+ ' matches)'
                #print '  ',si,'-->', smatchin, '   ',ci,'-->', cmatchin
                #print '  PA change', dpamode
                #print '       PAs:', dpa

                ndegeneracies = len(smatchin)-len(unique(smatchin)) + len(cmatchin)-len(unique(cmatchin))
                # this isn't quite accurate (overestimates if degeneracies are mixed up)

                mpa.append(dpamode)
                primarymatchs.append(si)
                primarymatchc.append(ci)
                smatch.append(smatchin)
                cmatch.append(cmatchin)
                nmatch.append(len(smatchin)-ndegeneracies)

                if (len(smatchin)-ndegeneracies > 6): countgreatmatches += 1
        if countgreatmatches > 16 and fastmatch == 1: break #save processing time

        #print '   ', sexdistarr
        #print '   ', catdistarr
        #print '   ', sexlist[si].ra, sexlist[si].dec, sexlist[si].mag
        #print '   ', catlist[ci].ra, catlist[ci].dec, catlist[ci].mag
        #print '   ', distance(sexlist[si], catlist[ci])

    #print 'All done (in', time.clock()-mtime0, 's)'

    nmatches = len(smatch)
    if (nmatches == 0):
        logger.debug('Found no potential matches of any sort (including pairs).')
        logger.debug('The algorithm is probably not finding enough real stars to solve the field.  Check seeing.')
        #print 'This is an unusual error.  Check that a parameter is not set to a highly nonstandard value?'
        return [], [], []

    # Kill the bad matches
    rejects = 0

    #Get rid of matches that don't pass the reqmatch cut
    #if nmatches > 10 and max(nmatch) >= reqmatch:
    for i in range(len(primarymatchs)-1,-1,-1):
        if nmatch[i] < reqmatch:
            del mpa[i]
            del primarymatchs[i]
            del primarymatchc[i]
            del smatch[i]
            del cmatch[i]
            del nmatch[i]
            #rejects += 1  no longer a "reject"

    if len(smatch) < 1:
        logger.debug(f'Found no matching clusters of reqmatch = {reqmatch}')
        return [], [], []

    #If we still have lots of matches, get rid of those with the minimum number of submatches
    #(that is, increase reqmatch by 1)
    minmatch = min(nmatch)
    countnotmin = 0
    for n in nmatch:
        if n > minmatch:
            countnotmin += 1
    if len(nmatch) > 16 and countnotmin > 3:
        logger.debug(f'Too many matches: increasing reqmatch to {reqmatch+1}')
        for i in range(len(primarymatchs)-1,-1,-1):
            if nmatch[i] == minmatch:
                del mpa[i]
                del primarymatchs[i]
                del primarymatchc[i]
                del smatch[i]
                del cmatch[i]
                del nmatch[i]
                #rejects += 1   no longer a "reject"

    nmatches = len(smatch) # recalculate with the new reqmatch and with prunes supposedly removed
    logger.debug(f'Found {nmatches} candidate matches.')

    #for i in range(len(primarymatchs)):
    #    si = primarymatchs[i]
    #    ci = primarymatchc[i]
    #    print si, 'matches', ci, ' (dPA = %.3f)' % mpa[i]
    #    if len(smatch[i]) < 16:
    #       print '  ', si, '-->', smatch[i]
    #       print '  ', ci, '-->', cmatch[i]

    # Use only matches with a consistent PA

    offpa = mode(mpa)

    if len(smatch) > 2:

        #Coarse iteration for anything away from the mode
        for i in range(len(primarymatchs)-1,-1,-1):
            if abs(mpa[i] - offpa) > patolerance:
                del mpa[i]
                del primarymatchs[i]
                del primarymatchc[i]
                del smatch[i]
                del cmatch[i]
                del nmatch[i]
                rejects += 1

        medpa = median(mpa)
        stdevpa = stdev(mpa)
        refinedtolerance = (2.2 * stdevpa)

        #Fine iteration to flag outliers now that we know most are reliable
        for i in range(len(primarymatchs)-1,-1,-1):
            if abs(mpa[i] - offpa) > refinedtolerance:
                del mpa[i]
                del primarymatchs[i]
                del primarymatchc[i]
                del smatch[i]
                del cmatch[i]
                del nmatch[i]
                rejects += 1  #these aren't necessarily bad, just making more manageable.

    # New verification step: calculate distances and PAs between central stars of matches
    ndistflags = [0]*len(primarymatchs)
    #npaflats = [0]*len(primarymatchs)
    for v in range(2):  #two iterations
        # find bad pairs
        if len(primarymatchs) == 0:
            break

        for i in range(len(primarymatchs)):
            for j in range(len(primarymatchs)):
                if i == j: continue
                si = primarymatchs[i]
                ci = primarymatchc[i]
                sj = primarymatchs[j]
                cj = primarymatchc[j]

                sexdistij = distance(sexlist[si], sexlist[sj])
                catdistij = distance(catlist[ci], catlist[cj])

                try:
                    if abs((sexdistij/catdistij)-1.0) > tolerance:
                        ndistflags[i] += 1
                except:  # (occasionally will get divide by zero)
                    pass

        # delete bad clusters
        ntestmatches = len(primarymatchs)
        for i in range(ntestmatches-1,-1,-1):
            if ndistflags[i] == ntestmatches-1:   #if every comparison is bad, this is a bad match
                del mpa[i]
                del primarymatchs[i]
                del primarymatchc[i]
                del smatch[i]
                del cmatch[i]
                del nmatch[i]
                rejects += 1

    logger.debug(f'Rejected {rejects} bad matches.')
    nmatches = len(primarymatchs)
    logger.debug(f'Found {nmatches} good matches.')

    if nmatches == 0:
        return [], [], []


    # check the pixel scale while we're at it
    pixscalelist = []
    if len(primarymatchs) >= 2:
        for i in range(len(primarymatchs)-1):
            for j in range(i+1,len(primarymatchs)):
                si = primarymatchs[i]
                ci = primarymatchc[i]
                sj = primarymatchs[j]
                cj = primarymatchc[j]
                try:
                    pixscalelist.append(distance(catlist[ci],catlist[cj])/imdistance(sexlist[si],sexlist[sj]))
                except:
                    pass
        pixelscale = median(pixscalelist)
        pixelscalestd = stdev(pixscalelist)

        if len(primarymatchs) >= 3:
            logger.debug('Refined pixel scale measurement: %.4f"/pix (+/- %.4f)' % (pixelscale, pixelscalestd))
        else:
            logger.debug('Refined pixel scale measurement: %.4f"/pix' % pixelscale)

    for i in range(len(primarymatchs)):
        si = primarymatchs[i]
        ci = primarymatchc[i]
        # logger.debug('%3i' % si, 'matches', '%3i' % ci, ' (dPA =%7.3f)' % mpa[i])
        logger.debug(f'{si} matches {ci} (dPA =%7.3f)' % mpa[i])
        if showmatches:
            if len(smatch[i]) < 16:
                logger.debug(f'  {si} --> {smatch[i]}')
                logger.debug(f'  {ci} --> {cmatch[i]}')
            else:
                logger.debug(f'  {si} --> {smatch[i][0:10]} {len(smatch[i])-10} more')
                logger.debug(f'  {ci} --> {cmatch[i][0:10]} {len(cmatch[i])-10} more')
            if i+1 >= 10 and len(primarymatchs)-10 > 0:
                logger.debug(f''
                             f'{(len(primarymatchs)-10)} additional matches not shown.')
                break
        else:
            logger.debug(f': {str(len(smatch[i])).strip()} rays')

    with open('matchlines.im.reg', 'w') as out:

        color = 'red'
        out.write('# Region file format: DS9 version 4.0\nglobal color='+color+' font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
        out.write('image\n')
        for i in range(len(primarymatchs)):
            si = primarymatchs[i]
            for j in range(len(smatch[i])):
                sj = smatch[i][j]
                out.write("line(%.3f,%.3f,%.3f,%.3f) # line=0 0\n" % (sexlist[si].x, sexlist[si].y, sexlist[sj].x, sexlist[sj].y))

    with open('matchlines.wcs.reg', 'w') as out:
        color='green'
        out.write('# Region file format: DS9 version 4.0\nglobal color='+color+' font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
        out.write('fk5\n')
        for i in range(len(primarymatchs)):
            ci = primarymatchc[i]
            for j in range(len(smatch[i])):
                cj = cmatch[i][j]
                out.write("line(%.5f,%.5f,%.5f,%.5f) # line=0 0\n" % (catlist[ci].ra, catlist[ci].dec, catlist[cj].ra, catlist[cj].dec))

    #future project: if not enough, go to the secondary offsets

    return (primarymatchs, primarymatchc, mpa)

############################################


def writetextfile(filename, objlist):
    out = open(filename,'w')
    for ob in objlist:
        out.write("%11.7f %11.7f %5.2f\n" % (ob.ra, ob.dec, ob.mag))
    out.close()


def writeregionfile(filename, objlist, color="green",sys=''):
    if sys == '': sys = 'wcs'
    out = open(filename,'w')
    i = -1
    out.write('# Region file format: DS9 version 4.0\nglobal color='+color+' font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    if sys == 'wcs':
        out.write('fk5\n')
        for ob in objlist:
            i += 1
            out.write("point(%.7f,%.7f) # point=boxcircle text={%i}\n" % (ob.ra, ob.dec, i))
    if sys == 'img':
        out.write('image\n')
        for ob in objlist:
            i += 1
            out.write("point(%.3f,%.3f) # point=boxcircle text={%i}\n" % (ob.x, ob.y, i))
    out.close()


############################################

def autoastrometry(
        filename: str,
        pixel_scale: float = None,
        pa: float = None,
        inv: bool = False,
        uncpa: float = None,
        userra: float = None,
        userdec: float = None,
        max_ellip: float = 0.5,
        box_size: float = None,
        max_rad: float = None,
        tolerance: float = default_tolerance,
        catalog: str = "",
        no_solve: bool = False,
        overwrite: bool = True,
        outfile: str = "",
        output_dir: str = base_output_dir,
        saturation: float = default_saturation,
        no_rot: bool = False,
        min_fwhm: float = defaultminfwhm,
        max_fwhm: float = defaultmaxfwhm,
) -> (int, ):
    """

    Parameters
    ----------
    filename: Path of file
    pixel_scale: The pixel scale in arcsec/pix.  Must be within ~1%. By default: ???
    pa: The position angle in degrees.  Not usually needed.
    uncpa: Uncertainty of the position angle (degrees)
    inv: Reverse(=positive) parity.
    userra: RA in deg
    userdec: Dec in deg
    max_ellip: Maximum elliptical something?
    box_size: Half-width of box for reference catalog query (arcsec)
    max_rad: Maximum distance to look for star pairs.
    tolerance: Amount of slack allowed in match determination
    catalog: Catalog to use (ub2, tmc, sdss, or file)
    no_solve: Do not attempt to solve astrometry; just write catalog.
    overwrite: Overwrite output files
    outfile: Output file
    output_dir: Directory for output file
    saturation: Saturation level; do not use stars exceeding.
    no_rot: Some kind of bool
    min_fwhm: Minimum fwhm
    max_fwhm: Maximum fwhm

    Returns
    -------

    """

    # Get some basic info from the header  
    try:
        fits = af.open(filename)
        fits.verify('silentfix')
    except KeyboardInterrupt:
        logger.debug('Error opening', filename)
        if os.path.isfile(filename)==False:
            logger.debug('File does not exist.')
        raise AstrometryException

    sciext = 0

    for nh in range(len(fits)): # iterate through extensions to find the data
        try:
            if len(fits[nh].data) > 1:
                sciext = nh
                break
        except:
            pass

    h = fits[sciext].header

    #h = fits[0].header  #ideally check for primary extension, or even iterate
    #fits.close()
    sfilename = filename

    if np.logical_and(pixel_scale is not None, pa is None):
        pa = 0

    # Check for old-style WCS header
    if pixel_scale is None:
        #   hkeys = h.ascardlist().keys()
        hkeys = list(h.keys())
        oldwcstype = False
        for hkey in hkeys:
            if hkey == 'CDELT1' or hkey == 'CDELT2':
                oldwcstype = True
        if oldwcstype:
            key = 'CDELT1'
            cdelt1 = h[key]
            key = 'CDELT2'
            cdelt2 = h[key]

            try:
                crot = 0
                key = 'CROTA1'
                crot = h[key]
                key = 'CROTA2'
                crot = h[key]
            except KeyError:
                pass

            if sqrt(cdelt1**2 + cdelt2**2) < 0.1:  # some images use CDELT to indicate nonstandard things
                h['CD1_1'] = cdelt1 * cos(crot*math.pi/180.)
                h['CD1_2'] = -cdelt2 * sin(crot*math.pi/180.)
                h['CD2_1'] = cdelt1 * sin(crot*math.pi/180.)
                h['CD2_2'] = cdelt2 * cos(crot*math.pi/180.)

    if np.logical_and(pixel_scale is not None, pa is not None):
        # Create WCS header information if pixel scale is specified
        parad = pa * math.pi / 180.
        pxscaledeg = pixel_scale / 3600.
        if inv > 0:
            parity = -1
        else:
            parity = 1

        if userra is not None:
            ra = userra
        else:
            ra = rasex2deg(h['CRVAL1'])

        if userdec is not None:
            dec = userdec
        else:
            dec = decsex2deg(h['CRVAL2'])

        try:
            epoch = float(h.get('EPOCH', 2000))
        except:
            #warn?
            epoch = 2000.

        try:
            equinox = float(h.get('EQUINOX', epoch)) #If RA and DEC are not J2000 then convert
        except:
            # warn??
            equinox = 2000. # could be 'J2000'; try to strip off first character?

        if abs(equinox-2000) > 0.5:
            logger.debug('Converting equinox from', equinox, 'to J2000')
            j2000 = ephem.Equatorial(ephem.Equatorial(str(ra/15), str(dec), epoch=str(equinox)),epoch=ephem.J2000)
            [ra, dec] = [rasex2deg(j2000.ra), decsex2deg(j2000.dec)]

        h["CD1_1"] = pxscaledeg * cos(parad)*parity
        h["CD1_2"] = pxscaledeg * sin(parad)
        h["CD2_1"] = -pxscaledeg * sin(parad)*parity
        h["CD2_2"] = pxscaledeg * cos(parad)
        h["CRPIX1"] = h['NAXIS1']/2
        h["CRPIX2"] = h['NAXIS2']/2
        h["CRVAL1"] = ra
        h["CRVAL2"] = dec
        h["CTYPE1"] = "RA---TAN"
        h["CTYPE2"] = "DEC--TAN"
        #h.update("EPOCH",2000.0)  does it matter?
        h["EQUINOX"] = 2000.0
        #print ra, dec

        try:
            os.remove('temp.fits')
        except FileNotFoundError:
            pass

        fits[sciext].header = h
        fits.writeto(os.path.join(output_dir, 'temp.fits'), output_verify='silentfix') #,clobber=True
        fits.close()
        fits = af.open(os.path.join(output_dir, 'temp.fits'))
        h = fits[sciext].header
        sfilename = os.path.join(output_dir, 'temp.fits')

    #Read the header info from the file.
    try:
        # no longer drawing RA and DEC from here.
        key = 'NAXIS1'
        nxpix = h[key]
        key = 'NAXIS2'
        nypix = h[key]
    except:
        logger.debug('Cannot find necessary WCS header keyword', key)
        sys.exit(1)
    try:
        key = 'CRVAL1'
        cra = float(h[key])
        key = 'CRVAL2'
        cdec = float(h[key])

        key = 'CRPIX1'
        crpix1 = float(h[key])
        key = 'CRPIX2'
        crpix2 = float(h[key])

        key = 'CD1_1'
        cd11 = float(h[key])
        key = 'CD2_2'
        cd22 = float(h[key])
        key = 'CD1_2'
        cd12 = float(h[key]) # deg / pix
        key = 'CD2_1'
        cd21 = float(h[key])

        equinox = float(h.get('EQUINOX', 2000.))
        if abs(equinox-2000.) > 0.2: logger.debug('Warning: EQUINOX is not 2000.0')
    except:
        if pixel_scale == -1:
            logger.debug('Cannot find necessary WCS header keyword', key)
            logger.debug('Must specify pixel scale (-px VAL) or provide provisional basic WCS info via CD matrix.')
            #Some images might use CROT parameters, could try to be compatible with this too...?
            sys.exit(1)

    # Wipe nonstandard fits info from the header (otherwise this will confuse verification)
    hkeys = list(h.keys())
    ctypechange = 0
    irafkeys = []
    highkeys = []
    oldkeys = []
    distortionkeys = []
    for hkey in hkeys:
        if hkey=='RADECSYS' or hkey == 'WCSDIM' or hkey.find('WAT')==0 or hkey.find('LTV')>=0 or hkey.find('LTM')==0:
            del h[hkey]
            irafkeys.append(hkey)
        if hkey.find('CO1_') == 0 or hkey.find('CO2_') == 0 or hkey.find('PV1_') == 0 or hkey.find('PV2_') == 0 or hkey.find('PC00') == 0:
            del h[hkey]
            highkeys.append(hkey)
        if hkey.find('CDELT1') == 0 or hkey.find('CDELT2') == 0 or hkey.find('CROTA1') == 0 or hkey.find('CROTA2') == 0:
            del h[hkey]
            oldkeys.append(hkey)
        if hkey.find('A_') == 0 or hkey.find('B_') == 0 or hkey.find('AP_') == 0 or hkey.find('BP_') == 0:
            del h[hkey]
            distortionkeys.append(hkey)

    if h['CTYPE1'] != 'RA---TAN':
        logger.debug('Changing CTYPE1 from', h['CTYPE1'], 'to', "RA---TAN")
        h["CTYPE1"] = "RA---TAN"
        ctypechange = 1

    if h['CTYPE2'] != 'DEC--TAN':
        if ctypechange:
            logger.debug('Changing CTYPE2 from', h['CTYPE2'], 'to', "DEC--TAN")
        h["CTYPE2"] = "DEC--TAN"
        ctypechange = 1

    wcskeycheck = ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_2', 'CD2_1', 'EQUINOX', 'EPOCH']
    headerformatchange = 0

    for w in wcskeycheck:
        if type(w) == type('0'):
            try:
                h[w] = float(h[w])
                headerformatchange = 1
            except:
                pass

    if len(irafkeys) > 0:
        logger.debug('Removed nonstandard WCS keywords: ')
        for key in irafkeys:
            logger.debug(key)
    if len(highkeys) > 0:
        logger.debug('Removed higher-order WCS keywords: ')
        for key in highkeys:
            logger.debug(key)
    if len(oldkeys) > 0:
        logger.debug('Removed old-style WCS keywords: ')
        for key in oldkeys:
            logger.debug(key)
    if len(distortionkeys) > 0:
        logger.debug('Removed distortion WCS keywords: ')
        for key in distortionkeys:
            logger.debug(key)

    if len(highkeys)+len(distortionkeys)+ctypechange+headerformatchange > 0:
        #Rewrite and reload the image if the header was modified in a significant way so sextractor sees the same thing that we do.
        try:
            os.remove(os.path.join(output_dir, 'temp.fits'))
        except FileNotFoundError:
            pass
        fits[sciext].header = h
        fits.writeto(os.path.join(output_dir, 'temp.fits'), output_verify='silentfix') #,clobber=True
        fits.close()
        fits = af.open(os.path.join(output_dir, 'temp.fits'))
        h = fits[sciext].header
        sfilename = os.path.join(output_dir, 'temp.fits')

    # Get image info from header (even if we put it there in the first place)
    if cd11 * cd22 < 0 or cd12 * cd21 > 0:
        parity = -1
    else:
        parity = 1

    xscale = sqrt(cd11**2 + cd21**2)
    yscale = sqrt(cd12**2 + cd22**2)
    initpa = -parity * np.arctan2(cd21 * yscale, cd22 * xscale) * 180 / math.pi
    xscale = abs(xscale)
    yscale = abs(yscale)
    fieldwidth = max(xscale * nxpix, yscale * nypix) * 3600.
    area_sqdeg = xscale * nxpix * yscale * nypix
    area_sqmin = area_sqdeg * 3600.
    centerx = nxpix/2
    centery = nypix/2
    centerdx = centerx - crpix1
    centerdy = centery - crpix2
    centerra  = cra - centerdx*xscale*cos(initpa*math.pi/180.) + centerdy*yscale*sin(initpa*math.pi/180.)
    centerdec = cdec + parity*centerdx*xscale*sin(-initpa*math.pi/180.) + centerdy*yscale*cos(initpa*math.pi/180.)

    # this has only been checked for a PA of zero.

    logger.debug('Initial WCS info:')
    logger.debug('   pixel scale:     x=%.4f"/pix,   y=%.4f"/pix' % (xscale*3600, yscale*3600))
    logger.debug('   position angle: PA=%.2f' % initpa)
    if parity == 1:
        logger.debug('   normal parity')
    if parity == -1:
        logger.debug('   inverse parity')
    logger.debug('   center:        RA=%10.6f, dec=%9.6f' % (centerra, centerdec))

    # Sextract stars to produce image star catalog

    goodsexlist = sextract(
        sfilename,
        nxpix,
        nypix,
        3,
        12,
        minfwhm=min_fwhm,
        maxfwhm=max_fwhm,
        maxellip=max_ellip,
        saturation=saturation,
        output_dir=output_dir
    )

    ngood = len(goodsexlist)

    if ngood < 4:

        err = f'Only {ngood} good stars were found in the image.  The image is too small or shallow, the ' \
              f'detection threshold is set too high, or stars and cosmic rays are being confused.'
        logger.error(err)
        writetextfile('det.init.txt', goodsexlist)
        writeregionfile('det.im.reg', goodsexlist, 'red', 'img')
        raise AstrometryException(err)

    density = len(goodsexlist) / area_sqmin
    logger.debug('Source density of %f4 /arcmin^2' % density)

    if no_solve is True:
        if catalog == '':
            catalog = 'det.ref.txt'
        writetextfile(catalog, goodsexlist)
        return

    # If no catalog specified, check availability of SDSS
    if catalog == '':
        trycats = ['sdss', 'ub2', 'tmc']
        for trycat in trycats:
            testqueryurl = f"http://tdc-www.harvard.edu/cgi-bin/scat?catalog={trycat}&ra={centerra}" \
                           f"&dec={centerdec}&system=J2000&rad=-90"
            check = urllib.request.urlopen(testqueryurl)
            checklines = check.readlines()
            check.close()
            if len(checklines) > 15:
                catalog = trycat
                logger.debug('Using catalog', catalog)
                break
        if catalog == '':
            err = 'No catalog is available.  Check your internet connection.'
            logger.error(err)
            raise AstrometryException(err)

    # Load in reference star catalog
    if box_size is None:
        box_size = fieldwidth

    catlist = getcatalog(catalog, centerra, centerdec, box_size)

    ncat = len(catlist)
    catdensity = ncat / (2 * box_size / 60.) ** 2
    logger.debug(f'{ncat} good catalog objects.')
    logger.debug(f'Source density of {catdensity} /arcmin^2')

    if ncat > 0 | ncat < 5:
        logger.warning(f'Only {ncat} catalog objects in the search zone.'
                       f'Increase the magnitude threshold or box size.')

    if ncat == 0:
        logger.error('No objects found in catalog.')
        logger.error('The web query failed, all stars were excluded by the FHWM clip, or the image')
        logger.error('is too small.  Check input parameters or your internet connection.')
        raise AstrometryException

    # If this image is actually shallower than reference catalog, trim the reference catalog down
    if ncat > 16 and catdensity > 3 * density:
        logger.debug('Image is shallow.  Trimming reference catalog...')
        while catdensity > 3 * density:
            catlist = catlist[0:int(len(catlist)*4/5)]
            ncat = len(catlist)
            catdensity = ncat / (2 * box_size / 60.) ** 2

    # If the image is way deeper than USNO, trim the image catalog down
    if ngood > 16 and density > 4 * catdensity and ngood > 8:
        logger.debug('Image is deep.  Trimming image catalog...')
        while density > 4 * catdensity and ngood > 8:
            goodsexlist = goodsexlist[0:int(len(goodsexlist)*4/5)]
            ngood = len(goodsexlist)
            density = ngood / area_sqmin

    # If too many objects, do some more trimming
    if ngood*ncat > 120*120*4:
        logger.debug('Image and/or catalog still too deep.  Trimming...')
        while ngood*ncat > 120*120*4:
            if density > catdensity:
                goodsexlist = goodsexlist[0:int(len(goodsexlist)*4/5)]
                ngood = len(goodsexlist)
                density = ngood / area_sqmin
            else:
                catlist = catlist[0:int(len(catlist)*4/5)]
                ncat = len(catlist)
                catdensity = ncat / (2 * box_size / 60.) ** 2

    # Remove fainter object in close pairs for both lists
    minsep = 3
    deletelist = []
    for i in range(len(goodsexlist)):
        for j in range(i+1, len(goodsexlist)):
            if i == j: continue
            dist = distance(goodsexlist[i], goodsexlist[j])
            if dist < minsep:
                if goodsexlist[i].mag > goodsexlist[j].mag:
                    deletelist.append(i)
                else:
                    deletelist.append(j)
    deletelist = unique(deletelist)
    deletelist.reverse()
    for d in deletelist:
        del goodsexlist[d]
    deletelist = []
    for i in range(len(catlist)):
        for j in range(i+1, len(catlist)):
            if i == j: continue
            dist = distance(catlist[i], catlist[j])
            if dist < minsep:
                if catlist[i].mag > catlist[j].mag:
                    deletelist.append(i)
                else:
                    deletelist.append(j)
    deletelist = unique(deletelist)
    deletelist.reverse()
    for d in deletelist:
        del catlist[d]

    writetextfile('det.init.txt', goodsexlist)
    writeregionfile('det.im.reg', goodsexlist, 'red', 'img')
    writetextfile('cat.txt', catlist)
    writeregionfile('cat.wcs.reg', catlist, 'green', 'wcs')

    # The catalogs have now been completed.


    # Now start getting into the actual astrometry.

    min_rad = 5.0
    #if (maxrad == -1): maxrad = 180
    if max_rad is None:
        max_rad = 60 * (15 / (math.pi * min(density, catdensity))) ** 0.5 # numcomp ~ 15 [look at 15 closest objects
        max_rad = max(max_rad, 60.0)                                 #               in sparse dataset]
        if max_rad == 60.0:
            min_rad = 10.0   # in theory could scale this up further to reduce #comparisons
        max_rad = min(max_rad, fieldwidth * 3. / 4)

        # note that density is per arcmin^2, while the radii are in arcsec, hence the conversion factor.

    circdensity = density * min([area_sqmin, (math.pi * (max_rad / 60.) ** 2 - math.pi * (min_rad / 60) ** 2)])
    circcatdensity = catdensity * (math.pi * (max_rad / 60.) ** 2 - math.pi * (min_rad / 60) ** 2)
    catperimage = catdensity * area_sqmin

    logger.debug('After trimming: ')
    logger.debug(f'{len(goodsexlist)} detected objects ({density:.2f}/arcmin^2, {circdensity:.1f}/searchzone)')
    logger.debug(f'{len(catlist)} catalog objects ({catdensity:.2f}/arcmin^2, {circcatdensity:.1f}/searchzone)')


    patolerance = defaultpatolerance
    expectfalsetrios = ngood * ncat * circdensity**2 * circcatdensity**2 * tolerance**2 * (patolerance/360.)**1

    overlap1 = 0.3 * min(1,catdensity/density) # fraction of stars in image that are also in catalog - a guess
    truematchesperstar = (circdensity * overlap1) # but how many matches >3 and >4?  some annoying binomial thing

    reqmatch = 3
    if expectfalsetrios > 30 and truematchesperstar >= 4:
        reqmatch = 4
    #should check that this will actually work for the catalog, too.
    if catperimage <= 6 or ngood <= 6:
        reqmatch = 2
    if catperimage <= 3 or ngood <= 3:
        reqmatch = 1
    #for an extremely small or shallow image

    logger.debug('Pair comparison search radius: %.2f"' % max_rad)
    logger.debug(f'Using reqmatch = {reqmatch}')
    (primarymatchs, primarymatchc, mpa) = distmatch(goodsexlist, catlist, max_rad, min_rad, tolerance, reqmatch, patolerance, uncpa)

    nmatch = len(primarymatchs)
    if nmatch == 0:
        err = (' No valid matches found!\n '
               'Possible issues: \n'
               '  - The specified pixel scale (or PA or parity) is incorrect.  Double-check the input value. \n'
               '  - The field is outside the catalog search region.  Check header RA/DEC or increase search radius. \n'
               ' - The routine is flooded by bad sources.  Specify or check the input seeing. \n'
               '  - The routine is flagging many real stars.  Check the input seeing. \n'
               ' You can display a list of detected/catalog sources using det.im.reg and cat.wcs.reg. \n'
               )
        logger.error(err)
        raise AstrometryException(err)

    if nmatch <= 2:
        logger.warning(f'Warning: only {nmatch} match(es).  Astrometry may be unreliable.')
        logger.warning('   Check the pixel scale and parity and consider re-running.')
        warning = 1

    #We now have the PA and a list of stars that are almost certain matches.
    offpa = median(mpa)  #get average PA from the excellent values
    stdevpa = stdev(mpa)

    skyoffpa = -parity*offpa # This appears to be necessary for the printed value to agree with our normal definition.

    #if abs(offpa) < 0.2: offpa = 0.0
    logger.debug('PA offset:')
    logger.debug('  dPA = %.3f  (unc. %.3f)' % (skyoffpa, stdevpa))

    if no_rot <= 0:
        # Rotate the image to the new, correct PA
        #  NOTE: when CRPIX don't match CRVAL this shifts the center and screws things up.
        #  I don't understand why they don't always match.  [[I think this was an equinox issue.
        #  should be solved now, but be alert for further problems.]]

        #Greisen et al.:
        #WCS_i = SUM[j] (CD_ij)(p_j - CRPIX_j)      i.e.
        # RA - CRVAL1 = CD1_1 (x - CRPIX1) + CD1_2 (y - CRPIX2)
        #dec - CRVAL2 = CD2_1 (x - CRPIX1) + CD2_2 (y - CRPIX2)   [times a projection scale...]

        #Rotate....
        rot = offpa * math.pi/180
        #...the image itself
        h["CD1_1"] = cos(rot)*cd11 - sin(rot)*cd21
        h["CD1_2"] = cos(rot)*cd12 - sin(rot)*cd22   # a parity issue may be involved here?
        h["CD2_1"] = sin(rot)*cd11 + cos(rot)*cd21
        h["CD2_2"] = sin(rot)*cd12 + cos(rot)*cd22
        #...the coordinates (so we don't have to resex)
        for i in range(len(goodsexlist)):  #do all of them, though this is not necessary
            #print "%11.7f %11.7f %5.2f " % (goodsexlist[i].ra, goodsexlist[i].dec, goodsexlist[i].mag),
            goodsexlist[i].rotate(offpa,cra,cdec)
            #print "%11.7f %11.7f %5.2f" % (goodsexlist[i].ra, goodsexlist[i].dec, goodsexlist[i].mag)
        #print
    else:
        rotwarn = ''
        if abs(skyoffpa) > 1.0: rotwarn = ' (WARNING: image appears rotated, may produce bad shift)'
        logger.debug('  Skipping rotation correction ')


    writetextfile('det.wcs.txt',goodsexlist)


    #print
    #imoffsets = []
    #imoffpas = []
    imraoffset = []
    imdecoffset = []
    for i in range(len(primarymatchs)):
        #print goodsexlist[primarymatchs[i]].ra, goodsexlist[primarymatchs[i]].dec, ' to ', catlist[primarymatchc[i]].ra, catlist[primarymatchc[i]].dec
        imraoffset.append(goodsexlist[primarymatchs[i]].ra - catlist[primarymatchc[i]].ra)
        imdecoffset.append(goodsexlist[primarymatchs[i]].dec - catlist[primarymatchc[i]].dec)
        #imoffsets.append(distance(goodsexlist[primarymatchs[i]], catlist[primarymatchc[i]]))
        #imoffpas.append( posangle(goodsexlist[primarymatchs[i]], catlist[primarymatchc[i]]))

    #for i in range(len(imraoffset)):
    #    print primarymatchs[i], ':', (imraoffset[i]-median(imraoffset))*3600*cos(cdec*math.pi/180), (imdecoffset[i]-median(imdecoffset))*3600


    raoffset = -median(imraoffset)
    decoffset = -median(imdecoffset)
    rastd = stdev(imraoffset)*cos(cdec*math.pi/180)  # all of these are in degrees
    decstd = stdev(imdecoffset)
    stdoffset = sqrt(rastd**2 + decstd**2)

    #finaloffset = median(imoffsets)
    #finalpa =    median(imoffpas)


    #PA is east of north.
    #raoffset  =  (1./3600) * finaloffset * sin(finalpa*math.pi/180) / cos(cdec*math.pi/180)
    #decoffset =  (1./3600) * finaloffset * cos(finalpa*math.pi/180)

    raoffsetarcsec = raoffset*3600*cos(cdec*math.pi/180)
    decoffsetarcsec = decoffset*3600
    totoffsetarcsec = (raoffsetarcsec**2 + decoffset**2)**0.5
    stdoffsetarcsec = stdoffset*3600

    logger.debug('Spatial offset:')
    #print '  %.2f" (PA = %.2f deg)' % (finaloffset, finalpa)
    logger.debug('  dra = %.2f",  ddec = %.2f"  (unc. %.3f")' % (raoffsetarcsec, decoffsetarcsec, stdoffsetarcsec))

    warning = 0
    if (stdoffset*3600 > 1.0):
        logger.debug('WARNING: poor solution - some matches may be bad.  Check pixel scale?')
        warning = 1

    h["CRVAL1"] = cra + raoffset
    h["CRVAL2"] = cdec + decoffset

    #h.update("ASTRMTCH", catalog)
    try:
        oldcat = h['ASTR_CAT']
        h["OLD_CAT"] = (oldcat, "Earlier reference catalog")
    except:
        pass
    h["ASTR_CAT"] = (catalog, "Reference catalog for autoastrometry")
    h["ASTR_UNC"] = (stdoffsetarcsec, "Astrometric scatter vs. catalog (arcsec)")
    h["ASTR_SPA"] = (stdevpa, "Measured uncertainty in PA (degrees)")
    h["ASTR_DPA"] = (skyoffpa, "Change in PA (degrees)")
    h["ASTR_OFF"] = (totoffsetarcsec, "Change in center position (arcsec)")
    h["ASTR_NUM"] = (len(primarymatchs), "Number of matches")

    #Write out a match list to allow doing a formal fit with WCStools.

    outmatch = open('match.list','w')
    for i in range(len(primarymatchs)):
        si = primarymatchs[i]
        ci = primarymatchc[i]
        outmatch.write("%s %s  %s %s\n" % (goodsexlist[si].x, goodsexlist[si].y, catlist[ci].ra, catlist[ci].dec))
        #goodsexlist[si].ra, goodsexlist[si].dec))
    outmatch.close()

    # Could repeat with scale adjustment


    # Could then go back to full good catalog and match all sources

    if overwrite:
        outfile = filename
    if outfile == '':
        slashpos = filename.rfind('/')
        dir = filename[0:slashpos+1]
        fil = filename[slashpos+1:]
        outfile = dir+'a'+fil # alternate behavior would always output to current directory
    try:
        os.remove(outfile)
    except:
        pass

    fits[sciext].header = h
    fits.writeto(outfile,output_verify='silentfix') #,clobber=True
    logger.info('Written to '+outfile)

    # Delete!

    if (warning == 0):
        try:
            #os.remove("temp.fits")
            pass
        except:
            pass

    #print scatpath+'imwcs '+outfile+' -u match.list -w'

    #imwcscommand = scatpath+'imwcs -h 500 -w -v -t 15 -u match.list -o '+'a'+outfile + ' ' + outfile
    #print imwcscommand
    #os.system(imwcscommand)   For some reason this is giving ridiculous results


    fits.close()

    return nmatch, skyoffpa, stdevpa, raoffsetarcsec, decoffsetarcsec, stdoffsetarcsec      #stdoffset*3600


# #####################################################
# def argument_help():
#     #          12345678901234567890123456789012345678901234567890123456789012345678901234567890
#     logger.debug("autoastrometry.py - A fast astrometric solver")
#     logger.debug("Command-line options:")
#     logger.debug("   Automatic WCS information:")
#     logger.debug("      -px pixelscale:  The pixel scale in arcsec/pix.  Must be within ~1%.")
#     logger.debug("      -pa PA:          The position angle in degrees.  Not usually needed.")
#     logger.debug("      -inv:            Reverse(=positive) parity.")
#     logger.debug("   If your image has provisional WCS already, you do not need this information.")
#     #print "   If your image has detailed WCS, you need to flush it for this code to work."
#     logger.debug("   Options:")
#     logger.debug("      -b  boxsize:     Half-width of box for reference catalog query (arcsec)")
#     logger.debug("      -s  seeing:      Approximate seeing in pixels for CR/star/galaxy ID'ing.")
#     logger.debug("      -x  saturation:  Saturation level; do not use stars exceeding.")
#     logger.debug("      -upa PAdiff:     Uncertainty of the position angle (degrees)")
#     logger.debug("   Additional options:")
#     logger.debug("      -c  catalog:     Catalog to use (ub2, tmc, sdss, or file: see --catalog)")
#     logger.debug("      -d  searchdist:  Maximum distance to look for star pairs.")
#     logger.debug("      -t  tolerance:   Amount of slack allowed in match determination")
#     logger.debug("      -o  output:      Specify output file (no argument overwrites input file)")
#     logger.debug("      -n:              Do not attempt to solve astrometry; just write catalog.")
#     logger.debug("   More help:")
#     logger.debug("      --examples:      Some examples for running at the command line.")
#     logger.debug("      --trouble:       Troubleshooting info if you can't get a field to solve.")
#     logger.debug("      --catalog:       More information on catalogs.")
#     logger.debug("      --algorithm:     Description of the program runtime steps.")
#     logger.debug("      --output:        Information on the output files.")
#     logger.debug("      --input:         Information on the input files.")


def algorithmhelp():
    logger.debug("Algorithm info:")
    logger.debug("  The code uses a combination of pair-distance matching and asterism matching to")
    logger.debug("  solve a field without knowing the position angle.  The full list of steps is:")
    logger.debug(" 1 - Extract all stars from the image.  Attempts to filter out cosmic rays,")
    logger.debug("     galaxies, bad columns, and so on.")
    logger.debug(" 2 - Query a catalog to get a list of known star positions.")
    logger.debug(" 3 - Trim the catalogs if necessary and try to optimize the search strategy")
    logger.debug("     based on the relative catalog densities and areas.")
    logger.debug(" 4 - Calculate the distance from every star to its neighbors, both in the")
    logger.debug("     image and the catalog.  Stars too far apart are ignored.")
    logger.debug("     If the distance between a star and a neighbor matches the distance between")
    logger.debug("     some catalog star and one of its neighbors, calculate that PA too and store")
    logger.debug("     as a potential match.")
    logger.debug(" 5 - Look for stars that seem to have a lot of distance matches with catalog")
    logger.debug("     stars.  If the PA's also match, store the list as an asterism match.")
    logger.debug(" 6 - Prune out asterisms whose matches seem least useful.")
    logger.debug(" 7 - Calculate analytically the PA shift and offset by averaging asterism")
    logger.debug("     centers.")
    logger.debug(" 8 - Update the header and write to disk.")


# def troublehelp():
#     #          12345678901234567890123456789012345678901234567890123456789012345678901234567890
#     logger.debug("Troubleshooting info:")
#     logger.debug("   Supplying the correct pixel scale (within 1%) and correct parity is critical")
#     logger.debug("   if the image does not already contain this information in the FITS header.")
#     logger.debug("   If you have difficulty solving a field correctly, double-check these values.")
#     logger.debug("   If still having trouble, try opening temp.fits and an archival image of the")
#     logger.debug("   field (from DSS, etc.) and loading the .reg files in DS9.  The problem might")
#     logger.debug("   be in the telescope pointing/header info (in this case, increase the boxsize)")
#     logger.debug("   or good matching stars may be thrown away or confused by artifacts (in this")
#     logger.debug("   case, specify a seeing value).  If the PA is known, restricting it can also")
#     logger.debug("   help (try -upa 0.5); by default all orientations are searched.")
#     logger.debug("   If still having issues, e-mail dperley@astro.berkeley.edu for help.")


def examplehelp():
    #         12345678901234567890123456789012345678901234567890123456789012345678901234567890
    logger.debug("Examples:")
    logger.debug(" For an image with provisional WCS header but possibly incorrect PA:")
    logger.debug("    autoastrometry image.fits   ")
    logger.debug(" For an image with provisional WCS, definitely correct PA:")
    logger.debug("    autoastrometry image.fits -upa 0.5")
    logger.debug(' For an image with no WCS (or bad WCS) and a pixel scale of 0.35"/pixel:')
    logger.debug("    autoastrometry image.fits -px 0.35")
    logger.debug(' For an image with no WCS, pixel scale of 0.35", and PA of 121 degrees:')
    logger.debug("    autoastrometry image.fits -px 0.35 -pa 121 -upa 0.5")
    logger.debug(" For an image with no WCS, pixel scale of 0.35, and positive parity:")
    logger.debug("    autoastrometry image.fits -px 0.35 -inv")
    logger.debug(" Use a catalog on disk instead of wcstools query:")
    logger.debug("    autoastrometry image.fits -c catalog.txt")
    logger.debug(" Widen the catalog search to 12x12 arcminutes if the pointing is bad:")
    logger.debug("    autoastrometry image.fits -b 720")
    logger.debug(" Specify seeing (7 pixels) to better exclude cosmic rays and galaxies:")
    logger.debug("    autoastrometry image.fits -s 7")
    logger.debug(" Combine lots of options for a maximally directed solution:")
    logger.debug("    autoastrometry image.fits -px 0.35 -pa 121 -upa 0.5 -inv -b 600 -s 7")
    logger.debug(" (Substitute 'autoastrometry' with 'python autoastrometry.py' if not aliased.)")


def outputhelp():
    logger.debug("Explanation of output files:")
    logger.debug(" DS9 region files -")
    logger.debug("   cat.wcs.reg        - WCS positions of all objects in the catalog.")
    logger.debug("   det.im.reg         - XY positions of all non-flagged detected objects.")
    logger.debug("   matchlines.im.reg  - spoke lines for matched asterisms from XY positions.")
    logger.debug("   matchlines.wcs.reg - spoke lines for matched asterisms from WCS positions.")
    logger.debug(" Text output -")
    logger.debug("   cat.txt            - objects in the catalog as a text file (RA, dec, mag)")
    logger.debug("   det.init.txt       - objects in the image as a text file ('RA', 'dec', mag)")
    logger.debug("   det.final.txt      - objects in the image as a text file (RA, dec, mag)")
    logger.debug("   match.list         - list of matched stars: (X Y RA dec)")
    logger.debug(" Image output -")
    logger.debug("   a[image.fits]      - output image with astrometry (if not specified with -o)")
    logger.debug("   temp.fits          - image with provisional WCS from your -px and -pa inputs")

# def inputhelp():
#     logger.debug("Sextractor input files:")
#     logger.debug("  You can modify the following sextractor inputs if you choose.  This should")
#     logger.debug("  rarely be necessary, if ever.")
#     logger.debug("     sex.config - Overall input file")
#     logger.debug("     sex.conv   - Convolution matrix for detection")
#     logger.debug("  These files are written to disk using internal defaults if not present.")


#Some other conceivable options for the future:
#noncircular PSF
#instrument defaults

######################################################################
def usage(args):
    (xdir,xname) = os.path.split(args[0])
    logger.error("Usage:  %s filename(s) [-px pixelscale -pa PA -inv -b boxsize -s seeing -upa PAunc]" % xname)
    logger.error("     or %s -help for instructions and more options." % xname)

######################################################################


def main(
        files: str | list,
        seeing: float = None,
        pixel_scale: float = None,
        pa: float = None,
        uncpa: float = None,
        inv: bool = False,
        userra: float = None,
        userdec: float = None,
        max_ellip: float = 0.5,
        box_size: float = None,
        max_rad: float = None,
        tolerance: float = default_tolerance,
        catalog: str = "",
        no_solve: bool = False,
        overwrite: bool = False,
        outfile: str = None,
        output_dir: str = base_output_dir,
        saturation: float = default_saturation,
        no_rot: bool = False,
):
    """Function based on 'autoastrometry.py' by Daniel Perley and Kishalay De.

    This function runs sextractor and automatically performs astrometry.

    Supplying the correct pixel scale (within 1%) and correct parity is critical
    if the image does not already contain this information in the FITS header.
    If you have difficulty solving a field correctly, double-check these values.
    If still having trouble, try opening temp.fits and an archival image of the
    field (from DSS, etc.) and loading the .reg files in DS9. The problem might
    be in the telescope pointing/header info (in this case, increase the boxsize)
    or good matching stars may be thrown away or confused by artifacts (in this
    case, specify a seeing value).  If the PA is known, restricting it can also
    help (try -upa 0.5); by default all orientations are searched.

    Catalog info:
    Leave the catalog field blank will use SDSS if available and USNO otherwise.
    The catalog query uses wcstools (tdc-www.harvard.edu/wcstools).  However, you
    can also use your own catalog file on disk if you prefer using -c [filename]
    The default format is a text file with the first three columns indicating
    ra, dec, magnitude.  However, you can change the order of the columns by
    adding, e.g.#:1,2,6 to the first line.  In this case, this would indicate that the RA is in the
    1st column, dec in the 2nd, and magnitude in the 6th.   The mag column can be
    omitted completely, although if the catalog is not the same depth as the
    image this may compromise the search results.

    Parameters
    ----------
    files: Files to reduce
    seeing: Approximate seeing in pixels for CR/star/galaxy ID'ing.
    pixel_scale: The pixel scale in arcsec/pix.  Must be within ~1%. By default: ???
    pa: The position angle in degrees.  Not usually needed.
    uncpa: Uncertainty of the position angle (degrees)
    inv: Reverse(=positive) parity.
    userra: RA in deg
    userdec: Dec in deg
    max_ellip: Maximum elliptical something?
    box_size: Half-width of box for reference catalog query (arcsec)
    max_rad: Maximum distance to look for star pairs.
    tolerance: Amount of slack allowed in match determination
    catalog: Catalog to use (ub2, tmc, sdss, or file)
    no_solve: Do not attempt to solve astrometry; just write catalog.
    overwrite: Overwrite output files
    outfile: Output file
    output_dir: Output directory
    saturation: Saturation level; do not use stars exceeding.
    no_rot: Some kind of bool

    Returns
    -------

    """

    if isinstance(files, str):
        files = [files]

    if len(files) == 0:
        err = 'No files selected!'
        logger.error(err)
        raise ValueError(err)

    if np.logical_and(overwrite, outfile is not None):
        err = f"An output file was specified ({outfile}), but the script was configured to overwrite the original file."
        logger.error(err)
        raise ValueError(err)

    # nickel : 0.371"/pix, PA=177.7
    # lick Geminicam:  A=0.697"/pix  B=0.672"/pix
    # Nov observations are all offset by a large amount in dec even after equinox
    # Keck: 0.135"/pix, chipPA=90+instPA

    # args = []

    # while i < len(args):
    #     arg = args[i]
    #     isarg = 0
    #
    #     if (arg.find("-norot") == 0): # some ambiguity with -n (nosolve)
    #         norot = 1
    #         isarg = 1
    #
    #     if (arg.find("-d") == 0) and (arg.find("-dec") != 0):
    #         max_rad = float(args[i + 1])
    #         i+=1
    #         isarg = 1
    #
    #     if (arg.find("-o") == 0):
    #         outfile = ''
    #         if len(args) > i+1:
    #             if (args[i+1]).find('-') != 0:
    #                 outfile = args[i+1] #output
    #                 i+=1
    #         if outfile == '':
    #             overwrite = 1
    #         isarg = 1
    #
    #     if (arg.find("-n") == 0 and arg.find("-norot") == -1):
    #         nosolve = 1
    #         if len(args) > i+1:
    #             if (args[i+1]).find('-') != 0:
    #                 catalog = args[i+1] #output
    #                 i+=1
    #         isarg = 1
    #     if not isarg:
    #         if len(files) > 0:
    #             files += ",%s" % arg
    #         else:
    #             files = arg
    #     i+=1

    if seeing is None:
        min_fwhm = defaultminfwhm #1.5
        max_fwhm = defaultmaxfwhm #40
    else:
        min_fwhm = 0.7 * seeing
        max_fwhm = 2. * seeing

    write_param_file()
    write_config_file(saturation=saturation)

    n_image = len(files)
    failures = []
    questionable = []
    multiinfo = []

    for filename in files:

        if len(files) > 1:
            logger.debug(f'Processing {filename}')

        if no_solve and catalog == '':
            catalog = filename+'.cat'

        fitinfo = autoastrometry(
            filename,
            pixel_scale=pixel_scale,
            pa=pa,
            inv=inv,
            uncpa=uncpa,
            min_fwhm=min_fwhm,
            max_fwhm=max_fwhm,
            max_ellip=max_ellip,
            box_size=box_size,
            max_rad=max_rad,
            userra=userra,
            userdec=userdec,
            tolerance=tolerance,
            catalog=catalog,
            no_solve=no_solve,
            overwrite=overwrite,
            outfile=outfile,
            output_dir=output_dir,
            saturation=saturation,
            no_rot=no_rot
        )

        if no_solve:
            continue

        # WTF?

        if isinstance(fitinfo, int):
            fitinfo = (0,0,0,0,0,0)

        multiinfo.append(fitinfo)

        if fitinfo[0] == 0:   #number of matches
            failures.append(filename)
        if fitinfo[5] > 2:    #stdev of offset
            questionable.append(filename)

    if np.logical_and(n_image > 1, no_solve is False):

        if len(failures) == 0 and len(questionable) == 0:
            logger.info('Successfully processed all images!')
        else:
            logger.warning(f'Finished processing all images, not all were successful.')

        if len(questionable) > 0:
            logger.warning('The following images solved but have questionable astrometry: \n')
            for f in questionable:
                logger.warning(f)
        if len(failures) > 0:
            logger.error('The following images failed to solve: \n')
            for f in failures:
                logger.error(f)

        logger.debug("%25s " %'Filename')
        logger.debug("%6s %8s (%6s)  %7s %7s (%6s)" % ('#match', 'dPA ', 'stdev', 'dRA', 'dDec', 'stdev'))
        for i in range(len(files)):
            info = multiinfo[i]
            logger.debug("%25s " % files[i])
            if info[0] > 0:
                logger.debug("%6d %8.3f (%6.3f)  %7.3f %7.3f (%6.3f)" % info)
            else:
                logger.debug("failed to solve")

    try:
        os.remove('temp.param')
    except:
        logger.debug('Could not remove temp.param for some reason')


######################################################################
# Running as executable
if __name__ == '__main__':
    main(*sys.argv)

######################################################################


# some possible future improvements:
# verification to assess low-confidence solutions
# full automatic retry mode (parity detection, etc.)
# dealing with unknown pixel scale
# run wcstools for distortion parameters
# merge catalog check with catalog search to save a query
# improve the CR rejection further... maybe think about recognizing elliptical "seeing"?


