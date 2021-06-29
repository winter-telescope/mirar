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


#Modified by Kishalay De (kde@astro.caltech.edu) for removing dependency on deprecated pyfits
#and making it compatible with astropy headers and python 3.6 (June 11, 2018)

sexpath = ''  # if "sex" works in any directory, leave blank

defaulttolerance = 0.01  # these defaults should generally not be altered.
defaultpatolerance = 1.4   
defaultminfwhm = 1.5
defaultmaxfwhm = 40

fastmatch = 1
showmatches = 0

import numpy
#import pyfits
from astropy.io import fits as af
import os, sys
import glob, math
from math import sin, cos, tan, asin, sqrt, pi
import time, datetime
import string
import urllib.request, urllib.parse, urllib.error
try:
   import ephem
except:
   pass

def writeparfile():
    params = '''X_IMAGE
Y_IMAGE
ALPHA_J2000
DELTA_J2000
MAG_AUTO
MAGERR_AUTO
ELLIPTICITY
FWHM_IMAGE
FLAGS'''
    pf = open('temp.param','w')
    pf.write(params)
    pf.close()


def writeconfigfile(satlevel=55000.):
    configs='''
#-------------------------------- Catalog ------------------------------------
 
CATALOG_NAME     temp.cat       # name of the output catalog
CATALOG_TYPE     ASCII_HEAD     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC
PARAMETERS_NAME  temp.param     # name of the file containing catalog contents
 
#------------------------------- Extraction ----------------------------------
 
DETECT_TYPE      CCD            # CCD (linear) or PHOTO (with gamma correction)
DETECT_MINAREA   5              # minimum number of pixels above threshold
DETECT_THRESH    3              # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH  3              # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
 
FILTER           Y              # apply filter for detection (Y or N)?
FILTER_NAME      sex.conv       # name of the file containing the filter
 
DEBLEND_NTHRESH  16             # Number of deblending sub-thresholds
DEBLEND_MINCONT  0.02           # Minimum contrast parameter for deblending
 
CLEAN            Y              # Clean spurious detections? (Y or N)?
CLEAN_PARAM      1.0            # Cleaning efficiency
 
MASK_TYPE        CORRECT        # type of detection MASKing: can be one of
                                # NONE, BLANK or CORRECT
 
#------------------------------ Photometry -----------------------------------
 
PHOT_APERTURES   5              # MAG_APER aperture diameter(s) in pixels
PHOT_AUTOPARAMS  2.5, 3.5       # MAG_AUTO parameters: <Kron_fact>,<min_radius>
PHOT_PETROPARAMS 2.0, 3.5       # MAG_PETRO parameters: <Petrosian_fact>,
                                # <min_radius>
 
 
 
MAG_ZEROPOINT    0.0            # magnitude zero-point
MAG_GAMMA        4.0            # gamma of emulsion (for photographic scans)
GAIN             0.0            # detector gain in e-/ADU
PIXEL_SCALE      1.0            # size of pixel in arcsec (0=use FITS WCS info)
 
#------------------------- Star/Galaxy Separation ----------------------------
 
SEEING_FWHM      1.2            # stellar FWHM in arcsec
STARNNW_NAME     default.nnw    # Neural-Network_Weight table filename
 
#------------------------------ Background -----------------------------------
 
BACK_SIZE        64             # Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE  3              # Background filter: <size> or <width>,<height>
 
BACKPHOTO_TYPE   GLOBAL         # can be GLOBAL or LOCAL
 
#------------------------------ Check Image ----------------------------------
 
CHECKIMAGE_TYPE  NONE           # can be NONE, BACKGROUND, BACKGROUND_RMS,
                                # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                # or APERTURES
CHECKIMAGE_NAME  check.fits     # Filename for the check-image
 
#--------------------- Memory (change with caution!) -------------------------
 
MEMORY_OBJSTACK  3000           # number of objects in stack
MEMORY_PIXSTACK  300000         # number of pixels in stack
MEMORY_BUFSIZE   1024           # number of lines in buffer
 
#----------------------------- Miscellaneous ---------------------------------
 
VERBOSE_TYPE     QUIET          # can be QUIET, NORMAL or FULL
WRITE_XML        N              # Write XML file (Y/N)?
XML_NAME         sex.xml        # Filename for XML output
'''
    #SATUR_LEVEL      '''+str(satlevel)+'''        # level (in ADUs) at which arises saturation
    pf = open('sex.config','w')
    pf.write(configs)
    pf.close()

    convol='''CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
'''
    if not os.path.exists('sex.conv'): 
        cf = open('sex.conv','w')
        cf.write(convol)
        cf.close()




class Obj:
    ra = 0.0
    dec = 0.0
    mag = 0.0
    
    ra_rad = 0.0
    dec_rad = 0.0
    
    def __init__(self, inra, indec, inmag):
        self.ra = inra
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
        x = (self.ra  - ra0 ) * rascale
        y = (self.dec - dec0)
        
        xrot = cosdpa * x - sindpa * y
        yrot = sindpa * x + cosdpa * y
        
        self.ra   = (xrot / rascale) + ra0
        self.dec  =  yrot + dec0 
        self.ra_rad  = self.ra  * math.pi/180
        self.dec_rad =  self.dec * math.pi/180

class SexObj(Obj):
    x = 0.
    y = 0.
    mag = 0.0
    magerr = 0.0
    ellip = 0.0
    fwhm = 0.0
    flag = 0
    
    def __init__(self, inline):
        inlinearg = inline.split()
    
        if len(inlinearg) < 8: return # maybe throw an error?
        self.x = float(inlinearg[0])
        self.y = float(inlinearg[1])
        self.ra = float(inlinearg[2])
        self.dec = float(inlinearg[3])
        self.mag = float(inlinearg[4])
        self.magerr = float(inlinearg[5])
        self.ellip = float(inlinearg[6])
        self.fwhm = float(inlinearg[7])
        if len(inlinearg) >= 9: self.flag = int(inlinearg[8])

        self.ra_rad = self.ra * math.pi/180
        self.dec_rad = self.dec * math.pi/180

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
    dra  = obj2.ra  - obj1.ra
    if dra > 180: dra = 360 - dra
    return 3600 * sqrt(ddec**2 + (cosdec*dra)**2)

#Calculate the (spherical) position angle between two objects.
def posangle(obj1, obj2):
    
    dra  = obj2.ra_rad - obj1.ra_rad
    pa_rad = numpy.arctan2(cos(obj1.dec_rad)*tan(obj2.dec_rad)-sin(obj1.dec_rad)*cos(dra), sin(dra));
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
   a = numpy.array(l)
   return numpy.median(a)

def stdev(l):
   a = numpy.array(l)
   return numpy.std(a)

def mode(l):
   if len(l) == 0: return
   s = numpy.array(sorted(l))
   d = s[1:] - s[:-1]
   nd = len(d)
   if   nd >= 32: g = nd/16
   elif nd >= 6: g = 2
   else:         g = 1
   #g = max(nd / 16,1)  #sensitive to clusters up to a little less than 1/16 of the data set
   minmean = d.sum()
   imean = nd / 2
   for i in range(nd):
       r = [int(max(i-g,0)),int(min(i+g,nd))]
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
    if len(dec) == 1: return float(decstr)
    sign=1
    if (decstr[0] == '-'): sign=-1
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

def sextract(sexfilename, nxpix, nypix, border=3, corner=12, minfwhm=1.5, maxfwhm=25, maxellip=0.5, saturation=-1):

    if maxellip == -1: maxellip = 0.5
    if saturation > 0: 
       sexsaturation = saturation
    else:
       sexsaturation = 1e10

    try:
       # Sextract the image !
       os.system(sexpath+"sex " + sexfilename + " -c sex.config -SATUR_LEVEL "+str(sexsaturation))
    except:
       print(' Error: Problem running sextractor')
       print(' Check that program is installed and runs at command line using ' + sexpath+'sex')
       sys.exit(1)

    # Read in the sextractor catalog
    try:
       cat = open("temp.cat",'r')
       catlines = cat.readlines()
       cat.close()
    except:
       print('Cannot load sextractor output file!')
       sys.exit(1)

    if len(catlines) == 0:
       print('Sextractor catalog is empty: try a different catalog?')
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
    while l < len(catlines)-1:
        l += 1
        if (len(catlines[l]) <= 1 or catlines[l][0] == '#'):
            continue
        
        iobj = SexObj(catlines[l]) #process the line into an object
        nsexinit += 1
        
        #Initial filtering
        if iobj.ellip > maxellip : continue
        if iobj.fwhm < minfwhm: continue
        if iobj.fwhm > maxfwhm: continue
        if iobj.x < minx: continue
        if iobj.y < miny: continue
        if iobj.x > maxx: continue
        if iobj.y > maxy: continue
        if iobj.x + iobj.y < corner: continue
        if iobj.x + (nypix-iobj.y) < corner: continue
        if (nxpix-iobj.x) < corner: continue
        if (nxpix-iobj.x) + (nypix-iobj.y) < corner: continue
        if saturation > 0:
           if iobj.flag > 0: continue  # this will likely overdo it for very deep fields.
        
        sexlist.append(iobj)
        xlist.append(iobj.x)
        ylist.append(iobj.y)
        fwhmlist.append(iobj.fwhm)
        elliplist.append(iobj.ellip)
        flaglist.append(iobj.flag)
        nsexpass += 1

    #print nsexinit, 'raw sextractor detections'
    #print nsexpass, 'pass initial critiera'

     # Remove detections along bad columns

    threshprob = 0.0001
    ctbadcol = 0
    for i in range(5):
        txp = 1.0
        xthresh = 1
        while txp > threshprob: 
          txp *= min((len(sexlist)*1.0/nxpix),0.8) # some strange way of estimating the threshold.
          xthresh += 1                          #what I really want is a general analytic expression for
        removelist = []                         #the 99.99% prob. threshold for value of n for >=n out 
        modex = mode(xlist)                     #of N total sources to land in the same bin (of NX total bins)
        for j in range(len(sexlist)):
           if (sexlist[j].x > modex-1) and (sexlist[j].x < modex+1):
             removelist.append(j)
        removelist.reverse()
        if len(removelist) > xthresh:
         #print removelist
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
    if ctbadcol > 0: print(' Removed ', ctbadcol, ' detections along bad columns.')

    
    # Remove galaxies and cosmic rays

    if len(fwhmlist) > 5:
       fwhmlist.sort()
       fwhm20 = fwhmlist[int(len(fwhmlist)/5)]
       fwhm25 = fwhmlist[int(len(fwhmlist)/4)]
       fwhm50 = fwhmlist[int(len(fwhmlist)/2)]     #percentile values
       fwhm75 = fwhmlist[int(len(fwhmlist)*3/4)]
       fwhmmode = mode(fwhmlist)
    else:
       fwhmmode = minfwhm
       fwhm20 = minfwhm
    #hifwhmlist = []
    #for f in fwhmlist:
    #   if f > fwhmmode*1.5: hifwhmlist.append(f)
    #fwhmhimode = mode(hifwhmlist)

 #   ellipmode = mode(elliplist)   # in theory, if this is large we could use a theta cut, too.  (THETA_IMAGE)
 #   ellipstdev = stdev(elliplist)
 #   elliptol = 1.0 #min(0.25, 2.5*ellipstdev)
             # this effectively disables the ellipticity filter.
    
    #print fwhmlist
    #print fwhmmode #commonly CR's have their own cluster which can be tigheter than the real cluster?

    #print fwhmlist
    #print fwhm25, fwhm50, fwhm75, fwhmmode, fwhmhimode

    #print elliplist
    #print ellipmode, elliptol
    
    #print fwhmmode, fwhm20, minfwhm

    # formerly a max, but occasionally a preponderance of long CR's could cause fwhmmode to be bigger than the stars
    refinedminfwhm = median([0.75*fwhmmode,0.9*fwhm20,minfwhm]) # if CR's are bigger and more common than stars, this is dangerous...
    print('Refined min FWHM:', refinedminfwhm, 'pix')
    #refinedmaxfwhm = 35
    
    #maxfwhm = min(1.5 * fwhm50, 1.3 * fwhm75)  #some hard-coded parameters
    #minfwhm = max(0.7 * fwhm25, 0.5 * fwhm50)
    #It would be better to calculate a mode somehow.
    
    #print 'Maximum FWHM: ', maxfwhm
    #print 'Minimum FWHM: ', minfwhm
    
    #Might also be good to screen for false detections created by bad columns/rows
    

    ngood = 0
    goodsexlist = []
    for sex in sexlist:
       #print sex.x, sex.y, sex.ra, sex.dec, sex.mag, sex.ellip, sex.fwhm, sex.flag
       if sex.fwhm > refinedminfwhm: # and sex.ellip < ellipmode +elliptol:
          goodsexlist.append(sex)
          #print ' o',
          ngood += 1
       #print    

    # Sort by magnitude
    goodsexlist.sort(key=magcomp)   

    #i = 0
    #for sex in goodsexlist:
    #      if i < 1000: print i, sex.mag, sex.fwhm, sex.ellip
    #      i += 1


    print(len(sexlist), 'objects detected in image ('+ str(len(sexlist)-len(goodsexlist)) +' discarded)')


    return goodsexlist


def getcatalog(catalog, ra, dec, boxsize, minmag=8.0, maxmag=-1, maxpm=60.):
    # Get catalog from USNO

    if maxmag == -1:
        maxmag = 999 #default (custom catalog)
        if catalog == 'ub2': maxmag = 21.0#19.5
        if catalog == 'sdss': maxmag = 22.0
        if catalog == 'tmc': maxmag = 20.0

    
    if (catalog =='ub2' or catalog=='sdss' or catalog=='tmc'):
        usercat = 0
        racolumn = 1
        deccolumn = 2
        magcolumn = 6
        if catalog=='tmc': magcolumn=3
        pmracolumn = 10
        pmdeccolumn = 11   
        queryurl = "http://tdc-www.harvard.edu/cgi-bin/scat?catalog=" + catalog +  "&ra=" + str(ra) + "&dec=" + str(dec) + "&system=J2000&rad=" + str(-boxsize) + "&sort=mag&epoch=2000.00000&nstar=6400"
        #print queryurl
        cat = urllib.request.urlopen(queryurl)
        catlines = cat.readlines()
        cat.close()
        if len(catlines) > 6400-20:
           print('WARNING: Reached maximum catalog query size.')
           print('         Gaps may be present in the catalog, leading to a poor solution or no solution.')
           print('         Decrease the search radius.')
    else:
        usercat = 1
        try:
           cat = open(catalog,'r')
           print('Reading user catalog ', catalog)
        except:
           print('Failed to open user catalog ', catalog)
           print('File not found or invalid online catalog.  Specify ub2, sdss, or tmc.')
           return []
        racolumn = 0
        deccolumn = 1   # defaults
        magcolumn = -1    #  (to override, specify in first line using format #:0,1,2)  
        catlines = cat.readlines()
        cat.close()
    #print '                ', maxmag, maxpm
    l = -1
    catlist = []
    fwhmlist = []

    while l < len(catlines)-1:
        l += 1
        inline = catlines[l].strip()
        if len(inline) <= 2: continue
        if inline[0:2] == '#:':
             inlinearg = inline[2:].split(',')
             racolumn = int(inlinearg[0])-1
             deccolumn = int(inlinearg[1])-1
             if len(inlinearg) > 2: magcolumn = int(inlinearg[2])-1
             continue

        if (inline[0] < ord('0') or inline[0] > ord('9')) and str(inline[0]) != '.': continue #this may be too overzealous about
        if (inline[1] < ord('0') or inline[1] > ord('9')) and str(inline[1]) != '.': continue # removing comments...

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
        if mag > maxmag: continue #don't believe anything this faint
        if mag < minmag: continue #ignore anything this bright
        if abs(pmra) > maxpm or abs(pmdec) > maxpm: continue
        #print ' OK',
        iobj = Obj(ra, dec, mag) #process the line into an object
        catlist.append(iobj)
        
    catlist.sort(key=magcomp)   

    #print
    return catlist




def distmatch(sexlist, catlist, maxrad=180, minrad=10, tolerance=0.010, reqmatch=3, patolerance=1.2,uncpa=-1):
    
    if tolerance <= 0:
       print('Tolerance cannot be negative!!!')
       tolerance = abs(tolerance)
    if reqmatch < 2:
       print('Warning: reqmatch >=3 suggested')
    if patolerance <= 0: 
       print('PA tolerance cannot be negative!!!')
       patolerance = abs(patolerance)
    if uncpa < 0: uncpa = 720

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
        #print i, ':',
        d = []
        p = []
        dj = []
        for j in range(len(sexlist)):
            if i == j: continue
            if abs(sexlist[i].dec - sexlist[j].dec) > maxrad: continue
            if rascale*abs(sexlist[i].ra - sexlist[j].ra) > maxrad: continue
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
        #print i, ':',
        d = []
        p = []
        dj = []
        for j in range(len(catlist)):
            if i == j: continue
            if abs(catlist[i].dec - catlist[j].dec) > maxrad: continue
            if rascale*abs(catlist[i].ra - catlist[j].ra) > maxrad: continue
            dist = quickdistance(catlist[i], catlist[j], rascale)
            if dist > minrad and dist < maxrad :
                #print "%.2f " % dist, '['+str(j).strip()+']  ',
                d.append(dist)
                dj.append(j)
        #print
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
    offset = []
    offpa = []
    nmatch = []
    
    primarymatchs = []
    primarymatchc = []

    for si in range(len(sexdists)):
        sexdistarr = sexdists[si]
        sexidarr = sexmatchids[si]
        if len(sexdistarr) < 2: continue
        for ci in range(len(catdists)):
            catdistarr = catdists[ci]
            catidarr = catmatchids[ci]
            if len(catdistarr) < 2: continue
            match = 0
            smatchin = []
            cmatchin = []
            for sj in range(len(sexdistarr)):
                sexdist = sexdistarr[sj]
                newmatch = 1
                for cj in range(len(catdistarr)):
                    catdist = catdistarr[cj]
                    if abs((sexdist/catdist)-1.0) < tolerance:
                        #print si, ci, sexidarr[sj], catidarr[cj]
                        #print catmatchids[ci]
                        #print catdists[ci]
                        #print sexdist, catdist
                        #print '----'
                        match += newmatch
                        newmatch = 0 #further matches before the next sj loop indicate degeneracies
                        smatchin.append(sexmatchids[si][sj])
                        cmatchin.append(catmatchids[ci][cj])
            if match >= reqmatch: 
                #if showallmatches == 1:
                #    print si, 'matches', ci, ' ('+ str(match)+ ' matches)'
                #    print '  ',si,'-->', smatchin, '   ',ci,'-->', cmatchin
                
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
                    
                if len(smatchin) < 2: continue
                    
                dpamode = mode(dpa)
                #print dpa
                #print '       ', dpa
                #print '        mode =', mode(dpa)
                
                #Remove deviant matches by PA
                for i in range(len(smatchin)-1,-1,-1):
                    if abs(dpa[i] - dpamode) > patolerance:
                        del smatchin[i]
                        del cmatchin[i]
                        del dpa[i]
                
                if len(smatchin) < 2: continue

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
        print('Found no potential matches of any sort (including pairs).')
        print('The algorithm is probably not finding enough real stars to solve the field.  Check seeing.')
        #print 'This is an unusual error.  Check that a parameter is not set to a highly nonstandard value?'
        return [], [], []
    

    # Kill the bad matches
    rejects = 0
    
    #Get rid of matches that don't pass the reqmatch cut
    #if nmatches > 10 and max(nmatch) >= reqmatch:
    for i in range(len(primarymatchs)-1,-1,-1):
        #if len(smatch[i]) < reqmatch:
        if nmatch[i] < reqmatch:
            del mpa[i]
            del primarymatchs[i]
            del primarymatchc[i]
            del smatch[i]
            del cmatch[i]
            del nmatch[i]
            #rejects += 1  no longer a "reject"

    if len(smatch) < 1:
        print('Found no matching clusters of reqmatch =', reqmatch)
        return [], [], []


    #If we still have lots of matches, get rid of those with the minimum number of submatches
    #(that is, increase reqmatch by 1)
    minmatch = min(nmatch)
    countnotmin = 0
    for n in nmatch:
       if n > minmatch: countnotmin += 1
    if len(nmatch) > 16 and countnotmin > 3:
        print('Too many matches: increasing reqmatch to', reqmatch+1)
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
    print('Found',nmatches,'candidate matches.')


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
        if len(primarymatchs) == 0: break

        for i in range(len(primarymatchs)):
            for j in range(len(primarymatchs)):
                if i == j: continue
                si = primarymatchs[i]
                ci = primarymatchc[i]
                sj = primarymatchs[j]
                cj = primarymatchc[j]
    
                sexdistij = distance(sexlist[si], sexlist[sj])
                catdistij = distance(catlist[ci], catlist[cj])
                #print i, j, sexdistij, catdistij
    
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

    print('Rejected', rejects, 'bad matches.')
    nmatches = len(primarymatchs)
    print('Found', nmatches, 'good matches.')

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
           print('Refined pixel scale measurement: %.4f"/pix (+/- %.4f)' % (pixelscale, pixelscalestd))
        else:
           print('Refined pixel scale measurement: %.4f"/pix' % pixelscale)



    for i in range(len(primarymatchs)):
        si = primarymatchs[i]
        ci = primarymatchc[i]
        print('%3i' % si, 'matches', '%3i' % ci, ' (dPA =%7.3f)' % mpa[i], end=' ')
        if showmatches:
           print()
           if len(smatch[i]) < 16:
              print('  ', si, '-->', smatch[i], end=' ') 
              if len(smatch[i]) >= 7: print()
              print('  ', ci, '-->', cmatch[i])
           else:
              print('  ', si, '-->', smatch[i][0:10], '+', len(smatch[i])-10, 'more')
              print('  ', ci, '-->', cmatch[i][0:10], '+')#, len(cmatch[i])-10, ' more'
           if i+1 >= 10 and len(primarymatchs)-10 > 0: 
              print((len(primarymatchs)-10), 'additional matches not shown.')
              break
        else:
           print(':', str(len(smatch[i])).strip(), 'rays')
    
    out = open('matchlines.im.reg','w')
    i = -1
    color='red'
    out.write('# Region file format: DS9 version 4.0\nglobal color='+color+' font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    out.write('image\n')
    for i in range(len(primarymatchs)):
        si = primarymatchs[i]
        for j in range(len(smatch[i])):
           sj = smatch[i][j] 
           out.write("line(%.3f,%.3f,%.3f,%.3f) # line=0 0\n" % (sexlist[si].x, sexlist[si].y, sexlist[sj].x, sexlist[sj].y))
    out.close()

    out = open('matchlines.wcs.reg','w')
    i = -1
    color='green'
    out.write('# Region file format: DS9 version 4.0\nglobal color='+color+' font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    out.write('fk5\n')
    for i in range(len(primarymatchs)):
        ci = primarymatchc[i]
        for j in range(len(smatch[i])):
           cj = cmatch[i][j]
           out.write("line(%.5f,%.5f,%.5f,%.5f) # line=0 0\n" % (catlist[ci].ra, catlist[ci].dec, catlist[cj].ra, catlist[cj].dec))
    out.close()




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

def autoastrometry(filename,pixelscale=-1,pa=-999,inv=0,uncpa=-1,userra=-999, userdec=-999, minfwhm=1.5,maxfwhm=20,maxellip=0.5,boxsize=-1,maxrad=-1,tolerance=0.010,catalog='',nosolve=0,
overwrite=False, outfile='', saturation=-1, quiet=False, norot=0):
  
    # Get some basic info from the header  
    try: 
       fits = af.open(filename)
       fits.verify('silentfix')
    except:
       print('Error opening', filename)
       if os.path.isfile(filename)==False: print('File does not exist.')
       return -1

    sciext = 0
    for nh in range(len(fits)): # iterate through extensions to find the data
       try:
          if len(fits[nh].data) > 1 : 
            sciext = nh
            break
       except:
          pass
    h = fits[sciext].header

    #h = fits[0].header  #ideally check for primary extension, or even iterate
    #fits.close()
    sfilename = filename
    
    if pixelscale > 0 and pa == -999: pa = 0

    # Check for old-style WCS header
    if pixelscale < 0:
    #   hkeys = h.ascardlist().keys()
       hkeys = list(h.keys())
       oldwcstype = 0
       for hkey in hkeys:
          if hkey=='CDELT1' or hkey=='CDELT2': oldwcstype=1
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
          except:
             pass
          if sqrt(cdelt1**2 + cdelt2**2) < 0.1:  # some images use CDELT to indicate nonstandard things
             h['CD1_1'] = cdelt1 * cos(crot*math.pi/180.)
             h['CD1_2'] = -cdelt2 * sin(crot*math.pi/180.)
             h['CD2_1'] = cdelt1 * sin(crot*math.pi/180.)
             h['CD2_2'] = cdelt2 * cos(crot*math.pi/180.)


    if pixelscale > 0 and pa > -360:
       # Create WCS header information if pixel scale is specified
       parad = pa * math.pi / 180.
       pxscaledeg = pixelscale / 3600.
       if inv > 0: 
          parity = -1
       else:
          parity = 1
       if userra >= 0. and userra < 360.:
          ra = userra
       else:
          try:
             ra = rasex2deg(h['RA'])
          except:
             ra = rasex2deg(h['CRVAL1'])
       if userdec >= -90. and userdec <= 90.:
          dec = userdec
       else:
          try:
             dec = decsex2deg(h['DEC'])
          except:
             dec = decsex2deg(h['CRVAL2'])

       try:
          epoch = float(h.get('EPOCH', 2000))
       except:
          epoch = 2000.
       try:
          equinox = float(h.get('EQUINOX', epoch)) #If RA and DEC are not J2000 then convert
       except:
          equinox = 2000. # could be 'J2000'; try to strip off first character?

       if abs(equinox-2000) > 0.5:
           print('Converting equinox from', equinox, 'to J2000')
           try:
              j2000 = ephem.Equatorial(ephem.Equatorial(str(ra/15), str(dec), epoch=str(equinox)),epoch=ephem.J2000)
              [ra, dec] = [rasex2deg(j2000.ra), decsex2deg(j2000.dec)]
           except:
              print('PyEphem is not installed but is required to precess this image.')
              return -1
       h["CD1_1"] =  pxscaledeg * cos(parad)*parity
       h["CD1_2"] = pxscaledeg * sin(parad)
       h["CD2_1"] = -pxscaledeg * sin(parad)*parity
       h["CD2_2"] =  pxscaledeg * cos(parad)
       h["CRPIX1"] = h['NAXIS1']/2
       h["CRPIX2"] = h['NAXIS2']/2
       h["CRVAL1"] = ra
       h["CRVAL2"] = dec
       h["CTYPE1"] = "RA---TAN"
       h["CTYPE2"] = "DEC--TAN"
       #h.update("EPOCH",2000.0)  does it matter?
       h["EQUINOX"] = 2000.0
       #print ra, dec
       
       if os.path.isfile('temp.fits'): os.remove('temp.fits')
       fits[sciext].header = h
       fits.writeto('temp.fits',output_verify='silentfix') #,clobber=True
       fits.close()
       fits = af.open('temp.fits')
       h = fits[sciext].header
       sfilename = 'temp.fits'


    #Read the header info from the file.
    try:
        # no longer drawing RA and DEC from here.
        key = 'NAXIS1'
        nxpix = h[key]
        key = 'NAXIS2'
        nypix = h[key]
    except:
        print('Cannot find necessary WCS header keyword', key)
        sys.exit(1)
    try:
        key = 'CRVAL1'
        cra =  float(h[key])
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
        if abs(equinox-2000.) > 0.2: print('Warning: EQUINOX is not 2000.0')
    except:
        if pixelscale == -1:
            print('Cannot find necessary WCS header keyword', key)
            print('Must specify pixel scale (-px VAL) or provide provisional basic WCS info via CD matrix.')
            #Some images might use CROT parameters, could try to be compatible with this too...?
            sys.exit(1)

    # Wipe nonstandard fits info from the header (otherwise this will confuse verification)
    #hkeys = h.ascardlist().keys()
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
       if hkey.find('CO1_')==0 or hkey.find('CO2_')==0 or hkey.find('PV1_')==0 or hkey.find('PV2_')==0 or hkey.find('PC00')==0:
          del h[hkey]
          highkeys.append(hkey)
       if hkey.find('CDELT1')==0 or hkey.find('CDELT2')==0 or hkey.find('CROTA1')==0 or hkey.find('CROTA2')==0:
          del h[hkey]
          oldkeys.append(hkey)
       if hkey.find('A_')==0 or hkey.find('B_')==0 or hkey.find('AP_')==0 or hkey.find('BP_')==0 :
          del h[hkey]
          distortionkeys.append(hkey)
    if h['CTYPE1'] != 'RA---TAN': 
       if quiet==False: print('Changing CTYPE1 from',h['CTYPE1'],'to',"RA---TAN")
       h["CTYPE1"] = "RA---TAN"
       ctypechange = 1
    if h['CTYPE2'] != 'DEC--TAN': 
       if ctypechange: print('Changing CTYPE2 from',h['CTYPE2'],'to',"DEC--TAN")
       h["CTYPE2"] = "DEC--TAN"
       ctypechange = 1
    wcskeycheck = ['CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_2','CD2_1','EQUINOX','EPOCH']
    headerformatchange = 0
    for w in wcskeycheck:
       if type(w)==type('0'): 
          try:
             h[w] = float(h[w])
             headerformatchange = 1
          except:
             pass
    if quiet == False:
       if len(irafkeys) > 0:
          print('Removed nonstandard WCS keywords: ', end=' ')
          for key in irafkeys: print(key, end=' ')
          print()
       if len(highkeys) > 0:
          print('Removed higher-order WCS keywords: ', end=' ')
          for key in highkeys: print(key, end=' ')
          print()
       if len(oldkeys) > 0:
          print('Removed old-style WCS keywords: ', end=' ')
          for key in oldkeys: print(key, end=' ')
          print()
       if len(distortionkeys) > 0:
          print('Removed distortion WCS keywords: ', end=' ')
          for key in distortionkeys: print(key, end=' ')
          print()    
    if len(highkeys)+len(distortionkeys)+ctypechange+headerformatchange > 0:
       #Rewrite and reload the image if the header was modified in a significant way so sextractor sees the same thing that we do.
       if os.path.isfile('temp.fits'): os.remove('temp.fits')
       fits[sciext].header = h
       fits.writeto('temp.fits',output_verify='silentfix') #,clobber=True
       fits.close()
       fits = af.open('temp.fits')
       h = fits[sciext].header
       sfilename = 'temp.fits'

    # Get image info from header (even if we put it there in the first place)
    if cd11 * cd22 < 0 or cd12 * cd21 > 0:
       parity = -1
    else:
       parity = 1
    xscale = sqrt(cd11**2 + cd21**2)
    yscale = sqrt(cd12**2 + cd22**2)
    initpa = -parity * numpy.arctan2(cd21 * yscale, cd22 * xscale) * 180 / math.pi
    xscale = abs(xscale)
    yscale = abs(yscale)
    fieldwidth = max(xscale * nxpix, yscale * nypix) * 3600.
    area_sqdeg = xscale * nxpix * yscale * nypix
    area_sqmin = area_sqdeg * 3600. 
    area_sqsec = area_sqmin * 3600. 
    pixscale = sqrt(xscale*yscale) * 3600.
  
  
    centerx = nxpix/2
    centery = nypix/2
    centerdx = centerx - crpix1
    centerdy = centery - crpix2
    centerra  = cra  -        centerdx*xscale*cos(initpa*math.pi/180.) + centerdy*yscale*sin(initpa*math.pi/180.)
    centerdec = cdec + parity*centerdx*xscale*sin(-initpa*math.pi/180.) + centerdy*yscale*cos(initpa*math.pi/180.)
    # this has only been checked for a PA of zero.
  
  
    if quiet == False:
       print('Initial WCS info:')
       print('   pixel scale:     x=%.4f"/pix,   y=%.4f"/pix' % (xscale*3600, yscale*3600))
       print('   position angle: PA=%.2f' % initpa)
       if parity == 1: print('   normal parity')
       if parity ==-1: print('   inverse parity')
       print('   center:        RA=%10.6f, dec=%9.6f' % (centerra, centerdec))
  
    #catalog = 'ub2' #ub2
    #mincatmag = 8.0
    #maxcatmag = 19.0
    #maxpm = 60.       
    #boxsize = 200 #600
  
    # Sextract stars to produce image star catalog

    goodsexlist = sextract(sfilename, nxpix, nypix, 3, 12, minfwhm=minfwhm, maxfwhm=maxfwhm, maxellip=maxellip, saturation=saturation)
    
    ngood = len(goodsexlist)
    if ngood < 4:
       print('Only', ngood, 'good stars were found in the image.  The image is too small or shallow, the detection')
       print('threshold is set too high, or stars and cosmic rays are being confused.')
       #print '  If that does not work, the image may be too shallow for this method (or may be a flatfield, etc.).'
       writetextfile('det.init.txt', goodsexlist)
       writeregionfile('det.im.reg', goodsexlist, 'red', 'img')
       return -1
       #sys.exit(1)

    #print ngood , 'good stars found in frame.'
    density = len(goodsexlist) / area_sqmin
    print('Source density of %f4 /arcmin^2' % density)
    
    if nosolve == 1: 
       if catalog == '': catalog = 'det.ref.txt'
       writetextfile(catalog, goodsexlist)
       return

    # If no catalog specified, check availability of SDSS
    if catalog == '':
        trycats = ['sdss', 'ub2', 'tmc']
        for trycat in trycats:
            testqueryurl = "http://tdc-www.harvard.edu/cgi-bin/scat?catalog=" + trycat +  "&ra=" + str(centerra) + "&dec=" + str(centerdec) + "&system=J2000&rad=" + str(-90)
            #print testqueryurl
            check = urllib.request.urlopen(testqueryurl)
            checklines = check.readlines()
            check.close()
            if len(checklines) > 15:
                catalog = trycat
                print('Using catalog', catalog)
                break
        if (catalog == ''):
            print('No catalog is available.  Check your internet connection.')
            return -1


    # Load in reference star catalog
    if (boxsize == -1):
        boxsize = fieldwidth
    catlist = getcatalog(catalog, centerra, centerdec, boxsize)
    
    ncat = len(catlist)
    catdensity = ncat / (2*boxsize/60.)**2
    print(ncat, 'good catalog objects.')
    print('Source density of %f4 /arcmin^2' % catdensity)
    
    if ncat < 5 and ncat > 0:
       print('Only ', ncat, ' catalog objects in the search zone.  Increase the magnitude threshold or box size.')

    if ncat == 0 :
       print()
       print('No objects found in catalog.')
       print('The web query failed, all stars were excluded by the FHWM clip, or the image')
       print('is too small.  Check input parameters or your internet connection.')
       return -1
    

    #If this image is actually shallower than reference catalog, trim the reference catalog down
    if ncat > 16 and catdensity > 3 * density:
        print('Image is shallow.  Trimming reference catalog...')
        while catdensity > 3 * density:
            catlist = catlist[0:int(len(catlist)*4/5)]
            ncat = len(catlist)
            catdensity = ncat / (2*boxsize/60.)**2
        
    #If the image is way deeper than USNO, trim the image catalog down
    if ngood > 16 and density > 4 * catdensity and ngood > 8:
        print('Image is deep.  Trimming image catalog...')
        while density > 4 * catdensity and ngood > 8:
            goodsexlist = goodsexlist[0:int(len(goodsexlist)*4/5)]
            ngood = len(goodsexlist)
            density = ngood / area_sqmin

    #If too many objects, do some more trimming
    if ngood*ncat > 120*120*4:
        print('Image and/or catalog still too deep.  Trimming...')
        while ngood*ncat > 120*120*4:
            if density > catdensity: 
                goodsexlist = goodsexlist[0:int(len(goodsexlist)*4/5)]
                ngood = len(goodsexlist)
                density = ngood / area_sqmin
            else:
                catlist = catlist[0:int(len(catlist)*4/5)]
                ncat = len(catlist)
                catdensity = ncat / (2*boxsize/60.)**2


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
           
    minrad = 5.0
    #if (maxrad == -1): maxrad = 180
    if (maxrad == -1):                                            
        maxrad = 60*(15 / (math.pi * min(density,catdensity)) )**0.5 # numcomp ~ 15 [look at 15 closest objects
        maxrad = max(maxrad, 60.0)                                 #               in sparse dataset]
        if maxrad == 60.0:                           
             minrad = 10.0   # in theory could scale this up further to reduce #comparisons
        maxrad = min(maxrad, fieldwidth*3./4)                

    #note that density is per arcmin^2, while the radii are in arcsec, hence the conversion factor.
    circdensity =       density * min([area_sqmin, (math.pi*(maxrad/60.)**2 - math.pi*(minrad/60)**2)])
    circcatdensity = catdensity *                  (math.pi*(maxrad/60.)**2 - math.pi*(minrad/60)**2)
    catperimage    = catdensity * area_sqmin

    print('After trimming: ')
    print('   ', len(goodsexlist), 'detected objects (%.2f/arcmin^2, %.1f/searchzone)' % (density, circdensity))
    print('   ', len(catlist),     'catalog objects (%.2f/arcmin^2, %.1f/searchzone)' % (catdensity, circcatdensity))


    patolerance = defaultpatolerance
    expectfalsepairs = ngood * ncat * circdensity**1 * circcatdensity**1 * tolerance**1 * (patolerance/360.)**0
    expectfalsetrios = ngood * ncat * circdensity**2 * circcatdensity**2 * tolerance**2 * (patolerance/360.)**1
    expectfalsequads = ngood * ncat * circdensity**3 * circcatdensity**3 * tolerance**3 * (patolerance/360.)**2
    expectfalsequint = ngood * ncat * circdensity**4 * circcatdensity**4 * tolerance**4 * (patolerance/360.)**3

    overlap1 = 0.3 * min(1,catdensity/density) # fraction of stars in image that are also in catalog - a guess
    truematchesperstar = (circdensity * overlap1) # but how many matches >3 and >4?  some annoying binomial thing
    #print expectfalsepairs, expectfalsetrios, expectfalsequads, expectfalsequint

    reqmatch = 3
    if expectfalsetrios > 30 and truematchesperstar >= 4: reqmatch = 4   
       #should check that this will actually work for the catalog, too.
    if catperimage <= 6 or ngood <= 6: reqmatch = 2 
    if catperimage <= 3 or ngood <= 3: reqmatch = 1
        #for an extremely small or shallow image

    print('Pair comparison search radius: %.2f"'%maxrad)
    print('Using reqmatch =', reqmatch)
    (primarymatchs, primarymatchc, mpa) = distmatch(goodsexlist, catlist, maxrad, minrad, tolerance, reqmatch, patolerance, uncpa)

    nmatch = len(primarymatchs)
    if nmatch == 0:
        print(' No valid matches found!')
        if quiet == False:
           print(' Possible issues:')
           print('  - The specified pixel scale (or PA or parity) is incorrect.  Double-check the input value.')
           print('  - The field is outside the catalog search region.  Check header RA/DEC or increase search radius.')
           print('  - The routine is flooded by bad sources.  Specify or check the input seeing.')
           print('  - The routine is flagging many real stars.  Check the input seeing.')
           print(' You can display a list of detected/catalog sources using det.im.reg and cat.wcs.reg.')
        return -1
    if nmatch <= 2:
        print('Warning: only', nmatch, 'match(es).  Astrometry may be unreliable.')
        if quiet == False:
           print('   Check the pixel scale and parity and consider re-running.')
        warning = 1

    #We now have the PA and a list of stars that are almost certain matches.
    offpa = median(mpa)  #get average PA from the excellent values
    stdevpa = stdev(mpa)
    
    skyoffpa = -parity*offpa # This appears to be necessary for the printed value to agree with our normal definition.

    #if abs(offpa) < 0.2: offpa = 0.0
    print('PA offset:')
    print('  dPA = %.3f  (unc. %.3f)' % (skyoffpa, stdevpa))

    if norot <= 0:
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
       print('  Skipping rotation correction ')


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
    
    print('Spatial offset:')
    #print '  %.2f" (PA = %.2f deg)' % (finaloffset, finalpa)
    print('  dra = %.2f",  ddec = %.2f"  (unc. %.3f")' % (raoffsetarcsec, decoffsetarcsec, stdoffsetarcsec))
    
    warning = 0
    if (stdoffset*3600 > 1.0):
        print('WARNING: poor solution - some matches may be bad.  Check pixel scale?')
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
    
    if overwrite: outfile = filename
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
    print('Written to '+outfile)
    
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
    
    return (nmatch, skyoffpa, stdevpa, raoffsetarcsec, decoffsetarcsec, stdoffsetarcsec)      #stdoffset*3600


#####################################################
def help():
#          12345678901234567890123456789012345678901234567890123456789012345678901234567890
    print("autoastrometry.py - A fast astrometric solver")
    print("Command-line options:")
    print("   Automatic WCS information:")
    print("      -px pixelscale:  The pixel scale in arcsec/pix.  Must be within ~1%.")
    print("      -pa PA:          The position angle in degrees.  Not usually needed.")
    print("      -inv:            Reverse(=positive) parity.")
    print("   If your image has provisional WCS already, you do not need this information.")
    #print "   If your image has detailed WCS, you need to flush it for this code to work."
    print("   Options:")
    print("      -b  boxsize:     Half-width of box for reference catalog query (arcsec)")
    print("      -s  seeing:      Approximate seeing in pixels for CR/star/galaxy ID'ing.")
    print("      -x  saturation:  Saturation level; do not use stars exceeding.")
    print("      -upa PAdiff:     Uncertainty of the position angle (degrees)")
    print("   Additional options:")
    print("      -c  catalog:     Catalog to use (ub2, tmc, sdss, or file: see --catalog)")
    print("      -d  searchdist:  Maximum distance to look for star pairs.")
    print("      -t  tolerance:   Amount of slack allowed in match determination")
    print("      -o  output:      Specify output file (no argument overwrites input file)")
    print("      -n:              Do not attempt to solve astrometry; just write catalog.")
    print("   More help:")
    print("      --examples:      Some examples for running at the command line.")
    print("      --trouble:       Troubleshooting info if you can't get a field to solve.")
    print("      --catalog:       More information on catalogs.")
    print("      --algorithm:     Description of the program runtime steps.")
    print("      --output:        Information on the output files.")
    print("      --input:         Information on the input files.")

def algorithmhelp():
    print("Algorithm info:")
    print("  The code uses a combination of pair-distance matching and asterism matching to")
    print("  solve a field without knowing the position angle.  The full list of steps is:")
    print(" 1 - Extract all stars from the image.  Attempts to filter out cosmic rays,")
    print("     galaxies, bad columns, and so on.")
    print(" 2 - Query a catalog to get a list of known star positions.")
    print(" 3 - Trim the catalogs if necessary and try to optimize the search strategy")
    print("     based on the relative catalog densities and areas.")
    print(" 4 - Calculate the distance from every star to its neighbors, both in the")
    print("     image and the catalog.  Stars too far apart are ignored.")
    print("     If the distance between a star and a neighbor matches the distance between")
    print("     some catalog star and one of its neighbors, calculate that PA too and store")
    print("     as a potential match.")
    print(" 5 - Look for stars that seem to have a lot of distance matches with catalog")
    print("     stars.  If the PA's also match, store the list as an asterism match.")
    print(" 6 - Prune out asterisms whose matches seem least useful.")
    print(" 7 - Calculate analytically the PA shift and offset by averaging asterism")
    print("     centers.")
    print(" 8 - Update the header and write to disk.")


def troublehelp():
#          12345678901234567890123456789012345678901234567890123456789012345678901234567890
    print("Troubleshooting info:")
    print("   Supplying the correct pixel scale (within 1%) and correct parity is critical")
    print("   if the image does not already contain this information in the FITS header.")
    print("   If you have difficulty solving a field correctly, double-check these values.")
    print("   If still having trouble, try opening temp.fits and an archival image of the")
    print("   field (from DSS, etc.) and loading the .reg files in DS9.  The problem might")
    print("   be in the telescope pointing/header info (in this case, increase the boxsize)")
    print("   or good matching stars may be thrown away or confused by artifacts (in this")
    print("   case, specify a seeing value).  If the PA is known, restricting it can also")
    print("   help (try -upa 0.5); by default all orientations are searched.")
    print("   If still having issues, e-mail dperley@astro.berkeley.edu for help.")

def cataloghelp():
   print("Catalog info:")
   print("   Leave the catalog field blank will use SDSS if available and USNO otherwise.")
   print("   The catalog query uses wcstools (tdc-www.harvard.edu/wcstools).  However, you")
   print("   can also use your own catalog file on disk if you prefer using -c [filename]")
   print("   The default format is a text file with the first three columns indicating")
   print("   ra, dec, magnitude.  However, you can change the order of the columns by")
   print("   adding, e.g.")
   print("              #:1,2,6")
   print("   to the first line.  In this case, thisw ould indicate that the RA is in the") 
   print("   1st column, dec in the 2nd, and magnitude in the 6th.   The mag column can be")
   print("   omitted completely, although if the catalog is not the same depth as the")
   print("   image this may compromise the search results.")

def examplehelp():
#         12345678901234567890123456789012345678901234567890123456789012345678901234567890
   print("Examples:")
   print(" For an image with provisional WCS header but possibly incorrect PA:")
   print("    autoastrometry image.fits   ")
   print(" For an image with provisional WCS, definitely correct PA:")
   print("    autoastrometry image.fits -upa 0.5")
   print(' For an image with no WCS (or bad WCS) and a pixel scale of 0.35"/pixel:')
   print("    autoastrometry image.fits -px 0.35")
   print(' For an image with no WCS, pixel scale of 0.35", and PA of 121 degrees:')
   print("    autoastrometry image.fits -px 0.35 -pa 121 -upa 0.5")
   print(" For an image with no WCS, pixel scale of 0.35, and positive parity:")    
   print("    autoastrometry image.fits -px 0.35 -inv")
   print(" Use a catalog on disk instead of wcstools query:")    
   print("    autoastrometry image.fits -c catalog.txt")
   print(" Widen the catalog search to 12x12 arcminutes if the pointing is bad:")
   print("    autoastrometry image.fits -b 720")
   print(" Specify seeing (7 pixels) to better exclude cosmic rays and galaxies:")
   print("    autoastrometry image.fits -s 7")
   print(" Combine lots of options for a maximally directed solution:")
   print("    autoastrometry image.fits -px 0.35 -pa 121 -upa 0.5 -inv -b 600 -s 7")
   print(" (Substitute 'autoastrometry' with 'python autoastrometry.py' if not aliased.)")

def outputhelp():
   print("Explanation of output files:")
   print(" DS9 region files -")
   print("   cat.wcs.reg        - WCS positions of all objects in the catalog.")
   print("   det.im.reg         - XY positions of all non-flagged detected objects.")
   print("   matchlines.im.reg  - spoke lines for matched asterisms from XY positions.")
   print("   matchlines.wcs.reg - spoke lines for matched asterisms from WCS positions.")
   print(" Text output -")
   print("   cat.txt            - objects in the catalog as a text file (RA, dec, mag)")
   print("   det.init.txt       - objects in the image as a text file ('RA', 'dec', mag)")
   print("   det.final.txt      - objects in the image as a text file (RA, dec, mag)")
   print("   match.list         - list of matched stars: (X Y RA dec)")
   print(" Image output -")
   print("   a[image.fits]      - output image with astrometry (if not specified with -o)")
   print("   temp.fits          - image with provisional WCS from your -px and -pa inputs")

def inputhelp():
   print("Sextractor input files:")
   print("  You can modify the following sextractor inputs if you choose.  This should")
   print("  rarely be necessary, if ever.")
   print("     sex.config - Overall input file")
   print("     sex.conv   - Convolution matrix for detection")
   print("  These files are written to disk using internal defaults if not present.")


#Some other conceivable options for the future:
  #noncircular PSF
  #instrument defaults

######################################################################
def usage():
    (xdir,xname) = os.path.split(sys.argv[0])
    print("Usage:  %s filename(s) [-px pixelscale -pa PA -inv -b boxsize -s seeing -upa PAunc]" % xname)
    print("     or %s -help for instructions and more options." % xname)

######################################################################
def main():
    
    files=[]

    if (len(sys.argv)==1):
        usage()
        sys.exit(1)

    i=1
    #defaults
    pixelscale = -1
    pa = -999
    uncpa = -1
    inv = 0
    boxsize = -1
    maxrad = -1
    maxellip = -1
    seeing = -1
    tolerance = defaulttolerance
    catalog = ''
    nosolve = 0
    norot = 0
    saturation = -1
    overwrite = False
    outfile = ''
    quiet = False
    userra = -999
    userdec = -999
    
    # nickel : 0.371"/pix, PA=177.7
    # lick Geminicam:  A=0.697"/pix  B=0.672"/pix
     # Nov observations are all offset by a large amount in dec even after equinox
    # Keck: 0.135"/pix, chipPA=90+instPA

    while (i<len(sys.argv)):
        arg=sys.argv[i]
        isarg=0
        if (arg.find("-h") == 0):
            help()
            sys.exit(1)
        if (arg.find("-examp") == 0):
            examplehelp()
            sys.exit(1)
        if (arg.find("-troub") == 0):
            troublehelp()
            sys.exit(1)
        if (arg.find("-catal") == 0):
            cataloghelp()
            sys.exit(1)
        if (arg.find("-output") == 0):
            outputhelp()
            sys.exit(1)
        if (arg.find("-input") == 0):
            inputhelp()
            sys.exit(1)
        if (arg.find("-algor") == 0):
            inputhelp()
            sys.exit(1)
        if (arg.find("-pi")  == 0 or arg.find("-px") == 0):
            pixelscale=float(sys.argv[i+1])
            i+=1
            isarg=1
        if (arg.find("-ra") == 0):
            userra = rasex2deg(sys.argv[i+1])
            i+=1
            isarg=1
        if (arg.find("-dec") == 0):
            userdec = decsex2deg(sys.argv[i+1])
            i+=1
            isarg=1
        if (arg.find("-pa") == 0):
            pa=float(sys.argv[i+1])
            i+=1
            isarg=1
        if (arg.find("-upa") == 0):
            uncpa=float(sys.argv[i+1])
            i+=1
            isarg=1
        if (arg.find("-inv") == 0):
            inv = 1
            isarg = 1
        if (arg.find("-norot") == 0): # some ambiguity with -n (nosolve)
            norot = 1
            isarg = 1
        if (arg.find("-b") == 0):
            boxsize = float(sys.argv[i+1])
            i+=1
            isarg = 1
        if (arg.find("-d") == 0) and (arg.find("-dec") != 0):
            maxrad = float(sys.argv[i+1])
            i+=1
            isarg = 1
        if (arg.find("-s") == 0):
            seeing = float(sys.argv[i+1])
            i+=1
            isarg = 1
        if (arg.find("-x") == 0):
            saturation = float(sys.argv[i+1])
            i+=1
            isarg = 1
        if (arg.find("-t") == 0):
            tolerance = float(sys.argv[i+1])
            i+=1
            isarg = 1
        if (arg.find("-e") == 0):
            maxellip = float(sys.argv[i+1])
            i+=1
            isarg = 1
        if (arg.find("-c") == 0):
            catalog = sys.argv[i+1]
            i+=1
            isarg = 1
        if (arg.find("-q") == 0):
            quiet = True
            isarg = 1
        if (arg.find("-o") == 0):
            outfile = ''
            if len(sys.argv) > i+1:
               if (sys.argv[i+1]).find('-') != 0: 
                  outfile = sys.argv[i+1] #output
                  i+=1
            if outfile == '': overwrite = 1
            isarg = 1
        if (arg.find("-n") == 0 and arg.find("-norot") == -1):
            nosolve = 1
            if len(sys.argv) > i+1:
               if (sys.argv[i+1]).find('-') != 0: 
                  catalog = sys.argv[i+1] #output
                  i+=1
            isarg = 1
        if (not isarg):
            if (len(files) > 0):
                files+=",%s" % arg
            else:
                files=arg
        i+=1

    if len(files) == 0:
        print('No files selected!')
        return

    filenames=files.split(",")

    if (seeing == -1):
        minfwhm = defaultminfwhm #1.5
        maxfwhm = defaultmaxfwhm #40
    else:
        minfwhm = 0.7 * seeing
        maxfwhm = 2 * seeing

    writeparfile()
    if not os.path.exists('sex.config'): writeconfigfile(saturation)

    nimage = len(filenames)
    failures = []
    questionable = []
    multiinfo = []    

    for filename in filenames:
        if len(filenames) > 1: print('Processing', filename)
        if nosolve and catalog=='': catalog = filename+'.cat'
        fitinfo = autoastrometry(filename,pixelscale=pixelscale,pa=pa,inv=inv,uncpa=uncpa,minfwhm=minfwhm,maxfwhm=maxfwhm,maxellip=maxellip,boxsize=boxsize, maxrad=maxrad, userra=userra, userdec=userdec, tolerance=tolerance, catalog=catalog, nosolve=nosolve, overwrite=overwrite, outfile=outfile, saturation=saturation, quiet=quiet, norot=norot)
        if nosolve: continue
        if type(fitinfo)==int: 
           fitinfo = (0,0,0,0,0,0)
        
        multiinfo.append(fitinfo)

        if (fitinfo[0] == 0):   #number of matches
            failures.append(filename)
        if (fitinfo[5] > 2):    #stdev of offset
            questionable.append(filename)
        if nimage > 1: print()
            
    if nimage > 1 and nosolve==0:

        if len(failures) == 0 and len(questionable) == 0:
            print('Successfully processed all images!')
        else:
            print('Finished processing all images.')
        
        if len(questionable) > 0:
            print('The following images solved but have questionable astrometry: ')
            print('    ', end=' ')
            for f in questionable: print(f, end=' ')
            print()
        if len(failures) > 0:
            print('The following images failed to solve: ')
            print('    ', end=' ')
            for f in failures: print(f, end=' ')
            print()

        print("%25s " %'Filename', end=' ') 
        print("%6s %8s (%6s)  %7s %7s (%6s)" % ('#match', 'dPA ', 'stdev', 'dRA', 'dDec', 'stdev'))
        for i in range(len(filenames)):
            info = multiinfo[i]
            print("%25s " % filenames[i], end=' ') 
            if info[0] > 0:
               print("%6d %8.3f (%6.3f)  %7.3f %7.3f (%6.3f)" % info)
            else:
               print("failed to solve")

    try:
       os.remove('temp.param')
    except:
       print('Could not remove temp.param for some reason')

######################################################################
# Running as executable
if __name__=='__main__':
    main()

######################################################################


# some possible future improvements:
# verification to assess low-confidence solutions
# full automatic retry mode (parity detection, etc.)
# dealing with unknown pixel scale
# run wcstools for distortion parameters
# merge catalog check with catalog search to save a query
# improve the CR rejection further... maybe think about recognizing elliptical "seeing"?


