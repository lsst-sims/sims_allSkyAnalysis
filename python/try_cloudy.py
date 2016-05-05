import numpy as np
import healpy as hp
import matplotlib.pylab as plt
from medDB import single_frame, medDB
from lsst.sims.skybrightness import stupidFast_RaDec2AltAz
from lsst.sims.utils import calcLmstLast, Site
from utils import robustRMS
import sys
from cloudy_stats import cloudyness
# Let's just load up some images and play with it 

# XXXXX
# OK, looks like there are issues with the data having been integers.  
# I can probably find the part of the bias level and subtract it off...

# XXXX
# Still have an interesting issue with 0.5 since np.median will mean values. Fantastic.
# I still kind of wish the 

def fixBias(frame):
    good = np.where((frame != hp.UNSEEN) & (frame > 0))
    biasLevel = frame[good].min()
    result = frame - biasLevel
    result[np.where(frame == hp.UNSEEN)] = hp.UNSEEN
    return result



umjd = medDB(full_select='select DISTINCT(mjd) from medskybrightness;', dtypes=float)

# nstart = 5506+67  # partly cloudy frame
nstart = 6447  # very cloudy frame
# nstart = 982  # clear, with moon
previous = single_frame(umjd[nstart-1])
nside = hp.npix2nside(previous.size)
mjd = umjd[nstart]


site = Site('LSST')
dec, ra = hp.pix2ang(nside, np.arange(previous.size))
dec = np.pi/2. - dec


lmst, last = calcLmstLast(mjd, site.longitude_rad)
lmst = lmst/12.*180.
alt, az = stupidFast_RaDec2AltAz(ra, dec, site.latitude_rad, site.longitude_rad, mjd)

frame = single_frame(mjd, filter_name='R')
seen = np.where((frame != hp.UNSEEN) & (previous != hp.UNSEEN))
unseen = np.where((frame == hp.UNSEEN) | (previous == hp.UNSEEN))

frame = fixBias(frame)
previous = fixBias(previous)


diff = frame - previous
diff[unseen] = hp.UNSEEN

diff_frac = diff/previous
diff_frac[unseen] = hp.UNSEEN
diff_frac[np.where(previous == 0)] = hp.UNSEEN


gdiff = np.where((diff != hp.UNSEEN) & (~np.isnan(diff)))[0]
if np.size(gdiff) > 0:
    cf = cloudyness(diff_frac) / (gdiff.size*hp.nside2pixarea(nside)*(180./np.pi)**2)
    print 'cloudy fraction, 5 deg scale = %.2f' % cf 
    cf = cloudyness(diff_frac, fwhm=30.) / (gdiff.size*hp.nside2pixarea(nside)*(180./np.pi)**2)
    print 'cloudy fraction, 30 deg scale = %.2f' % cf 