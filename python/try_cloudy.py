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



umjd =  medDB(full_select='select DISTINCT(mjd) from medskybrightness;', dtypes=float)
nstart = 150
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



#frame[seen] = frame[seen] - np.median(frame[seen])
#previous[seen] = previous[seen]  - np.median(previous[seen])


simple_diff = frame - previous
simple_diff[unseen] = hp.UNSEEN

diff = frame / previous - 1.
out = np.where((np.isnan(frame)) | (np.isnan(previous)) | (frame == hp.UNSEEN) | 
		      (previous == hp.UNSEEN) | (alt < np.radians(10.)))

diff[out] = hp.UNSEEN
frame[out] = hp.UNSEEN
# maybe rotate based on LMST and latitude?
gdiff = np.where( (diff != hp.UNSEEN) & (np.isnan(diff) == False))[0]
if np.size(gdiff) > 0:
	rms = robustRMS(diff[gdiff])
	median_value = np.median(diff[gdiff])
	#nout = np.size(np.where( (np.abs(diff[gdiff] - np.median(diff[gdiff]) ) > outlier_mag) & (alt[gdiff] > alt_limit))[0])
	#nabove = float(np.size(np.where(alt[gdiff] > alt_limit)[0]))*100
	cf = cloudyness(diff[gdiff], sigma_max = 0.06)

