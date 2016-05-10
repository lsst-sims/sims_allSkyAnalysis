import numpy as np
import healpy as hp
import matplotlib.pylab as plt
from medDB import single_frame, medDB
from lsst.sims.skybrightness import stupidFast_RaDec2AltAz
from lsst.sims.utils import calcLmstLast, Site
from utils import robustRMS
import sys
from cloudy_stats import cloudyness

# Try and find a good curve for the RMS as a function of median sky brightness


def fixBias(frame):
    good = np.where((frame != hp.UNSEEN) & (frame > 0))
    biasLevel = frame[good].min()
    result = frame - biasLevel
    result[np.where(frame == hp.UNSEEN)] = hp.UNSEEN
    return result


skyMaps = np.load('sky_maps.npz')
umjd = skyMaps['umjd'].copy()

# only do every few frames
nskip = 100

skyMedians = []
skyRMSs = []

previous = single_frame(umjd[0])
nside = hp.npix2nside(previous.size)

site = Site('LSST')
dec, ra = hp.pix2ang(nside, np.arange(previous.size))
dec = np.pi/2. - dec

imax = 10000  # umjd.size

for i in np.arange(0, imax, nskip):
    previous = single_frame(umjd[i-1])
    frame = single_frame(umjd[i])
    mjd = umjd[i]

    # lmst, last = calcLmstLast(mjd, site.longitude_rad)
    # lmst = lmst/12.*180.
    # alt, az = stupidFast_RaDec2AltAz(ra, dec, site.latitude_rad, site.longitude_rad, mjd)
    frame = single_frame(mjd, filter_name='R')
    seen = np.where((frame != hp.UNSEEN) & (previous != hp.UNSEEN))
    unseen = np.where((frame == hp.UNSEEN) | (previous == hp.UNSEEN))
    if unseen[0].size != frame.size:

        frame = fixBias(frame)
        previous = fixBias(previous)

        diff = frame - previous
        diff[unseen] = hp.UNSEEN

        diff_frac = diff/previous
        diff_frac[unseen] = hp.UNSEEN
        diff_frac[np.where(previous == 0)] = hp.UNSEEN

        gdiff = np.where((diff_frac != hp.UNSEEN) & (~np.isnan(diff_frac)))[0]

        skyMedians.append(np.median(frame[gdiff]))
        skyRMSs.append(np.std(diff_frac[gdiff]))

