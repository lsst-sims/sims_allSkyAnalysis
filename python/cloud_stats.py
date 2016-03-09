import numpy as np
import healpy as hp
from medBD import medDB, single_frame
from lsst.sims.utils import Site
import ephem
import lsst.sims.skybrightness as sb

def robustRMS(arr):
    iqr = np.percentile(arr,75)-np.percentile(arr,25)
    rms = iqr/1.349 #approximation
    return rms

def two_comp_model(x,map1=None,map2=None):
    """
    Function to generate the expected sky
    """
    result = x[0]*map1+x[1]*map2
    return result


def residMap(sky_frame, median_map, sm, moon_mask=10.):
    """
    Compute the residuals in a sky frame.

    sky_frame:  Healpix map of the sky
    median_map: The expected dark-time brightness map
    sm: Skymodel object that can return a mag map
    moon_mask: How far around the moon to mask (degrees)
    """

    # build a decent mask.

    # Expect the total flux to be
    # a*median_map_flux + b*skyModel_flux




# Loop through each frame, maybe fit a scaled background frame and a background model? Make some threshold and
# say what fraction of the sky is outside the threshold?

# Not sure what to do about airmass.  it looks like there's a good correlation within a night, but not clear if it varies night-to-night

# Load up the median sky map
mm = np.load('sky_maps.npz')
mm = mm['sky_maps'].copy()

# everything is in nside = 32
nside = 32

umjd = medDB(full_select='select DISTINCT(mjd) from medskybrightness;', dtypes=float)

filter_name = 'R'

# Load up the sky templates, but not with the components that should be in the median background.
sm = sb.SkyModel(mags=True, airglow=False, mergedSpec=False)

dec, ra = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
dec = np.pi/2 - dec

names = ['median_diff', 'rrms', 'frac_outliers']
types = [float]*3
frame_stats = np.zeros(umjd.size, dtype=zip(names,types))

am_limit = 10.
outlier_thresh = 3.
for i in np.arange(0,np.size(umjd), 20):

    mjd = umjd[i]
    # use the previous exposure as a template, unless there's a big gap
    if mjd - umjd[i-1] < 0.0009:
        mjd_template = umjd[i-1]
        frame = single_frame(mjd, filter_name=filter_name)
        template_frame = single_frame(mjd_template, filter_name=filter_name)

        sm.setRaDecMjd(ra,dec,mjd)

        good = np.where( (frame != hp.UNSEEN) & (template_frame != hp.UNSEEN)
                         & (sm.airmass >= 1.) & (sm.airmass <= am_limit) &
                         (~np.isnan(frame)) & (~np.isnan(template_frame)))

        #mags = sm.returnMags()

        diff = frame[good] - template_frame[good]
        frame_stats['median_diff'][i] = np.median(diff)
        frame_stats['rrms'][i] = robustRMS(diff)
        outliers = np.where(np.abs(diff) > outlier_thresh *frame_stats['rrms'][i])
        frame_stats['frac_outliers'][i] = np.size(outliers[0])/float(np.size(diff))
