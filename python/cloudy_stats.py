import numpy as np
from utils import robustRMS
import healpy as hp
import matplotlib.pylab as plt

# Make a function that calculates how cloudy the sky is based on a difference image


def cloudyness(diff_image, fwhm=5., sigma_cut=3., 
               airmass_map=None, airmass_limit=None, skyRMS_max = None):
    """
    Parameters
    ----------
    diff_image: np.array
        An array that is the difference of two all-sky images

    sigma_cut: float
        blah

    """

    # XXX.  Maybe the solution is to set a minimum angular scale, but then filter on the
    # scale that has the most power (if that's larger)?

    # Or, there should be some relation between the brightness level of the frame and the
    # RMS of the frame--such that if the median of the frame - dark_sky is ~0, the sigma ~= 0.06,
    # but if frame - dark >> dark, then sigma should scale down.  So that's how I could set skyRMS_max.

    # Also look at imposing some continuity

    unmasked = np.where(diff_image != hp.UNSEEN)[0]
    skyRMS = robustRMS(diff_image[unmasked])
    if skyRMS_max is not None:
        skyRMS = np.min([skyRMS, skyRMS_max])

    # Try to find the angular scale of any clouds?

    smooth_map = hp.sphtfunc.smoothing(diff_image, fwhm=np.radians(fwhm),
                                       verbose=False, iter=1)
    smooth_map[unmasked] = smooth_map[unmasked] - np.median(smooth_map[unmasked])

    outliers = np.where(np.abs(smooth_map[unmasked]) > sigma_cut*skyRMS)[0]
    # Can think about going back to the original map and growing the region that got flagged

    cloud_mask = np.zeros(diff_image.size, dtype=int)
    cloud_mask[unmasked[outliers]] = 1

    nside = hp.npix2nside(np.size(diff_image))
    pix_area = hp.nside2pixarea(nside)

    out_area = outliers.size*pix_area*(180./np.pi)**2

    return out_area




