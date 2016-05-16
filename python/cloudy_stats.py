import numpy as np
from utils import robustRMS
import healpy as hp
import matplotlib.pylab as plt

# Make a function that calculates how cloudy the sky is based on a difference image


def cloudyness(diff_image, fwhm=5., sigma_cut=3., sigma_second = 3.,
               airmass_map=None, airmass_limit=None, skyRMS_max = None,
               grow_iter=3, grow_fwhm=5., grow_lower_limit=0.1):
    """
    Parameters
    ----------
    diff_image: np.array
        An array that is the difference of two all-sky images

    sigma_cut: float
        blah

    Returns
    -------

    out_area:  float
        area (sq degrees) that are flagged as cloudy
    cloud_mask:  np.array
        Array same size as diff_image with values of -1, 0, 1 to show which 
        pixels have been flagged as negative or positive.
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

    # outliers = np.where(np.abs(smooth_map[unmasked]) > sigma_cut*skyRMS)[0]
    # Can think about going back to the original map and growing the region that got flagged

    cloud_mask = np.zeros(diff_image.size, dtype=int)
    highOutliers = np.where((smooth_map > sigma_cut*skyRMS) & (diff_image > sigma_second * skyRMS) &
                            (diff_image != hp.UNSEEN))
    lowOutliers = np.where((smooth_map < -1*sigma_cut*skyRMS) & (diff_image < -1*sigma_second * skyRMS) &
                           (diff_image != hp.UNSEEN))

    cloud_mask[highOutliers] = 1
    cloud_mask[lowOutliers] = -1

    # rather than loop, let's just use smoothing in clever ways!
    if np.max(np.abs(cloud_mask)) != 0:
        for ack in np.arange(grow_iter):
            expanded_mask = hp.sphtfunc.smoothing(cloud_mask, fwhm=np.radians(grow_fwhm),
                                                  verbose=False, iter=1)
            cloud_mask[np.where((expanded_mask > grow_lower_limit) & (diff_image > 0))] = 1
            cloud_mask[np.where((expanded_mask < -1.*grow_lower_limit) & (diff_image < 0))] = -1

    # What is the most efficient way to loop over this stuff? 
   # if np.max(np.abs(cloud_mask)) != 0:
   #     changed_pix =1
   #     while changed_pix != 0:
   #         clear_pix = np.where(cloud_mask == 0)[0]
   #         for i in clear_pix:

    nside = hp.npix2nside(np.size(diff_image))
    pix_area = hp.nside2pixarea(nside)
    outliers = np.where(cloud_mask != 0)[0]
    out_area = outliers.size*pix_area*(180./np.pi)**2

    return out_area, cloud_mask




