import numpy as np
import healpy as hp
from medBD import medDB, single_frame
from lsst.sims.utils import haversine
import lsst.sims.skybrightness as sb
from scipy.optimize import curve_fit
import sys

def robustRMS(arr):
    if np.size(arr) == 0:
        return 0
    iqr = np.percentile(arr,75)-np.percentile(arr,25)
    rms = iqr/1.349 #approximation
    return rms

def two_comp_model(x, a, b):
    """
    Function to generate the expected sky

    x : np.array
      two arrays with the maps to be scaled
    """
    result = a*x[0]+b*x[1]
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
    good = np.where((~np.isnan(sm)) & (sky_frame != hp.UNSEEN) & (median_map != hp.UNSEEN)
                    & (~np.isnan(median_map) & (~np.isinf(median_map))))

    # Expect the total flux to be
    # a*median_map_flux + b*skyModel_flux

    median_map_flux = 10.**(-0.8*median_map[good])
    skyModel_flux = 10.**(-0.8*sm[good])
    sky_frame_flux = 10.**(-0.8*sky_frame[good])

    p0 = np.array([1., np.median(sky_frame_flux/skyModel_flux)])
    flux_resid = np.zeros(sky_frame.size, dtype=float)+hp.UNSEEN
    try:
        best_fit, pcov = curve_fit(two_comp_model, [median_map_flux, skyModel_flux],
                                   sky_frame_flux, p0=p0)
        flux_resid[good] = -2.5*np.log10( sky_frame_flux/two_comp_model([median_map_flux,
                                                                        skyModel_flux],*best_fit))
    except:
        best_fit = p0*0

    
    return flux_resid, best_fit



if __name__ == '__main__':

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
    sm = sb.SkyModel(mags=True) #sb.SkyModel(mags=True, airglow=False, mergedSpec=False)

    dec, ra = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    dec = np.pi/2 - dec

    names = ['median_diff', 'rrms', 'frac_outliers']
    types = [float]*3
    frame_stats = np.zeros(umjd.size, dtype=zip(names,types))
    model_stats = np.zeros(umjd.size, dtype=zip(names,types))
    moon_alts = np.zeros(umjd.size, dtype=float)
    sun_alts = np.zeros(umjd.size, dtype=float)
    am_limit = 10.
    outlier_thresh = 3.
    moonLimit = 30. # Degrees

    maxj = np.size(umjd)

    for i in np.arange(0, maxj, 1):

        progress = i/float(maxj)*100
        text = "\rprogress = %.1f%%"%progress
        sys.stdout.write(text)
        sys.stdout.flush()

        mjd = umjd[i]
        frame = single_frame(mjd, filter_name=filter_name)
        sm.setRaDecMjd(ra,dec,mjd)
        moon_alts[i] += sm.moonAlt
        sun_alts[i] += sm.sunAlt
        # mask out anything too close to the moon
        dist2moon = haversine(sm.azs, sm.alts, sm.moonAz, sm.moonAlt)
        frame[np.where(dist2moon < np.radians(moonLimit))] = hp.UNSEEN

        # use the previous exposure as a template, unless there's a big gap
        if mjd - umjd[i-1] < 0.0009:
            mjd_template = umjd[i-1]
            template_frame = single_frame(mjd_template, filter_name=filter_name)
            good = np.where( (frame != hp.UNSEEN) & (template_frame != hp.UNSEEN)
                             & (sm.airmass >= 1.) & (sm.airmass <= am_limit) &
                             (~np.isnan(frame)) & (~np.isnan(template_frame)))

            diff = np.zeros(frame.size, dtype=float) + hp.UNSEEN
            diff[good] = frame[good] - template_frame[good]
            frame_stats['median_diff'][i] = np.median(diff[good])
            frame_stats['rrms'][i] = robustRMS(diff[good])
            outliers = np.where(np.abs(diff[good]) > outlier_thresh *frame_stats['rrms'][i])
            if np.size(diff[good]) != 0:
                frame_stats['frac_outliers'][i] = np.size(outliers[0])/float(np.size(diff[good]))


        model_mags = sm.returnMags()
        resid,fit_params = residMap(frame, mm['median'+filter_name], model_mags['r'])
        if np.max(np.abs(fit_params) != 0):
            unmasked = np.where(resid != hp.UNSEEN)
            if np.size(unmasked[0]) > 1:
                model_stats['median_diff'][i] = np.median(resid[unmasked])
                model_stats['rrms'][i] = robustRMS(resid[unmasked])
                outliers = np.where(np.abs(resid[unmasked]) > outlier_thresh *model_stats['rrms'][i])
                model_stats['frac_outliers'][i] = np.size(outliers[0])/float(np.size(unmasked[0]))


    # Save the output for later
    np.savez('cloud_stats.npz', model_stats=model_stats, frame_stats=frame_stats, 
             moon_alts=moon_alts, sun_alts=sun_alts, umjd=umjd)
