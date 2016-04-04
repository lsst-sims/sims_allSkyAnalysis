import numpy as np
import healpy as hp
from medBD import medDB, single_frame
from lsst.sims.utils import haversine
import lsst.sims.skybrightness as sb


nside = 32
filter_name = 'R'
sm = sb.SkyModel(mags=True) #sb.SkyModel(mags=True, airglow=False, mergedSpec=False)
dec, ra = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
dec = np.pi/2 - dec
moonLimit = 30. # Degrees
am_limit = 3.

def diff_image(mjd, mjd2, filter_name='R'):
	"""
	Let's just load up a single image and difference it with the one taken before
	"""
	sm.setRaDecMjd(ra,dec,mjd)
	frame = single_frame(mjd, filter_name=filter_name)
	dist2moon = haversine(sm.azs, sm.alts, sm.moonAz, sm.moonAlt)
	frame[np.where(dist2moon < np.radians(moonLimit))] = hp.UNSEEN
	template_frame = single_frame(mjd2, filter_name=filter_name)
	good = np.where( (frame != hp.UNSEEN) & (template_frame != hp.UNSEEN)
	                & (sm.airmass >= 1.) & (sm.airmass <= am_limit) &
	                (~np.isnan(frame)) & (~np.isnan(template_frame)))
	diff = np.zeros(frame.size, dtype=float) + hp.UNSEEN
	diff[good] = frame[good] - template_frame[good]
	return diff



data = np.load('cloud_stats.npz')
moon_alt = data['moon_alts'].copy()
sun_alt = data['sun_alts'].copy()
model_stats = data['model_stats'].copy()
frame_stats = data['frame_stats'].copy()
umjd = data['umjd'].copy()
data.close()


#umjd = medDB(full_select='select DISTINCT(mjd) from medskybrightness;', dtypes=float)

