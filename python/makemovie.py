import numpy as np
import healpy as hp
import matplotlib.pylab as plt
from medBD import medDB, single_frame
from lsst.sims.skybrightness import stupidFast_RaDec2AltAz
from lsst.sims.utils import calcLmstLast, Site

# Let's try making some simple movies to see what the frames and different images look like


# Set up side-by-side axes

# umjd = medDB(full_select='select DISTINCT(mjd) from medskybrightness;', dtypes=float)

nframes = 30
nstart = 500

outdir = 'MoviePlots'

data = np.load('cloud_stats.npz')
umjd = data['umjd'].copy()
data.close()


previous = single_frame(umjd[nstart-1])
nside = hp.npix2nside(previous.size)
site = Site('LSST')
dec, ra = hp.pix2ang(nside, np.arange(previous.size))
dec = np.pi/2. - dec

for i,mjd in enumerate(umjd[nstart:nstart+nframes]):

	lmst, last = calcLmstLast(mjd, site.longitude_rad)
	lmst = lmst/12.*180.
	alt, az = stupidFast_RaDec2AltAz(ra, dec, site.latitude_rad, site.longitude_rad, mjd)
	fig = plt.figure()
	frame = single_frame(mjd)
	# Need to crop off high airmass pixels. use stupid_fast
	diff = frame - previous
	previous = frame.copy()
	out = np.where((frame == hp.UNSEEN) | (previous == hp.UNSEEN) | (alt < np.radians(10.)))
	diff[out] = hp.UNSEEN
	frame[out] = hp.UNSEEN
	# maybe rotate based on LMST and latitude?
	hp.mollview(frame, sub=(1,2,1), rot=(lmst, site.latitude,0))
	hp.mollview(diff, sub=(1,2,2), min=-0.5, max=0.5,  rot=(lmst, site.latitude,0))
	fig.savefig('%s/%i_.png' % (outdir, i))
	plt.close(fig)
	

# then just call ffmpeg like here: http://user.astro.columbia.edu/~robyn/ffmpeghowto.html