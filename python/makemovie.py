#!/usr/bin/env python
import numpy as np
import healpy as hp
import matplotlib.pylab as plt
from medBD import single_frame
from lsst.sims.skybrightness import stupidFast_RaDec2AltAz
from lsst.sims.utils import calcLmstLast, Site
from utils import robustRMS

# Let's try making some simple movies to see what the frames and different images look like


if __name__ == '__main__':

	nframes = 400
	nstart = 500

	outdir = 'MoviePlots'

	data = np.load('cloud_stats.npz')
	umjd = data['umjd'].copy()
	frame_stats = data['frame_stats'].copy()
	data.close()

	RdBu = plt.get_cmap('RdBu')
	RdBu.set_bad('gray')
	RdBu.set_under('w')


	previous = single_frame(umjd[nstart-1])
	nside = hp.npix2nside(previous.size)
	site = Site('LSST')
	dec, ra = hp.pix2ang(nside, np.arange(previous.size))
	dec = np.pi/2. - dec

	median_map = np.load('sky_maps.npz')
	median_r = median_map['sky_maps']['medianR'].copy()
	median_map.close()


	# Load up the stellar density
	# Arrg, the healpixed database doesn't overlap with the stellar database.
	# starmap = starD()

	outlier_mag = 0.1
	alt_limit = np.radians(23.) 

	for i,(mjd, stats) in enumerate(zip(umjd[nstart:nstart+nframes], 
	                                    frame_stats[nstart:nstart+nframes])):

		lmst, last = calcLmstLast(mjd, site.longitude_rad)
		lmst = lmst/12.*180.
		alt, az = stupidFast_RaDec2AltAz(ra, dec, site.latitude_rad, site.longitude_rad, mjd)
		fig = plt.figure()
		frame = single_frame(mjd)
		# Need to crop off high airmass pixels. use stupid_fast
		diff = frame - previous
		previous = frame.copy()
		out = np.where((np.isnan(frame)) | (np.isnan(previous)) | (frame == hp.UNSEEN) | 
		               (previous == hp.UNSEEN) | (alt < np.radians(10.)))
		diff[out] = hp.UNSEEN
		frame[out] = hp.UNSEEN
		# maybe rotate based on LMST and latitude?
		gdiff = np.where( (diff != hp.UNSEEN) & (np.isnan(diff) == False))[0]
		rms = robustRMS(diff[gdiff])
		nout = np.size(np.where( (np.abs(diff[gdiff]) > outlier_mag) & (alt[gdiff] > alt_limit))[0])
		nout = nout/float(np.size(gdiff))*100
		hp.mollview(frame, sub=(2,2,1), rot=(lmst, site.latitude,0), min=-7, max=-2.5, unit='mag', 
		            title='%.2f' % mjd)
		hp.mollview(diff, sub=(2,2,2), min=-0.25, max=0.25,  rot=(lmst, site.latitude,0), 
		            cmap=RdBu, unit='(mag)', 
		            title=r'$\sigma$=%.2f, Nout=%i' % (rms, nout))
		diff = frame - median_r
		out = np.where((frame == hp.UNSEEN) | (median_r == hp.UNSEEN) | (alt < np.radians(10.)))
		diff[out] = hp.UNSEEN
		hp.mollview(diff, sub=(2,2,3), min=-2, max=2, rot=(lmst, site.latitude,0), 
		            cmap=RdBu, unit='(mag)', title='')

		#sm, sm_diff = starmap(mjd)
		#hp.mollview(sm, sub=(2,2,3), rot=(lmst, site.latitude,0), unit='N stars', title='')
		#hp.mollview(sm_diff, sub=(2,2,4), rot=(lmst, site.latitude,0), unit='N stars', title='')

		fig.savefig('%s/%04i_.png' % (outdir, i))
		plt.close(fig)
		

	# then just call ffmpeg like here: https://trac.ffmpeg.org/wiki/Create%20a%20video%20slideshow%20from%20images
	# open with vLC
	#  ffmpeg -framerate 10 -pattern_type glob -i '*.png'  out.mp4
