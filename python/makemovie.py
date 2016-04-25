#!/usr/bin/env python
import numpy as np
import healpy as hp
import matplotlib.pylab as plt
from medDB import single_frame, medDB
from lsst.sims.skybrightness import stupidFast_RaDec2AltAz
from lsst.sims.utils import calcLmstLast, Site
from utils import robustRMS

# Let's try making some simple movies to see what the frames and different images look like


if __name__ == '__main__':

	cannonFilter = 'R'
	

	outdir = 'MoviePlots'

	#data = np.load('cloud_stats.npz')
	#umjd = data['umjd'].copy()
	#sun_alt = data['sun_alts'].copy()
	#moon_alt = data['moon_alts']
	#frame_stats = data['frame_stats'].copy()
	#data.close()

	umjd =  medDB(full_select='select DISTINCT(mjd) from medskybrightness;', dtypes=float)

	RdBu = plt.get_cmap('RdBu')
	RdBu.set_bad('gray')
	RdBu.set_under('w')

	nframes = umjd.size -1 
	
	print 'making %i frames' % nframes
	nstart = 1
	

	previous = single_frame(umjd[nstart-1])
	nside = hp.npix2nside(previous.size)
	site = Site('LSST')
	dec, ra = hp.pix2ang(nside, np.arange(previous.size))
	dec = np.pi/2. - dec

	median_map = np.load('sky_maps.npz')
	median_filt = median_map['sky_maps']['median%s' % cannonFilter].copy()
	median_map.close()
	median_filt[np.isnan(median_filt)] = hp.UNSEEN


	# Load up the stellar density
	# Arrg, the healpixed database doesn't overlap with the stellar database.
	# starmap = starD()

	outlier_mag = 0.1
	alt_limit = np.radians(23.) 
	fracs_out = []
	med_diff_frame = []
	rms_diff_frame = []
	med_diff_med = []

	for i,mjd in enumerate(umjd[nstart:nstart+nframes]):

		lmst, last = calcLmstLast(mjd, site.longitude_rad)
		lmst = lmst/12.*180.
		alt, az = stupidFast_RaDec2AltAz(ra, dec, site.latitude_rad, site.longitude_rad, mjd)
		fig = plt.figure()
		frame = single_frame(mjd, filter_name=cannonFilter)
		# Need to crop off high airmass pixels. use stupid_fast
		diff = frame / previous - 1. #(10.**(-0.4*frame) / 10.**(-.4*previous)) - 1.
		previous = frame.copy()
		out = np.where((np.isnan(frame)) | (np.isnan(previous)) | (frame == hp.UNSEEN) | 
		               (previous == hp.UNSEEN) | (alt < np.radians(10.)))
		diff[out] = hp.UNSEEN
		frame[out] = hp.UNSEEN
		# maybe rotate based on LMST and latitude?
		gdiff = np.where( (diff != hp.UNSEEN) & (np.isnan(diff) == False))[0]
		if np.size(gdiff) > 0:
			rms = robustRMS(diff[gdiff])
			median_value = np.median(diff[gdiff])
			nout = np.size(np.where( (np.abs(diff[gdiff] - np.median(diff[gdiff]) ) > outlier_mag) & (alt[gdiff] > alt_limit))[0])
			nabove = float(np.size(np.where(alt[gdiff] > alt_limit)[0]))*100
			if nabove != 0:
				nout = nout/nabove
			else:
				nout = -666
		else:
			rms = hp.UNSEEN
			median_value = hp.UNSEEN
			nout = -666
		rms_diff_frame.append(rms)
		med_diff_frame.append(median_value)

		fracs_out.append(nout)
		hp.mollview(frame, sub=(2,2,1), rot=(lmst, site.latitude,0), unit='counts', 
		            title='%.2f' % mjd, min=5., max=2000., norm='log')
		hp.mollview(diff, sub=(2,2,2), min=-.3, max=.3,  rot=(lmst, site.latitude,0), 
		            cmap=RdBu, unit='(frame-prev)/prev', 
		            title=r'$\sigma$=%.2f, percent out=%i' % (rms, nout))
		diff2 = (10.**(-.4*frame)/10.**(-.4*median_filt)) -1.
		out = np.where((frame == hp.UNSEEN) | (median_filt == hp.UNSEEN) | (alt < np.radians(10.)))
		diff2[out] = hp.UNSEEN
		good = np.where(diff2 != hp.UNSEEN)
		hp.mollview(diff2, sub=(2,2,3), min=-2, max=2, rot=(lmst, site.latitude,0), 
		            cmap=RdBu, unit='(frame-median)/median (flux)', title='median = %.1f' % np.median(diff2[good]))
		med_diff_med.append(np.median(diff2[good]))

		#diff_correlation = (diff*diff2)**2
		#out = np.where((diff == hp.UNSEEN) | (diff2 == hp.UNSEEN))
		#diff[out] = hp.UNSEEN
		#diff2[out] = hp.UNSEEN
		#good = np.where(diff != hp.UNSEEN)
		#diff_correlation = ((diff-np.median(diff[good]))*(diff2-np.median(diff2[good])))**2
		#diff_correlation[out] = hp.UNSEEN
		#hp.mollview(diff_correlation, sub=(2,2,4), min=0, max=1, rot=(lmst, site.latitude,0), 
		#            unit='(power)', title='correlation')
		#sm, sm_diff = starmap(mjd)
		#hp.mollview(sm, sub=(2,2,3), rot=(lmst, site.latitude,0), unit='N stars', title='')
		#hp.mollview(sm_diff, sub=(2,2,4), rot=(lmst, site.latitude,0), unit='N stars', title='')

		fig.savefig('%s/%05i_.png' % (outdir, i))
		plt.close(fig)
		
	# then just call ffmpeg like here: https://trac.ffmpeg.org/wiki/Create%20a%20video%20slideshow%20from%20images
	# open with vLC
	#  ffmpeg -framerate 10 -pattern_type glob -i '*.png'  out.mp4

	# Save the stats that came out
	np.savez('movie_stats.npz', med_diff_med=med_diff_med, med_diff_frame=med_diff_frame, 
	         moon_alt=moon_alt, sun_alt=sun_alt, umjd=umjd, fracs_out=fracs_out)
