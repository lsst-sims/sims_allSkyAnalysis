import numpy as np
import healpy as hp

# read in the stellar density maps and be able to return hp maps as a function of MJD

class starD(object):
	"""
	Read in the saved stellar density maps from the all-sky camera
	"""

	def __init__(self, nside=8):
		filename = 'starDensity_nside%i.npz' % nside
		ack = np.load(filename)
		self.data = ack['starmap'].copy()
		ack.close()
		self.nside = nside
		self.umjd = np.unique(self.data['mjd'])
		self.empty = np.empty(hp.nside2npix(nside))
		self.empty.fill(hp.UNSEEN)
		self._gen_median()

	def _gen_median(self):
		"""
		Loop over all the healpixels and make a median image
		"""
		self.median_im = self.empty.copy()
		self.data.sort(order='hpid')
		uhpid = np.unique(self.data['hpid'])
		lefts = np.searchsorted(self.data['hpid'], uhpid)
		rights = np.searchsorted(self.data['hpid'], uhpid, side='right')
		for hpid, left, right in zip(uhpid, lefts, rights):
			self.median_im[hpid] = np.percentile(self.data['nstars'][left:right], 90) #np.median(self.data['nstars'][left:right])

		self.data.sort(order='mjd')

	def single_frame(self, mjd):
		if mjd not in self.umjd:
			return self.empty

		left = np.searchsorted(self.data['mjd'], mjd)
		right = np.searchsorted(self.data['mjd'], mjd, side='right')
		frame = self.data[left:right]

		hpmap = self.empty.copy()
		hpmap[frame['hpid']] = frame['nstars']
		return hpmap

	def __call__(self, mjd):

		frame_now = self.single_frame(mjd)
		diff = self.median_im - frame_now
		mask = np.where((frame_now == hp.UNSEEN) | (self.median_im == hp.UNSEEN))
		diff[mask] = hp.UNSEEN
		return frame_now, diff