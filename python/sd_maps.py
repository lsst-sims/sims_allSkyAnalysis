import numpy as np
import healpy as hp

# read in the stellar density maps and be able to return hp maps as a function of MJD

class starD(object):
	"""
	Read in the saved stellar density maps from the all-sky camera
	"""

	def __init__(self, nside=8):
		filename = 'starDensity_%i.dat' % nside
		names = ['mdj', 'hpid', 'nstars']
		types=[float, int, int]
		self.nside = nside
		self.empty = np.empty(hp.nside2npix(nside), dtype=int)
		self.empty.fill(hp.UNSEEN)

		self.data = np.loadtxt(filename, dtype=zip(names, types), delimiter=',', skiprows=1)
		self.umjd = np.unique(self.data['mjd'])


	def __call__(self, mjd):
		if mjd not in self.umjd:
			return 0

		left = np.searchsorted(self.data['mjd'], mjd)
		right = np.searchsorted(self.data['mjd'], mjd, side='right')
		frame = self.data[left:right]

		hpmap = self.empty.copy()
		hpmap[frame['hpid']] = frame['nstars']
		return hpmap
