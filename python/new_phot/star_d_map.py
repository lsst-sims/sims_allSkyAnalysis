import numpy as np
import healpy as hp
from lsst.sims.utils import healbin
import matplotlib.pylab as plt

# Just look at the number of stars we're getting on the sky

nside = 8

data = np.load('phot.npz')
phot = data['phot'].copy()
data.close()

umjd = np.unique(phot['mjd'])
phot.sort(order='mjd')
left = np.searchsorted(phot['mjd'], umjd)
right = np.searchsorted(phot['mjd'], umjd, side='right')

nstars = right-left

goodMJD = umjd[np.where(nstars == nstars.max())]
good = np.where(phot['mjd'] == goodMJD)
hmap = healbin(phot['RA'][good], phot['dec'][good], phot['dec'][good], nside=nside, reduceFunc=np.size)

fig = plt.figure(1)
hp.mollview(hmap, fig=1, title='%i stars' % np.size(good[0]), unit='N stars')
fig.savefig('stellar_density.png')
