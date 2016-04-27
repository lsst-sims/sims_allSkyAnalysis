import numpy as np 
import matplotlib.pylab as plt
from photutils import CircularAperture
from photutils import SkyCircularAperture, CircularAnnulus
from photutils import aperture_photometry
from astropy.io import fits
import os
from astropy.table import hstack


# I want to look at the simple statistics of apperture photometry of a single star over a few frames


fits_path = '/Users/yoachim/Scratch/allSkyCamera/fits/lsst-web.ncsa.illinois.edu/~coughlin/allsky/data/FITS/ut011316/R/'

n_list = [528, 530, 532, 534, 536, 538]
positions_list = [[(1900,1383), (1981,1340)], [(1902, 1382), (1984,1339)], [(1908,1381),(1988,1338)], 
				  [(1912, 1380), (1992,1337)], [(1916,1380), (1996, 1335)], [(1920,1377), (2000,1334)] ]

fluxes = []
for n, positions in zip(n_list, positions_list):

	filename = 'ut011316.0'+str(n)+'.long.R.fits'


	hdulist = fits.open(os.path.join(fits_path, filename))
	data = hdulist[0].data

	# Let's do a cude background
	bg = np.median(data[1500:2000, 1000:1500])
	data = data - bg

	# Create appertures in pixel coords
	# bright star and fainter star
	apertures = CircularAperture(positions, r=3.)
	annulus_apertures = CircularAnnulus(positions, r_in=6., r_out=8.)

	raw_flux = aperture_photometry(data, apertures)
	bkgflux_table = aperture_photometry(data, annulus_apertures)
	phot_table = hstack([raw_flux, bkgflux_table], table_names=['raw', 'bkg'])

	bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area()

	bkg_sum = bkg_mean * apertures.area()
	final_sum = phot_table['aperture_sum_raw'] - bkg_sum
	phot_table['residual_aperture_sum'] = final_sum
	print(phot_table['residual_aperture_sum'])
	ack = np.array(positions).T
	#plt.plot(ack[0], ack[1], 'o')
	fluxes.append(phot_table['residual_aperture_sum'].data.tolist())

fluxes = np.array(fluxes)
for i in np.arange(2):
	print 'mean, std = %f, %f' % (fluxes[:,i].mean(), fluxes[:,i].std())
