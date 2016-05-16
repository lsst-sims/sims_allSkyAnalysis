from photutils import datasets
from astropy.stats import sigma_clipped_stats
import matplotlib.pylab as plt
from astropy.io import fits
from photutils import daofind, CircularAperture
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize


# trying things out from http://photutils.readthedocs.io/en/latest/photutils/detection.html

fitsfile = '/Users/yoachim/Scratch/allSkyCamera/fits/lsst-web.ncsa.illinois.edu/~coughlin/allsky/data/FITS/ut011316/R/ut011316.0650.long.R.fits'

hdulist = fits.open(fitsfile)
hdu = hdulist[0]

data = hdu.data[800:900, 800:900]
mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5)
print(mean, median, std)


data = hdu.data
sources = daofind(data - median, fwhm=3.0, threshold=5.*std)
print 'Found %i sources' % len(sources)

positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r=4.)
norm = ImageNormalize(stretch=SqrtStretch(), vmin=2000, vmax=3000)
plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
plt.title('%i Sources from a single all-sky frame' % len(sources))
apertures.plot(color='blue', lw=1.5, alpha=0.5)
plt.savefig(filename='Plots/example_phot.png')
