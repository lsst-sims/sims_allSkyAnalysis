import numpy as np
from medDB import medDB
import healpy as hp
import matplotlib.pylab as plt
from lsst.sims.utils import Site
import ephem
from utils import robustRMS



# Generate the median of the median maps.

nside = 32
airmass_limit = 3.

names = ['medianR','medianG','medianB', 'stdR', 'stdG', 'stdB', 'numberR']
types = [float]*len(names)
sky_maps = np.zeros(hp.nside2npix(nside), dtype=zip(names, types))
sky_maps.fill(hp.UNSEEN)

umjd = medDB(full_select='select DISTINCT(mjd) from medskybrightness;', dtypes=float)
site = Site('LSST')
sun = ephem.Sun()
moon = ephem.Moon()

obs = ephem.Observer()
obs.lat, obs.lon, obs.elevation = site.latitude_rad, site.longitude_rad,site.height
doff = ephem.Date(0)-ephem.Date('1858/11/17')
udjd = umjd - doff
moonAlts = umjd*0.
sunAlts = umjd*0.

print 'Computing moon and sun altitudes'
for i,mjd in enumerate(umjd):
    obs.date = udjd[i]
    moon.compute(obs)
    sun.compute(obs)
    moonAlts[i] = moon.alt
    sunAlts[i] = sun.alt

goodDates = umjd[np.where((moonAlts < 0) & (sunAlts < np.radians(-18.)))]

print 'looping over each healpixel'
for i in np.arange(hp.nside2npix(nside)):
    data = medDB(where_clause='hpindex = %i' % i)
    data = data[np.in1d(data['mjd'], goodDates)]
    data = data[np.where(data['airmass'] <= airmass_limit)]
    if np.size(data) > 0:
        for key in ['R','G','B']:
            sky_maps['median'+key][i] = np.median(data[key][~np.isnan(data[key])])
            if np.size(data[key][~np.isnan(data[key])]) > 3:
                sky_maps['std'+key][i] = robustRMS(data[key][~np.isnan(data[key])])
        sky_maps['numberR'][i] = np.size(data['R'][~np.isnan(data[key])])

print 'Finished generating map'
np.savez('sky_maps.npz', sky_maps=sky_maps)
hp.mollview(sky_maps['medianR'])
plt.savefig('median_r.png')

#plt.show()
