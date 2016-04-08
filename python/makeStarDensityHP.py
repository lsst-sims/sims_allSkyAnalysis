import numpy as np
import lsst.sims.skybrightness as sb
import healpy as hp
from lsst.sims.utils import healbin


nside = 4
filt = 'R'

mjds = []
hpids = []
nstars = []

# Crop off junky early things
minMJD = 56900
minID,mjd = sb.allSkyDB(0,'select min(ID) from dates where mjd > %i;' % minMJD, dtypes='int')
maxID,mjd = sb.allSkyDB(0,'select max(ID) from dates where mjd > %i;' % minMJD, dtypes='int')


# Temp pad for testing
# maxID = minID+100

full_range = float(maxID - minID)
for i,dateID in enumerate(np.arange(minID.max(),maxID.max())):
	print i, mjd
	data,mjd = sb.allSkyDB(dateID)#  , sqlQ=sqlQ, dtypes=dtypes)
	starmap = healbin(data['ra'], data['dec'], data['dec'], 
	                  nside=nside, reduceFunc=np.size)
	good = np.where(starmap != hp.UNSEEN)[0]
	mjds.extend([mjd]*good.size)
	hpids.extend(good.tolist())
	nstars.extend(starmap[good].tolist())
	

	
print 'ack'
names = ['mjd', 'hpid', 'nstars']
types = [float, int, int]
result = np.array(zip(mjds, hpids, nstars), dtype=zip(names, types))

np.savez('starDensity_nside%i.npz' % nside, starmap=result, nside=nside)
