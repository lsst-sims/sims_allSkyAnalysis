import numpy as np
import lsst.sims.skybrightness as sb
import healpy as hp
from lsst.sims.utils import healbin


nside = 8
filt = 'R'

density_maps =[]
mjds = []

# Crop off junky early things
minMJD = 56900
minID,mjd = sb.allSkyDB(0,'select min(ID) from dates where mjd > %i;' % minMJD, dtypes='int')
maxID,mjd = sb.allSkyDB(0,'select max(ID) from dates where mjd > %i;' % minMJD, dtypes='int')


# Temp pad for testing
# maxID = minID+100

outfile = open('starDensity_%i.dat' % nside,'w')

print >>outfile, '## mjd, healpixID (nside=%i), Nstars' % nside


full_range = float(maxID - minID)
for i,dateID in enumerate(np.arange(minID.max(),maxID.max())):
	print i, mjd
	data,mjd = sb.allSkyDB(dateID)#  , sqlQ=sqlQ, dtypes=dtypes)
	starmap = healbin(data['ra'], data['dec'], data['dec'], 
	                  nside=nside, reduceFunc=np.size)
	good = np.where(starmap != hp.UNSEEN)[0]
	for i in good:
		print >>outfile, '%f, %i, %i' % (mjd, i, starmap[i])

	#progress = i/full_range*100
	#text = "\rprogress = %.3f%%"%progress
	#sys.stdout.write(text)
	#sys.stdout.flush()
	
print 'ack'
outfile.close()
