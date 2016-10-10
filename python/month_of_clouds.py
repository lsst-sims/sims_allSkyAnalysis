import numpy as np
import healpy as hp
from medDB import medDB, single_frame
import sqlalchemy as sqla
import matplotlib.pylab as plt
# import lsst.sims.skybrightness_pre as sb
import lsst.sims.skybrightness as sb
from cloudy_stats import cloudyness
import lsst.sims.utils as utils

# Make a month of cloud maps

sm = sb.SkyModel(mags=True)
hpindx = np.arange(hp.nside2npix(32))
ra, dec = utils._hpid2RaDec(32, hpindx)


do_query = False
dbAddress = 'sqlite:///meddata.sqlite'
engine = sqla.create_engine(dbAddress)
connection = engine.raw_connection()
cursor = connection.cursor()

# Check all the mjd values in the db
if do_query:
    query = 'select distinct(mjd) from medskybrightness;'
    cursor.execute(query)
    data = cursor.fetchall()
    data_clean = [ack[0] for ack in data]
    data_clean = np.array(data_clean)


# Let's go 57200 as a start
mjd_start = 57200
mjd_length = 40.

query = 'select distinct(mjd) from medskybrightness where mjd > %f and mjd < %f;' % (mjd_start,
                                                                                     mjd_start+mjd_length)

cursor.execute(query)
data = cursor.fetchall()
data_clean = np.array([ack[0] for ack in data])

sf = single_frame(data_clean[0])
sf2 = single_frame(data_clean[1])
diff = (sf-sf2)/sf
mask = np.where((sf == hp.UNSEEN) | (sf2 == hp.UNSEEN))

diff[mask] = hp.UNSEEN

cloud_maps = []
mjd_final = []

for i, mjd in enumerate(data_clean[1:10]):
    if (data_clean[i]-data_clean[i-1]) < (2./60./24.):
        sf = single_frame(data_clean[0])
        sf2 = single_frame(data_clean[1])
        diff = (sf-sf2)/sf
        mask = np.where((sf == hp.UNSEEN) | (sf2 == hp.UNSEEN))
        diff[mask] = hp.UNSEEN
        cf, cloud_mask = cloudyness(diff, skyRMS_max=0.05)
        sm.setRaDecMjd(ra, dec, mjd)
        skymap = sm.returnMags()['r']
        mask = np.where((skymap == hp.UNSEEN) | np.isnan(skymap))[0]
        cloud_mask[mask] = 2
        cloud_maps.append(cloud_mask)
        # put in some masking near the moon
        mjd_final.append(mjd)
