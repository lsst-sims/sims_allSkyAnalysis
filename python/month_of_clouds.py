import numpy as np
import healpy as hp
from medDB import medDB, single_frame
import sqlalchemy as sqla
import matplotlib.pylab as plt
# import lsst.sims.skybrightness_pre as sb
import lsst.sims.skybrightness as sb
from cloudy_stats import cloudyness
import lsst.sims.utils as utils
import sys

# Make a month of cloud maps


def gen_cloud_maps():
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
        mjd_clean = [ack[0] for ack in data]
        mjd_clean = np.array(mjd_clean)

    # Let's go 57200 as a start
    mjd_start = 57200
    mjd_length = 40.

    query = 'select distinct(mjd) from medskybrightness where mjd > %f and mjd < %f;' % (mjd_start,
                                                                                         mjd_start+mjd_length)

    cursor.execute(query)
    data = cursor.fetchall()
    mjd_clean = np.array([ack[0] for ack in data])

    cloud_maps = []
    mjd_final = []

    moon_limit = np.radians(30.)
    am_limit = 2.5

    maxI = mjd_clean[1:].size
    for i, mjd in enumerate(mjd_clean[1:]):
        if (mjd_clean[i]-mjd_clean[i-1]) < (2./60./24.):
            sf = single_frame(mjd_clean[i-1])
            sf2 = single_frame(mjd_clean[i])
            if (sf.max() != hp.UNSEEN) & (sf2.max() != hp.UNSEEN):
                diff = (sf-sf2)/sf
                mask = np.where((sf == hp.UNSEEN) | (sf2 == hp.UNSEEN))
                diff[mask] = hp.UNSEEN
                if diff.max() != hp.UNSEEN:
                    cf, cloud_mask = cloudyness(diff, skyRMS_max=0.05)
                    sm.setRaDecMjd(ra, dec, mjd)
                    skymap = sm.returnMags()['r']
                    mask = np.where((skymap == hp.UNSEEN) | np.isnan(skymap))[0]
                    cloud_mask[mask] = 2
                    moon_dist = utils.haversine(ra, dec, sm.moonRA, sm.moonDec)
                    mask = np.where((moon_dist < moon_limit) | (sm.airmass < 1.) | (sm.airmass > am_limit))
                    cloud_mask[mask] = 2
                    cloud_maps.append(cloud_mask)
                    # put in some masking near the moon
                    mjd_final.append(mjd)

                    progress = i/float(maxI)*100
                    text = "\r progress = %.1f%%" % progress
                    sys.stdout.write(text)
                    sys.stdout.flush()

    cloud_maps = np.array(cloud_maps)
    map_key = {'good': 0, 'masked': 2, 'bright_in_current': -1, 'bright_in_last': 1}
    np.savez('month_o_clouds.npz', cloud_maps=cloud_maps, mjds=mjd_final, map_key=map_key)


if __name__ == '__main__':

    gen_cloud_maps()
