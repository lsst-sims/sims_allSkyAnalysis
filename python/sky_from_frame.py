import numpy as np
import healpy as hp
import ephem

# Ideas on how to read in the all-sky frames and measure the sky values at various alt-az coords

# LSST site in radians
lsst_latitude = -0.52786436029017303
lsst_longitude = -1.2348099738104761
lsst_elev = 2650.0
# Stupid ephem uses Dublin Julian Date.
doff = ephem.Date(0)-ephem.Date('1858/11/17')


# Copying out some functions so this can be run outside the stack
def stupidFast_RaDec2AltAz(ra,dec,mjd):
    """
    Coordinate transformation is killing performance. Just use simple equations to speed it up
    and ignore abberation, precesion, nutation, nutrition, etc.
    """
    obs = ephem.Observer()
    obs.lon = lsst_longitude
    obs.lat = lsst_latitude
    obs.elevation = lsst_elev
    obs.date = mjd-doff
    lmst = obs.sidereal_time()
    ha = lmst-ra
    sindec = np.sin(dec)
    sinlat = np.sin(lat)
    coslat = np.cos(lat)
    alt = np.arcsin(sindec*sinlat+np.cos(dec)*coslat*np.cos(ha))
    az = np.arccos( (sindec-np.sin(alt)*sinlat)/(np.cos(alt)*coslat) )
    signflip = np.where(np.sin(ha) > 0)
    az[signflip] = 2.*np.pi-az[signflip]
    return alt,az


def haversine(long1, lat1, long2, lat2):
    """
    Return the angular distance between two points in radians

    @param [in] long1 is the longitude of point 1 in radians

    @param [in] lat1 is the latitude of point 1 in radians

    @param [in] long2 is the longitude of point 2 in radians

    @param [in] lat2 is the latitude of point 2 in radians

    @param [out] the angular separation between points 1 and 2 in radians

    From http://en.wikipedia.org/wiki/Haversine_formula
    """
    t1 = numpy.sin(lat2/2.-lat1/2.)**2
    t2 = numpy.cos(lat1)*numpy.cos(lat2)*numpy.sin(long2/2. - long1/2.)**2
    return 2*numpy.arcsin(numpy.sqrt(t1 + t2))


box_size = 10 # box size

nside = 32
# alt and az in radians
hpids = np.arange(hp.nside2npix(nside))
azs,alts = hp.pix2ang(nside,hpids)
alts = np.pi-alt

# only interested in things ~10 degrees above the horizon
alt_limit = np.radians(10.)

good = np.where(alts >= alt_limit)
hpids = hpids[good]
azs = azs[good]
alts = alts[good]

# Figure out the x-y points on the image that correspond to those alzs and alts
# XXX-make an array with the RA at every point in the image and an array with all the Decs. And the MJD of the frame

frame_alt, frame_az = stupidFast_RaDec2AltAz(ra,dec,mjd)
xs = []
ys = []
# This is the silly slow way to find the coordinates, but I'm too lazy to code up the kd-tree
for alt,az in zip(alts,azs):
    distances = haversine(frame_az, frame_alt, az, alt)
    good = np.where(distances == distances.min())[0]
    xs.append[good[0]]
    ys.append[good[1]]


# Print a header to the output file
f = open(filename, 'w')
print >>f, '# MJD, hpid, R_sky_mag, G_sky_mag, B_sky_mag ## nside=%i' % nside


for mjd in mjds:
    # Read in the RGB images for the given mjd
    images = [r_image, g_image, b_image]
        for x,y,hpid in zip(xs,ys,hpids):
            results = []
            for filterName, image in zip(['R','G','B'],images):
                # in case the point + boxsize goes off the frame
                try:
                    results.append(np.median(image[x-box_size/2:x+box_size/2,
                                                   y-box_size/2:y+box_size/2]))
                except:
                    results.append(0)
            # Convert to mags
            results = -2.5*np.log10(results)
            print f>> '%f, %i, %f, %f, %f' % (mjd, hpid, results[0],
                                              results[1], results[2])

f.close()
