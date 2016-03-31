import numpy as np
import glob
from lsst.sims.utils import calcLmstLast, Site
from scipy.spatial import cKDTree as kdtree

from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u

# RA        Dec      HA      xcalc  ycalc    BT     VT     V    (B-V)  (V-I)     m      x       y    HAcalc Decalc  alt    rad   minst   dm   sky  dvig


def readAndSave():
    names = ['mjd','RA','dec','HA','dunno','HAcalc','Decalc','dAngle','Azi','Alt','ompix','xcalc','ycalc','r','theta','<dx>','<dy>','resid']
    types = [float]*len(names)

    files = glob.glob('ut*.txt')

    phot = np.empty(1, dtype = zip(names,types))
    for filename in files:
        tmp = np.genfromtxt(filename, dtype=phot.dtype)
        phot = np.concatenate((phot,tmp))

    np.savez('phot.npz', phot=phot[1:])


def treexyz(ra, dec):
    """Calculate x/y/z values for ra/dec points, ra/dec in radians."""
    # Note ra/dec can be arrays.
    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)
    return x, y, z




if __name__ == '__main__':

    #readAndSave()

    # So what am I going to do here?  I think compute airmasses for everything, then look at stellar mag vs airmass?
    # Fit a line (rejecting outliers) to all the stars
    # predict mags for each star the next night, see if I can do it, maybe bin by healpixels to see if I can get a transparency.
    # Looks like I need to match stars on my own, and compute RA's

    data = np.load('phot.npz')
    phot = data['phot'].copy()
    data.close()

    umjd = np.unique(phot['mjd'])
    phot.sort(order='mjd')
    left = np.searchsorted(phot['mjd'], umjd)
    right = np.searchsorted(phot['mjd'], umjd, side='right')

    nstars = right-left

    goodMJD = umjd[np.where(nstars == nstars.max())]

    # create a kdtree for all the ra,dec points
    x, y, z = treexyz(phot['RA'], phot['dec'])
    tree = kdtree(zip(x,y,z), leafsize=100)

    starIDs = np.zeros(phot.size, dtype=int)-666
    # Array to say if the star has matched multiple times
    multiFlag = np.zeros(phot.size, dtype=bool)



    radius = 300./3600. # return everything within 300 arcsec?
    x0, y0, z0 = (1, 0, 0)
    x1, y1, z1 = treexyz(np.radians(radius), 0)
    rad = np.sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)

    good = np.where(phot['mjd'] == goodMJD)
    onID = 1
    for i, blah in enumerate(phot[good]):
        indices = tree.query_ball_point((x[i], y[i], z[i]), rad)
        multiFlag[np.where(starIDs[indices] != -666)] = True
        starIDs[indices] = onID
        onID += 1







    # How to match up the data? Should I just go frame-by-frame?

    #starID = np.zeros(phot.size, dtype=int)-666
    #starID[left[0]:right[0]] = np.arange(starID[left[0]:right[0]].size)

    # Need to assign a star ID to each point.  Loop over each star in a frame that has a lot of stars, query the kdtree to see what's nearby and match it?




#    for i, mjd in enumerate(umjd[50:]):
#        c1 = SkyCoord(ra=phot['RA'][left[i-1]:right[i-1]]*u.degree,
#                      dec= phot['dec'][left[i-1]:right[i-1]]*u.degree)
#        c2 = SkyCoord(ra=phot['RA'][left[i]:right[i]]*u.degree,
#                      dec= phot['dec'][left[i]:right[i]]*u.degree)

#        idx, d2d, d3d = c1.match_to_catalog_sky(c2)
