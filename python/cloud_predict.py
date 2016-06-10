import numpy as np
import healpy as hp
from cloudy_stats import cloudyness
from medDB import single_frame
from lsst.sims.skybrightness import stupidFast_altAz2RaDec, stupidFast_RaDec2AltAz
from lsst.sims.utils import _raDec2Hpid, _hpid2RaDec, Site
from scipy import interpolate
import warnings

site = Site('LSST')
lat = site.latitude_rad
lon = site.longitude_rad

# Let's load up a few frames where there are some easy streaks of clouds that we want to try and dodge

def diff_hp(frame1,frame2, norm=True):
    """
    Take the difference of two healpy maps
    """
    mask = np.where( (frame1==hp.UNSEEN) | (frame2 == hp.UNSEEN))
    result = frame1 - frame2
    if norm:
        result = result/frame1
    result[mask] = hp.UNSEEN

    return result


def screen2hp(nside, mjd, npix=1000, height=500.):

    # generate a screen
    xx, yy = np.meshgrid(np.arange(-npix, npix, 1), np.arange(-npix, npix, 1), indexing='ij')
    r = (xx**2+yy**2)**0.5
    az = np.arctan2(yy,xx)
    alt = np.arctan(height/r)
    ra, dec = stupidFast_altAz2RaDec(alt, az, lat, lon, mjd)
    # Ah, I can convert these alt,az coords to ra,dec, then there's no problem using an ra,dec cloud map.
    hpids = _raDec2Hpid(nside, ra, dec)

    return hpids




def hp2screen(inmap, mjd, height=500, alt_limit=10., npts=600, grid_x=None, grid_y=None):
    """
    Convert a healpy map to a flat screen at height h
    """
    nside = hp.npix2nside(inmap.size)
    unmasked = np.where(inmap != hp.UNSEEN)[0]
    ra, dec = _hpid2RaDec(nside, unmasked)
    alt, az = stupidFast_RaDec2AltAz(ra, dec, lat, lon, mjd)
    good = np.where(np.degrees(alt) > alt_limit)
    r = height/np.tan(alt[good])
    x = r*np.cos(az[good])
    y = r*np.sin(az[good])
    z = inmap[unmasked][good]

    if grid_x is None:
        grid_x, grid_y = np.mgrid[x.min():x.max():600j, y.min():y.max():600j]
    screen_grid = interpolate.griddata(np.vstack((x, y)).T, z, (grid_x, grid_y), method='nearest')

    # Maybe compute the alt-az of the grid_x, grid_y so that it'll be fast to convert back to ra, dec? 
    # Pretty clear this should be a class rather than a bunch of functions.

    return grid_x, grid_y, screen_grid



# Maybe class this up, so that I can pre-compute the nsides and alt az conversions?
def propigateClouds(cloudmask, height, velocity, dt):
    """
    Given a cloud map, use the input height and velocity to propigate which pixels will be cloudy in the next frame
    """

    higher_cm = hp.pixelfunc.ud_grade(cloudmask, nside_out=nside*2)
    cloud_blobs = np.where(higher_cm != 0)
    # find the alt and az of each cloud blob

    # get x,y,z coords of each blob

    # 
    xnew = x*velocity[0]*dt
    ynew = y*velocity[1]*dt

    

def predictCloudMotion(cloudMask, diffTime):
    """
    Given a cloud mask, try to fit an altitude and velocity to all the cloudy pixels and 
    propigate it forward in time to the next cloudy frame

    Inputs
    ------
    cloudMask : healpy mask
        Mask with +/-1 and 0 values.  Assumes positive pixels propigate to negative pixels
    diffTime : float
        The time gap between exposures.  Seconds.
    """



class CloudMotion(object):
    """
    Take two maps and predict where the clouds will be in the future
    """

    def __init__(self, nside=64):
        """
        Set up the resolution of various things?
        """

        self.predictLimit = 20.  # minutes

    def set_maps(self, map1,map2, mjd1,mjd2):
        """
        set the two cloud maps to use to generate cloud motion predictions
        """

        # Put things in correct order
        if mjd2 < mjd1:
            self.mjd1 = mjd2
            self.mjd2 = mjd1
            self.map2 = map1
            self.map1 = map2
        else:
            self.mjd1 = mjd1
            self.mjd2 = mjd2
            self.map2 = map2
            self.map1 = map1

        self.dt = self.mjd2-self.mjd1

        # Diff the maps, project to a screen.
        # hmm, should probably change this to be a diff on screen rather than on sphere. 
        diff = diff_hp(single_frame(umjd[startnum+1]), single_frame(umjd[startnum]))
        cloudyArea, cloudmask = cloudyness(diff)

        # Find the best shift, dx, dy


    def forecast(self, mjd, usemap='both'):
        """
        Return a healpix map that is an ra,dec cloud mask prediction at input mjd

        mjd : float
            The Modified Julian Date to predict the cloud cover too
        usemap : str (map1, map2, both)
            Which input map should be propigated forward in time.  Default is to be conservative 
            and propigate all pixels flagged as cloudy in both maps.
        """

        dt2 = (mjd-self.mjd2)*24*60
        if  dt2 > self.predictLimit:
            warnings.warn('Time gap of %f minutes is greater than limit %f minutes' % (dt2, self.predictLimit))





        # What is the best way to test the forecast? Maybe say what fraction of sky is cloudy vs what fraction of 
        # predicted clear space is cloudy?  Above some airmass limit.  Also a stat on how much clear sky was lost
        # to falsly predicted cloudy area.


if __name__ == '__main__':

    night_of_interest = 57413

    skyMaps = np.load('sky_maps.npz')
    umjd = skyMaps['umjd'].copy()
    good = np.where((np.floor(umjd) == night_of_interest) & (skyMaps['sunAlts'] < np.radians(-12)))
    umjd = umjd[good]
    sun_alts = skyMaps['sunAlts'][good].copy()
    moon_alts = skyMaps['moonAlts'][good].copy()

    startnum = 300

    diff = diff_hp(single_frame(umjd[startnum+1]), single_frame(umjd[startnum]))
    cloudyArea, cloudmask = cloudyness(diff)
    # ah, cloudmask is in ra,dec.  Need to rotate os that we're in alt-az!

    nside = hp.npix2nside(cloudmask.size)
    #hpids = screen2hp(nside)
    #screen = cloudmask[hpids]

    grid_x, grid_y, screen_grid = hp2screen(cloudmask, umjd[startnum])



    #H, xe, ye = np.histogram2d(x,y, bins=1000, weights=z)

    # OK, I can interpolate this to a meshgrid, then cross-correlate the plus and minus.  
    #screenInterp = interpolate.interp2d(x, y, z, kind='linear')

    #npix = 3000
    #xx, yy = np.meshgrid(np.arange(-npix, npix, 1), np.arange(-npix, npix, 1), indexing='ij')

    #zz = screenInterp(xx,yy)


    # Ah, this is how to change resolution on maps
    # high = hp.pixelfunc.ud_grade(cloudmask, nside_out=nside*2)
    # lower = hp.pixelfunc.ud_grade(high, nside_out=nside)


    # General plan:  flag a bunch of pixels as cloudy.  subsample to higher res.  find the best fit altitude and vel that moves cloudy pixels around.