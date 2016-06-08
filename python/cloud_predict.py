import numpy as np
import healpy as hp
from cloudy_stats import cloudyness
from medDB import single_frame
from lsst.sims.skybrightness import stupidFast_altAz2RaDec, stupidFast_RaDec2AltAz
from lsst.sims.utils import _raDec2Hpid, _hpid2RaDec, Site
from scipy import interpolate

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




def hp2screen(inmap, mjd, height=500, alt_limit=10.):
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
    return x,y,z



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

    x, y, z = hp2screen(cloudmask, umjd[startnum])

    H, xe, ye = np.histogram2d(x,y, bins=1000, weights=z)

    # OK, I can interpolate this to a meshgrid, then cross-correlate the plus and minus.  
    #screenInterp = interpolate.interp2d(x, y, z, kind='linear')

    #npix = 3000
    #xx, yy = np.meshgrid(np.arange(-npix, npix, 1), np.arange(-npix, npix, 1), indexing='ij')

    #zz = screenInterp(xx,yy)


    # Ah, this is how to change resolution on maps
    # high = hp.pixelfunc.ud_grade(cloudmask, nside_out=nside*2)
    # lower = hp.pixelfunc.ud_grade(high, nside_out=nside)


    # General plan:  flag a bunch of pixels as cloudy.  subsample to higher res.  find the best fit altitude and vel that moves cloudy pixels around.