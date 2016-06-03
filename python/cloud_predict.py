import numpy as np
import healpy as hp
from cloudy_stats import cloudyness
from medDB import single_frame
import healpy as hp

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


def screen2hp(nside, npix=1000, height=500.):

    # generate a screen
    xx, yy = np.meshgrid(np.arange(-npix, npix, 1), np.arange(-npix, npix, 1), indexing='ij')
    r = (xx**2+yy**2)**0.5
    az = np.arctan2(yy,xx)
    alt = np.arctan(height/r)

    # Ah, I can convert these alt,az coords to ra,dec, then there's no problem using an ra,dec cloud map.

    hpids = hp.ang2pix(nside, alt, az)
    return hpids




def hp2screen(nside, hpInd, height):
    """
    Convert a healpy map to a flat screen
    """
    
    # XXX--check that this is right!!
    alt, az = hp.pix2ang(nside, hpInd)
    radius = height/np.tan(alt)
    


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
    hpids = screen2hp(nside)

    screen = cloudmask[hpids]


    high = hp.pixelfunc.ud_grade(cloudmask, nside_out=nside*2)
    lower = hp.pixelfunc.ud_grade(high, nside_out=nside)


    # General plan:  flag a bunch of pixels as cloudy.  subsample to higher res.  find the best fit altitude and vel that moves cloudy pixels around.