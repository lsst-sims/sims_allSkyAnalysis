import numpy as np
import healpy as hp


# Let's load up a few frames where there are some easy streaks of clouds that we want to try and dodge

if __name__ == '__main__':

    night_of_interest = 57413

    skyMaps = np.load('sky_maps.npz')
    umjd = skyMaps['umjd'].copy()
    good = np.where(np.floor(umjd) == night_of_interest)
    umjd = umjd[good]
    sun_alts = skyMaps['sunAlts'][good].copy()
    moon_alts = skyMaps['moonAlts'][good].copy()

    