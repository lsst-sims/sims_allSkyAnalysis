import numpy as np
import healpy as hp
from medBD import medDB, single_frame
from lsst.sims.utils import Site
import ephem
import lsst.sims.skybrighntess as sb

# Loop through each frame, maybe fit a scaled background frame and a background model? Make some threshold and
# say what fraction of the sky is outside the threshold?

# Not sure what to do about airmass.  it looks like there's a good correlation within a night, but not clear if it varies night-to-night

mm = np.load('sky_maps.npz')
mm = mm['sky_maps'].copy()
