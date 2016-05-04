import numpy as np
from utils import robustRMS
from scipy.stats import norm

# Make a function that calculates how cloudy the sky is based on a difference image



def cloudyness(diff_image, sigma_limit=None, airmass_map=None, airmass_limit=None):
	"""
	Parameters
	----------
	diff_image: np.array
	    An array that is the difference of two all-sky images

	sigma_limit: float
	    Assume the difference standard devation should not be above this
	"""

	skyRMS = robustRMS(diff_image)
	if sigma_limit is not None:
		skyRMS = np.min([skyRMS, sigma_limit])	

	im_x = np.sort((diff_image - np.median(diff_image)/skyRMS))
	im_cdf = np.arange(1, im_x.size+1, dtype=float)/(im_x.size)
	expected_cdf = norm.cdf(im_x)

	diff = np.abs(im_cdf - expected_cdf)

	# What is the fraction of pixels that are distributed very non-gaussian like (to within a factor of 2ish)
	return np.max(diff)




