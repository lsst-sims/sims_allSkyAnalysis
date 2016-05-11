import numpy as np
import matplotlib.pylab as plt
from medDB import single_frame, medDB
import healpy as hp
from cloudy_stats import cloudyness

# Load up the statistics from running over all the nights and see how often it is cloudy

stats = np.load('movie_stats.npz')
cloudy_frac = stats['cloudy_frac'].copy()
umjd = stats['umjd'].copy().ravel()
sun_alts = stats['sun_alts'].copy().ravel()


good = np.where((cloudy_frac != -666) & (sun_alts[1:] < np.radians(-12.)))

plt.hist(cloudy_frac[good], bins=100)

# check a frame

# kinda cloudy? 
# kinda = np.where((cloudy_frac > .05) & (cloudy_frac < 0.1) & (sun_alts[1:] < np.radians(-12.)))[0]

kinda = np.where((cloudy_frac > .02) & (cloudy_frac < 0.05) & (sun_alts[1:] < np.radians(-12.)))[0]
clear = np.where((cloudy_frac < 0.02) & (sun_alts[1:] < np.radians(-12.)))[0]
cloudy = np.where((cloudy_frac >= 0.06) & (sun_alts[1:] < np.radians(-12.)))[0]

i = 501
f2 = single_frame(umjd[1:][kinda[i]])
f1 = single_frame(umjd[1:][kinda[i]-1])
unseen = np.where((f2 == hp.UNSEEN) | (f1 == hp.UNSEEN))
diff = (f2-f1)/f1
diff[unseen] = hp.UNSEEN
hp.mollview(diff, title=' %.2f' % cloudy_frac[kinda[i]], min=-.5, max=.5)

nObs = float(np.where((sun_alts[1:] < np.radians(-12.)))[0].size)

print 'clear = %f, kinda= %f, cloudy=%f' % (clear.size/nObs, kinda.size/nObs, cloudy.size/nObs)

cloudy_area, cloud_mask = cloudyness(diff)