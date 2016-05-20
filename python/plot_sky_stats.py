import numpy as np
import matplotlib.pylab as plt
from medDB import single_frame, medDB
import healpy as hp
from cloudy_stats import cloudyness
import ephem 
from lsst.sims.utils import Site

RdBu = plt.get_cmap('RdBu')
RdBu.set_bad('gray')
RdBu.set_under('w')

# Load up the statistics from running over all the nights and see how often it is cloudy

stats = np.load('movie_stats.npz')
cloudy_frac = stats['cloudy_frac'].copy()
umjd = stats['umjd'].copy().ravel()
sun_alts = stats['sun_alts'].copy().ravel()
rms_diff_frame = stats['rms_diff_frame'].copy().ravel()


good = np.where((cloudy_frac != -666) & (sun_alts[1:] < np.radians(-12.)) & (rms_diff_frame != hp.UNSEEN))

bins = np.arange(0, .21, .01).tolist()

# truncate large values
cloudy_frac[np.where(cloudy_frac > 0.2)] = 0.2

plt.hist(cloudy_frac[good], bins=bins, normed=True)
plt.xlabel('Cloudy Fraction')
plt.ylabel('Fraction of Observations')

plt.savefig('Plots/cloudy_hist.png')
plt.close('all')

# check a frame

# kinda cloudy? 
# kinda = np.where((cloudy_frac > .05) & (cloudy_frac < 0.1) & (sun_alts[1:] < np.radians(-12.)))[0]

# XXX--also make sure RMS is not hp.UNSEEN

fig = plt.figure(1)

clear_thesh = 0.02
cloudy_thresh = 0.05

kinda = np.where((cloudy_frac > clear_thesh) & (cloudy_frac < cloudy_thresh) & (sun_alts[1:] < np.radians(-12.)))[0]
clear = np.where((cloudy_frac < clear_thesh) & (sun_alts[1:] < np.radians(-12.)))[0]
cloudy = np.where((cloudy_frac >= cloudy_thresh) & (sun_alts[1:] < np.radians(-12.)))[0]

i = 501
f2 = single_frame(umjd[1:][kinda[i]])
f1 = single_frame(umjd[1:][kinda[i]-1])
unseen = np.where((f2 == hp.UNSEEN) | (f1 == hp.UNSEEN))
diff = (f2-f1)/f1
diff[unseen] = hp.UNSEEN
hp.mollview(diff, title='%.2f, Cloudy fraction %.2f' % (umjd[1:][kinda[i]], cloudy_frac[kinda[i]]),
            min=-.1, max=.1, cmap=RdBu, fig=1, sub=(2, 2, 1))


f2 = single_frame(umjd[1:][clear[i]])
f1 = single_frame(umjd[1:][clear[i]-1])
unseen = np.where((f2 == hp.UNSEEN) | (f1 == hp.UNSEEN))
diff = (f2-f1)/f1
diff[unseen] = hp.UNSEEN
hp.mollview(diff, title='%.2f, Cloudy fraction %.2f' % (umjd[1:][clear[i]], cloudy_frac[clear[i]]),
            min=-.1, max=.1, cmap=RdBu, fig=1, sub=(2, 2, 3))


f2 = single_frame(umjd[1:][cloudy[i]])
f1 = single_frame(umjd[1:][cloudy[i]-1])
unseen = np.where((f2 == hp.UNSEEN) | (f1 == hp.UNSEEN))
diff = (f2-f1)/f1
diff[unseen] = hp.UNSEEN
hp.mollview(diff, title='%.2f, Cloudy fraction %.2f' % (umjd[1:][cloudy[i]], cloudy_frac[cloudy[i]]),
            min=-.1, max=.1, cmap=RdBu, fig=1, sub=(2, 2, 4))


i = 800
f2 = single_frame(umjd[1:][kinda[i]])
f1 = single_frame(umjd[1:][kinda[i]-1])
unseen = np.where((f2 == hp.UNSEEN) | (f1 == hp.UNSEEN))
diff = (f2-f1)/f1
diff[unseen] = hp.UNSEEN
hp.mollview(diff, title='%.2f, Cloudy fraction %.2f' % (umjd[1:][kinda[i]], cloudy_frac[kinda[i]]),
            min=-.1, max=.1, cmap=RdBu, fig=1, sub=(2, 2, 2))

fig.savefig('Plots/cloudy_examples.png')


nObs = float(np.where((sun_alts[1:] < np.radians(-12.)) & (cloudy_frac != -666))[0].size)

print 'clear = %f, kinda= %f, cloudy=%f' % (clear.size/nObs, kinda.size/nObs, cloudy.size/nObs)


# Let's look at things per-night.  How many nights are totally clear, how many are totally clouded out, how many have clouds roll in.

# I've already clipped off the twilight time.
nights = np.zeros(umjd.size, dtype=int)
# maybe loop through and find all the sunrise times?
telescope = Site(name='LSST')
site = ephem.Observer()
site.lat = telescope.latitude_rad
site.lon = telescope.longitude_rad
site.elevation = telescope.height
doff = ephem.Date(0)-ephem.Date('1858/11/17')
udjd = umjd-doff

S = ephem.Sun()
site.date = udjd[0]
next_rise = site.next_rising(S)
current_night = 1

for i, djd in enumerate(udjd):
    site.date = djd
    nr = site.next_rising(S)
    if nr > next_rise + 0.1:
        current_night += 1
    nights[i] += current_night
    next_rise = nr

unights = np.unique(nights)

# make an array to describe how cloudy it is each night
names = ['clear', 'kinda', 'cloudy']
cloudy_nights = np.zeros(unights.size, dtype=zip(names, [float]*3))
for i, night in enumerate(unights[1:]):
    good = np.where(nights[1:] == night)[0]
    nframes = float(good.size)
    nclear = np.size(np.where(cloudy_frac[good] < clear_thesh)[0])
    ncloudy = np.size(np.where(cloudy_frac[good] > cloudy_thresh)[0])
    nkinda = nframes - nclear - ncloudy
    cloudy_nights[i]['clear'] = nclear/nframes
    cloudy_nights[i]['kinda'] = nkinda/nframes
    cloudy_nights[i]['cloudy'] = ncloudy/nframes

plt.close('all')

fig, ax = plt.subplots()
cloudy_nights.sort(order=['clear', 'cloudy', 'kinda'])
cloudy_nights = cloudy_nights[::-1]


ax.stackplot(np.arange(cloudy_nights.size)/float(cloudy_nights.size)*100, cloudy_nights['clear'],
             cloudy_nights['kinda'],
             cloudy_nights['cloudy'])


#for name in names:
#    ax.plot(unights, cloudy_nights[name], label=name)

#ax.legend()
ax.set_xlabel('%'+' of Nights')
ax.set_ylabel('Fraction')
fig.savefig('Plots/cloudy_nights.png')

nnights_all_clear = np.size(np.where(cloudy_nights['clear'] > 0.8)[0])
nnight_all_cloudy = np.size(np.where(cloudy_nights['cloudy'] > 0.8)[0])
nnights_want_to_dodge = cloudy_nights.size - nnights_all_clear - nnight_all_cloudy
print '    all_clear   all_cloudy    dodge'
print 'nights  %i    %i    %i' % (nnights_all_clear, nnight_all_cloudy, nnights_want_to_dodge)
print 'frac    %.2f  %.2f   %.2f ' % (nnights_all_clear/float(cloudy_nights.size),
                                      nnight_all_cloudy/float(cloudy_nights.size),
                                      nnights_want_to_dodge/float(cloudy_nights.size))

