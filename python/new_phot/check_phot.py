import numpy as np
import glob

# RA        Dec      HA      xcalc  ycalc    BT     VT     V    (B-V)  (V-I)     m      x       y    HAcalc Decalc  alt    rad   minst   dm   sky  dvig


def readAndSave():
    names = ['mjd','x','y','HA','Dec','HAcalc','Decalc','dAngle','Azi','Alt','ompix','xcalc','ycalc','r','theta','<dx>','<dy>','resid']
    types = [float]*len(names)

    files = glob.glob('ut*.txt')

    phot = np.empty(1, dtype = zip(names,types))
    for filename in files:
        tmp = np.genfromtxt(filename, dtype=phot.dtype)
        phot = np.concatenate((phot,tmp))

    np.savez('phot.npz', phot=phot[1:])


if __name__ == '__main__':

    #readAndSave()

    # So what am I going to do here?  I think compute airmasses for everything, then look at stellar mag vs airmass?
    # Fit a line (rejecting outliers) to all the stars
    # predict mags for each star the next night, see if I can do it, maybe bin by healpixels to see if I can get a transparency.
    # Looks like I need to match stars on my own, and compute RA's
