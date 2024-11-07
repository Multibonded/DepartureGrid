import os
import pickle
import sys
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from scipy.interpolate import interpn
from scipy.interpolate import interp1d
import warnings

from astropy.io import fits

atomname = r"/home/jama2357/Documents/TiFeAtoms/DepartureGrid/GALAH_DR3_main_allspec_v2.fits"
atom = fits.open(atomname)
atom1 = atom[1].data
print(repr(atom[1].header))

element = "Ti1"
reduced = True
all = True



feh = atom1['fe_h']
logg = atom1['logg']
teff = atom1['teff']
vmic = atom1['vmic']
# Indexes of stars close to our sun's parameters
twindex = np.where(
    np.logical_and(
    np.logical_and(
    np.logical_and(
        abs(teff-5772) <= 60,
        abs(logg-4.438) <= 0.05),
        abs(vmic-1) <= 0.2),
        abs(feh) <= 0.05)

)[0]
for i in twindex:

    if np.isnan(atom1[i]['ind_Ti4849_fe']):
        continue
    else:
        print(atom1[i]['ind_Ti4849_fe'])
exit()
# we have indexes of solar twins,n ot the star number which is whbat we use in the dictioanry so now we convert
