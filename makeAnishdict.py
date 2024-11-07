import pickle

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

anish = {}

ti1 = pickle.load(open("Ti1_Galah_impacts.pkl", "rb"))
ti2 = pickle.load(open("Ti2_Galah_impacts.pkl", "rb"))
counter = 0
for i in (ti1.keys()):
    if counter % 500 == 0:
        print(counter, "/", len(ti1.keys()))
    counter += 1
    if i not in anish:
        anish[i] = {}
    teff = atom1[i]['teff']
    e_teff = atom1[i]['e_teff']
    logg = atom1[i]['logg']
    e_logg = atom1[i]['e_logg']
    feh = atom1[i]['fe_h']
    e_feh = atom1[i]['e_fe_h']
    vmic = atom1[i]['vmic']
    c1 = atom1[i]['snr_c1_iraf']
    c2 = atom1[i]['snr_c2_iraf']
    c3 = atom1[i]['snr_c3_iraf']

    c4 = atom1[i]['snr_c4_iraf']
    anish[i]['params'] = [teff, e_teff, logg, e_logg, feh, e_feh, vmic, c1, c2, c3, c4]
    for wl in ti1[i].keys():
        anish[i][wl] = [ti1[i][wl][-3], ti1[i][wl][-2]]
    if i in ti2:
        for wl in ti2[i].keys():
            anish[i][wl] = [ti2[i][wl][-3], ti2[i][wl][-2]]
pickle.dump(anish, open("AnishGalahDict.pkl", "wb"))
print(anish)