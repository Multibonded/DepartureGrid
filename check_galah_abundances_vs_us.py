

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
full_galah_nlte = {}
element = "Ti1"
reduced = True
loadall = True
import pickle
if reduced:
    if loadall:
        full_galah_nlte = pickle.load(open(element + "_Galah_impacts2_Ati_reduced_all.pkl", "rb"))
    else:
        full_galah_nlte = pickle.load(open(element + "_Galah_impacts2_Ati_reduced.pkl", "rb"))

else:
    if loadall:
        full_galah_nlte = pickle.load(open(element + "_Galah_impacts2_Ati_all.pkl", "rb"))
    else:
        full_galah_nlte = pickle.load(open(element + "_Galah_impacts2_Ati.pkl", "rb"))


for starn in full_galah_nlte:
    for line in full_galah_nlte[starn]:
        lteab = full_galah_nlte[starn][line][-3] - full_galah_nlte[starn][line][2] - 4.9
        try:
            galah_ab = atom1['ind_Ti' + str(line) + "_fe"][starn]
        except KeyError:
            galah_ab = atom1['ind_Ti' + str(int(line)+1) + "_fe"][starn]
        diff = (round(lteab,4) - round(galah_ab, 4))
        if abs(diff) > 0.001:
            print(diff, round(lteab,4) - round(galah_ab, 4))
