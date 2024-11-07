"""We find the indexes for the most reliable stars in Galah. We also make sure to add all reasonable stars below -1.5 metallicity due to a very low number of them appearing in the MOST reliable ones."""
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
"""
print(atom1['snr_c1_iraf'][10000])
print(atom1['snr_c4_iraf'][10000])
a=np.argsort(atom1['snr_c1_iraf'])[:10000]
b=np.argsort(atom1['snr_c2_iraf'])[:10000]
c=np.argsort(atom1['snr_c3_iraf'])[:10000]
d=np.argsort(atom1['snr_c4_iraf'])[:10000]
print(a,b,c,d)
print(len(a))
e = (set(np.concatenate((a, b, c, d))))
print(len(e))"""

element = "Ti1"
for element in ["Ti1", "Ti2"]:
    # Allowing all reasonable low metallicity stars in
    if element =="Ti2":
        tiflag = "flag_"+element+"_fe"
    elif element == "Ti1":
        tiflag = "flag_ti_fe"
    print(atom1['flag_ti2_fe'])
    # Allowing all reasonable low metallicity stars in
    alow = np.where(


                                            np.logical_and(
        np.logical_and(
        np.logical_and(
        np.logical_and(
        np.logical_and(atom1['flag_sp']==0, atom1['flag_fe_h'] == 0)
            , atom1['fe_h'] <= -2.5)
            , atom1['snr_c3_iraf'] > 30)
            , atom1[tiflag] == 0),
    atom1['flag_ti2_fe'] <= 0.5),

    )[0]

    print(atom1['snr_c3_iraf'][alow])


    print(alow, len(alow))
    # We find indexes of reliable stars
    a = np.where(


        np.logical_and(
        np.logical_and(
        np.logical_and(
            np.logical_and(
                np.logical_and(
                atom1['flag_sp']==0,
                atom1['flag_fe_h'] == 0),
            atom1[tiflag] == 0),
        atom1['snr_c3_iraf'] > 30),
        atom1['fe_h'] <= 0.5),
            atom1['flag_ti2_fe'] <= 0.5),

    )[0]
    # Limit the numbero f stars to a vague number which is reduced further to match the indexes in the most reliable stars
    # that considers more than just snrc3
    accurate_params = np.argsort(atom1['snr_c3_iraf'])[-102000:]
    print("H", accurate_params)
    accurate_params = list(np.where(atom1['snr_c3_iraf'] > 70)[0])
    print(accurate_params)
    acc_final = [i for i in accurate_params if i in a]
    for lowindex in alow:
        if lowindex not in acc_final:
            acc_final.append(lowindex)

    print("Number of final reliable stars:", len(acc_final), "for element", element)

    print("min", min(atom1['snr_c3_iraf'][accurate_params]))

    pickle.dump(acc_final, open(element+"Accurate_param_indexes_galah.pkl", "wb"))
