
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




element = "Ti1"


reliable     = pkl.load(open(element+"Accurate_param_indexes_galah.pkl", "rb"))
atomname = r"/home/jama2357/Documents/TiFeAtoms/DepartureGrid/GALAH_DR3_main_allspec_v2.fits"
atom = fits.open(atomname)
atom1 = atom[1].data
all = False
if all:
    lines = ["4758", "4759", "4778", "4781", "4797", "4801", "4820", "5689", "5716", "5720", "5739", "5866", "6716",
             "7852"]
else:
    lines = ["4758", "4759", "5689", "5716", "5739", "6599"]
    lines = ["4758", "4759", "5689", "5739"]
    # lines = ['4758','4759','4781','4801','4820','5739'] # galah ones
ti_ab = {}

for line_ab in lines:
    if line_ab == "4781":
        ti_ab[(line_ab)] = atom1["ind_Ti4782_fe"]
    elif line_ab == "4797":
        ti_ab[(line_ab)] = atom1["ind_Ti4798_fe"]
    elif line_ab == "5739_single":
        ti_ab[(line_ab)] = atom1["ind_Ti5739_fe"]
    elif line_ab == "4801":
        ti_ab[(line_ab)] = atom1["ind_Ti4802_fe"]
    elif line_ab == "6716":
        ti_ab[(line_ab)] = atom1["ind_Ti6717_fe"]
    elif line_ab == "7852":
        ti_ab[(line_ab)] = atom1["ind_Ti7853_fe"]
    else:
        ti_ab[(line_ab)] = atom1['ind_Ti' + str(line_ab) + "_fe"]
ti_ab["GALAHavrg"] = atom1['ti_fe']
if all:
    lines = ["4874", "4865", "4849", "4798", "4764", "4719"]

else:
    lines = ["4798", "4719", "4874"]
ti_ab2 = {}

for line_ab in lines:
    if line_ab == "4865":
        ti_ab2[(line_ab)] = atom1["ind_Ti4866_fe"]
    elif line_ab == "4798":
        ti_ab2[(line_ab)] = atom1["ind_Ti4799_fe"]
    elif line_ab == "4764":
        ti_ab2[(line_ab)] = atom1["ind_Ti4765_fe"]
    elif line_ab == "4719":
        ti_ab2[(line_ab)] = atom1["ind_Ti4720_fe"]
    else:
        ti_ab2[(line_ab)] = atom1['ind_Ti' + str(line_ab) + "_fe"]
    ti_ab["GALAHavrg"] = atom1['ti2_fe']
if element =="Ti2":
    tiflag = "flag_"+element+"_fe"
    lines =         [4719, 4798, 4874]

elif element == "Ti1":
    lines =         [4758, 4759, 5689, 5739]

    tiflag = "flag_ti_fe"

sunstars = np.where(

    np.logical_and(

    np.logical_and(
    np.logical_and(
    np.logical_and(
        np.logical_and(
            np.logical_and(
            atom1['r_est_dr2']<500,
            atom1['flag_sp'] == 0),
        atom1[tiflag] == 0),
    atom1['snr_c2_iraf'] > 30),
    abs(4.44 - atom1['logg']) <= 0.1),
        abs(5773-atom1['teff']) <= 100),
        abs(0 - atom1['fe_h']) <= 0.1),

)[0]

a = np.where(

    np.logical_and(
        np.logical_and(
            np.logical_and(
                np.logical_and(
                    np.logical_and(
                        atom1['flag_sp'] == 0,
                        atom1['flag_fe_h'] == 0),
                    atom1[tiflag] == 0),
                atom1['snr_c3_iraf'] > 30),
            atom1['fe_h'] <= 0.5),
        atom1['flag_ti2_fe'] <= 0.5),

)[0]


almostall = np.where(atom1['flag_ti_fe'] == 0)

for ab in ti_ab:
    print(ab)
    if ab == "GALAHavrg":
        plt.scatter(atom1['fe_h'][almostall], ti_ab[ab][almostall]-ti_ab2[ab][almostall], s=2)
        plt.title(ab)
        plt.show()

