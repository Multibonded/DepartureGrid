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
def sundex():
    for element in ["Ti1", "Ti2"]:
        # Allowing all reasonable low metallicity stars in
        if element =="Ti2":
            tiflag = "flag_"+element+"_fe"
            lines =         [4719, 4798, 4874]

        elif element == "Ti1":
            lines =         [4758, 4759, 5689, 5739]

            tiflag = "flag_ti_fe"

        suns = np.where(

            np.logical_and(

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
                atom1['snr_c3_iraf'] > 30),

        )[0]

        suns = [i for i in suns]
        print(type(suns[0]))

        print("Number of final reliable solar twins:", len(suns), "for element", element)
        exit()
        for line_ab in lines:
            linenum = line_ab
            if str(line_ab) == "4781":
                line_ab = atom1["ind_Ti4782_fe"]
            elif str(line_ab) == "4797":
                line_ab = atom1["ind_Ti4798_fe"]
            elif str(line_ab) == "5739_single":
                line_ab = atom1["ind_Ti5739_fe"]
            elif str(line_ab) == "4801":
                line_ab = atom1["ind_Ti4802_fe"]
            elif str(line_ab) == "6716":
                line_ab = atom1["ind_Ti6717_fe"]
            elif str(line_ab) == "7852":
                line_ab = atom1["ind_Ti7853_fe"]

            elif str(line_ab) == "4865":
                line_ab = atom1["ind_Ti4866_fe"]
            elif str(line_ab) == "4798":
                line_ab = atom1["ind_Ti4799_fe"]
            elif str(line_ab) == "4764":
                line_ab = atom1["ind_Ti4765_fe"]
            elif str(line_ab) == "4719":
                line_ab = atom1["ind_Ti4720_fe"]
            else:
                line_ab = atom1['ind_Ti' + str(line_ab) + "_fe"]
            avrg = np.nanmean(line_ab[suns])
            print("GALAH:", linenum, avrg)

        b = pkl.load(open(element + "_Galah_impacts2_Ati_reduced_all.pkl", "rb"))
        c = pkl.load(open(element + "_Galah_impacts_lineless2_Ati_reduced_all.pkl", "rb"))


        feh = np.asarray([i[2] for i in c.values()])
        logg = np.asarray([i[1] for i in c.values()])
        teff = np.asarray([i[0] for i in c.values()])
        vmic = np.asarray([i[3] for i in c.values()])
        # Indexes of stars close to our sun's parameters
        suns = np.where(
            np.logical_and(
            np.logical_and(
            np.logical_and(
                abs(teff-5772) <= 60,
                abs(logg-4.438) <= 0.05),
                abs(vmic-1) <= 0.2),
                abs(feh) <= 0.05)

        )[0]
        for line_ab in lines:
            linenum = line_ab
            if str(line_ab) == "4781":
                line_ab = atom1["ind_Ti4782_fe"]
            elif str(line_ab) == "4797":
                line_ab = atom1["ind_Ti4798_fe"]
            elif str(line_ab) == "5739_single":
                line_ab = atom1["ind_Ti5739_fe"]
            elif str(line_ab) == "4801":
                line_ab = atom1["ind_Ti4802_fe"]
            elif str(line_ab) == "6716":
                line_ab = atom1["ind_Ti6717_fe"]
            elif str(line_ab) == "7852":
                line_ab = atom1["ind_Ti7853_fe"]

            elif str(line_ab) == "4865":
                line_ab = atom1["ind_Ti4866_fe"]
            elif str(line_ab) == "4798":
                line_ab = atom1["ind_Ti4799_fe"]
            elif str(line_ab) == "4764":
                line_ab = atom1["ind_Ti4765_fe"]
            elif str(line_ab) == "4719":
                line_ab = atom1["ind_Ti4720_fe"]
            else:
                line_ab = atom1['ind_Ti' + str(line_ab) + "_fe"]
            avrg = np.nanmean(line_ab[suns])
            print("OURS:", linenum, avrg)

for element in ["Ti1", "Ti2"]:
    # Allowing all reasonable low metallicity stars in
    if element == "Ti2":
        tiflag = "flag_" + element + "_fe"
        lines = [4719, 4798, 4874]

    elif element == "Ti1":
        lines = [4758, 4759, 5689, 5739]

        tiflag = "flag_ti_fe"


reliablestars = np.where(


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
print("argsort")
accurate_params = np.argsort(atom1['snr_c3_iraf'])[-102000:]
print("listing")
accurate_params = list(np.where(atom1['snr_c3_iraf'] > 70)[0])
print("i for i")
acc_final = [i for i in accurate_params if i in reliablestars]
print(len(acc_final))
print("min:", max(atom1['snr_c2_iraf'][acc_final]))
accurate = pickle.load(open(element + "Accurate_param_indexes_galah.pkl", "rb"))

print("on to suns")
suns = np.where(


        np.logical_and(
            np.logical_and(
                np.logical_and(
                    np.logical_and(
                        np.logical_and(
                            np.logical_and(
                                atom1['r_est_dr2'][acc_final] < 500,
                                atom1['flag_sp'][acc_final] == 0),
                            atom1[tiflag][acc_final] == 0),
                        atom1['snr_c2_iraf'][acc_final] > 30),
                    abs(4.44 - atom1['logg'][acc_final]) <= 0.1),
                abs(5773 - atom1['teff'][acc_final]) <= 100),
            abs(0 - atom1['fe_h'][acc_final]) <= 0.1),

)[0]
print(len(suns))

b = pkl.load(open(element + "_Galah_impacts2_Ati_reduced_all.pkl", "rb"))
c = pkl.load(open(element + "_Galah_impacts_lineless2_Ati_reduced_all.pkl", "rb"))

keylist = []
index_list = (np.asarray(suns))

# lists to contain the lkte and nlte titanium abundances (lte from galah, nlte from us) for solar twins for us to
# use line by line
lte = {}
nlte = {}
mean = {}
errors = {}
ltesquare = {}

if "1" in element:
    lines = ["4758", "4759", "5689", "5739"]
else:
    lines = ["4798", "4719", "4874"]

for starn in index_list:
    teffstar = atom1[starn]['teff']
    fehstar = atom1['fe_h'][starn]
    gstar = atom1['logg'][starn]
    for line in lines:
        """if int(line) not in linekeep:
            continue"""
        # print("b", b[starn])
        if line in lte:
            lte[line].append(b[starn][line][-3])
            nlte[line].append(b[starn][line][-2])
            errors[line].append(b[starn][line][-4])
            if line == "5739_single":
                ltesquare[line].append(atom1[starn]['ind_Ti5739_fe'] + atom1[starn]['fe_h'] + 4.9)
            else:
                key = "ind_Ti" + str(line) + "_fe"
                key2 = "ind_Ti" + str(int(line) + 1) + "_fe"
                try:
                    ltesquare[line].append(atom1[starn][key] + atom1[starn]['fe_h'] + 4.9)
                except KeyError:
                    ltesquare[line].append(atom1[starn][key2] + atom1[starn]['fe_h'] + 4.9)

        else:
            lte[line] = [b[starn][line][-3]]
            nlte[line] = [b[starn][line][-2]]
            errors[line] = [b[starn][line][-4]]
            if line == "5739_single":
                ltesquare[line] = ([atom1[starn]['ind_Ti5739_fe'] + atom1[starn]['fe_h'] + 4.9])
            else:
                try:
                    ltesquare[line] = [atom1[starn]['ind_Ti' + str(line) + '_fe'] + atom1[starn]['fe_h'] + 4.9]
                except KeyError:
                    ltesquare[line] = [atom1[starn]['ind_Ti' + str(int(line) + 1) + '_fe'] + atom1[starn]['fe_h'] + 4.9]

total = []
totalnlte = []
totalerrors = []
lenlist = []
meanerrors = []
for line in lte:
    total.extend(lte[line])

    print("length of ", line, ":", len(lte[line]))
    lenlist.append(len(lte[line]))
    totalnlte.extend(nlte[line])
    totalerrors.extend(errors[line])
    meanerrors.append(np.mean(errors[line]))
    print("std", line, np.std(lte[line]))

    if len(lte[line]) > 2:
        mean[line] = [np.average(lte[line], weights=errors[line]), np.average(nlte[line], weights=errors[line])]
        meangalah = np.average(ltesquare[line], weights=errors[line])
        print("square", line, meangalah)

    else:
        mean[line] = [4.9, 4.9]
    # print(line, mean[line])
print(total)
print(
    "\nLTE and NLTE averaged from EACH line in EACH star, weighted as usual (But remember there are many more 4758/4759 lines and are thus biased)"
    , np.average(np.asarray(total), weights=np.asarray(totalerrors)),
    np.average(np.asarray(totalnlte), weights=np.asarray(totalerrors)))
print("unweighted total avrg", np.mean(total), np.mean(totalnlte))
print("error total:", np.std(total), np.std(totalnlte))
lte = (np.average([i[0] for i in mean.values()], weights=meanerrors))
nlte = (np.average([i[1] for i in mean.values()], weights=meanerrors))
print("\n Following are lte and nlte values for averaged of each line, weighted by each line's mean error:")
print("LTE:", lte, "NLTE:", nlte)
print("LTE errors:", np.std([i[0] for i in mean.values()]))
print("NLTE errors:", np.std([i[1] for i in mean.values()]))
print("AVerage based on number of lines:", np.average(lte))
corrections = (nlte - lte)
print("Corrections", corrections)
print("mean", mean)
