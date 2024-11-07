

import csv
import pandas as pd
import numpy as np


for element in ["Ti1", "Ti2"]:
    lines = ["4533", "5193"]

    feh = []
    logg = []
    teff = []
    abundance = {}
    star_dict= {}
    df = pd.read_csv("Feros_abundances(1).csv")
    for column in df:
        # Means it's a star and not "nan" which is actualyl a float
        if (isinstance(df[column][0], str)):

            starn = df[column][0]
            teff = df[column][2]
            logg = df[column][3]
            feh = df[column][5]
            star_dict[starn] = {}

            if element == "Ti1":
                abund4533 = np.float32(df[column][31])
                abund5193 = np.float32(df[column][34])

                star_dict[starn][4533] = [teff, logg, feh, 1.5, abund4533]
                star_dict[starn][5193] = [teff, logg, feh, 1.5, abund5193]
            else:
                abund4443 = np.float32(df[column][44])
                abund4468 = np.float32(df[column][47])
                abund4501 = np.float32(df[column][48])
                abund4571 = np.float32(df[column][51])

                star_dict[starn][4443] = [teff, logg, feh, 1.5, abund4443]
                star_dict[starn][4468] = [teff, logg, feh,  1.5,abund4468]
                star_dict[starn][4501] = [teff, logg, feh,  1.5,abund4501]
                star_dict[starn][4571] = [teff, logg, feh,  1.5,abund4571]


    import pickle

    print([star_dict[i].keys() for i in star_dict.keys()])
    df = pd.read_csv("Magellan_Keck_abundances(1).csv")

    for column in df:

        # Means it's a star and not "nan" which is actualyl a float
        if (isinstance(df[column][0], str)):

            starn = df[column][0]
            teff = df[column][2]
            logg = df[column][3]
            feh = df[column][4]
            if starn not in star_dict:
                star_dict[starn] = {}
            if element == "Ti1":

                abund4533 = np.float32(df[column][31])
                abund5193 = np.float32(df[column][35])

                star_dict[starn][4533] = [teff, logg, feh,  1.5,abund4533]
                star_dict[starn][5193] = [teff, logg, feh,  1.5,abund5193]
            else:
                abund4443 = np.float32(df[column][46])
                abund4468 = np.float32(df[column][49])
                abund4501 = np.float32(df[column][50])
                abund4571 = np.float32(df[column][53])
                print(abund4571)
                star_dict[starn][4443] = [teff, logg, feh,  1.5,abund4443]
                star_dict[starn][4468] = [teff, logg, feh,  1.5,abund4468]
                star_dict[starn][4501] = [teff, logg, feh,  1.5,abund4501]
                star_dict[starn][4571] = [teff, logg, feh,  1.5,abund4571]

    import pickle
    print(len(star_dict))
    pickle.dump(star_dict, open(element+"yong_abundances.pkl", "wb"))


    print(star_dict)
