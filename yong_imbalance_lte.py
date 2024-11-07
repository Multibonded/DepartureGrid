

import csv
import pandas as pd
import numpy as np

star_dict = {}
star_dict2 = {}

for element in ["Ti1", "Ti2"]:
    lines = ["4533", "5193"]

    feh = []
    logg = []
    teff = []
    abundance = {}
    df = pd.read_csv("Feros_abundances(1).csv")
    for column in df:
        # Means it's a star and not "nan" which is actualyl a float
        if (isinstance(df[column][0], str)):

            starn = df[column][0]
            teff = df[column][2]
            logg = df[column][3]
            feh = df[column][5]
            if starn not in star_dict:
                star_dict[starn] = {}
            if starn not in star_dict2:
                star_dict2[starn] = {}


            if element == "Ti1":
                abund4533 = np.float32(df[column][31])
                abund5193 = np.float32(df[column][34])
                if abund4533 != -666.0:

                    star_dict[starn][4533] = [teff, logg, feh, 1, abund4533]
                if abund5193 != -666.0:

                    star_dict[starn][5193] = [teff, logg, feh, 1, abund5193]
            else:
                abund4443 = np.float32(df[column][44])
                abund4468 = np.float32(df[column][47])
                abund4501 = np.float32(df[column][48])
                abund4571 = np.float32(df[column][51])

                if abund4443 != -666.0:
                    star_dict2[starn][4443] = [teff, logg, feh, 1, abund4443]
                if abund4468 != -666.0:

                    star_dict2[starn][4468] = [teff, logg, feh,  1,abund4468]
                if abund4501 != -666.0:


                    star_dict2[starn][4501] = [teff, logg, feh,  1,abund4501]
                if abund4571 != -666.0:

                    star_dict2[starn][4571] = [teff, logg, feh,  1,abund4571]


    import pickle

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
            if starn not in star_dict2:
                star_dict2[starn] = {}
            if element == "Ti1":

                abund4533 = np.float32(df[column][31])
                abund5193 = np.float32(df[column][35])

                if abund4533 != -666.0:

                    star_dict[starn][4533] = [teff, logg, feh, 1, abund4533]
                if abund5193 != -666.0:

                    star_dict[starn][5193] = [teff, logg, feh, 1, abund5193]
            else:
                abund4443 = np.float32(df[column][46])
                abund4468 = np.float32(df[column][49])
                abund4501 = np.float32(df[column][50])
                abund4571 = np.float32(df[column][53])
                if abund4443!= -666.0:


                    star_dict2[starn][4443] = [teff, logg, feh,  1,abund4443]
                if abund4468!= -666.0:

                    star_dict2[starn][4468] = [teff, logg, feh,  1,abund4468]
                if abund4501 != -666.0:

                    star_dict2[starn][4501] = [teff, logg, feh,  1,abund4501]
                if abund4571 != -666.0:

                    star_dict2[starn][4571] = [teff, logg, feh,  1,abund4571]

    import pickle

    print([star_dict[i].keys() for i in star_dict.keys()])
    print(1,star_dict)

    print(2, star_dict2)

imblist = []
import matplotlib.pyplot as plt
for starn in star_dict:
    if starn in star_dict2:
        lte1 = np.asarray([i[-1] for i in star_dict[starn].values()])
        lte2 = np.asarray([i[-1] for i in star_dict2[starn].values()])
        if len(lte1) < 2 or len(lte2) < 2:
            continue

        feh = [i[2] for i in star_dict[starn].values()]
        logg = [i[1] for i in star_dict[starn].values()]
        feh = float(feh[0])
        logg = float(logg[0])
        lte1 = lte1 - 4.9 - feh
        lte2 = lte2 - 4.9 - feh


        imb = np.mean(lte1)-np.mean(lte2)
        if np.isnan(imb):
           continue
        print(logg)
        imblist.append(imb)

        plt.scatter(feh, imb)
print(len(imblist))
plt.show()