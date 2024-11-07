# Finding stars like the sun
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl

import seaborn as sns

from scipy import stats

def solar_twin_lines(element):
    b = pkl.load(open(element + "_Galah_impacts.pkl", "rb"))
    c = pkl.load(open(element + "_Galah_impacts_lineless.pkl", "rb"))
    feh = np.asarray([i[2] for i in c.values()])
    logg = np.asarray([i[1] for i in c.values()])
    teff = np.asarray([i[0] for i in c.values()])
    # Indexes of stars close to our sun's parameters
    twindex = np.where(
        np.logical_and(
        np.logical_and(
            abs(teff-5772) <= 60,
            abs(logg-4.438) <= 0.05),
            abs(feh) <= 0.05)

    )[0]

    # we have indexes of solar twins,n ot the star number which is whbat we use in the dictioanry so now we convert
    keylist = []
    for key in c.keys():
        keylist.append(key)
    index_list = (np.asarray(keylist)[twindex])
    # lists to contain the lkte and nlte titanium abundances (lte from galah, nlte from us) for solar twins for us to
    # use line by line
    lte = {}
    nlte = {}
    mean = {}
    for starn in index_list:
        for line in b[starn]:
            if line in lte:
                lte[line].append(b[starn][line][-3])
                nlte[line].append(b[starn][line][-2])
            else:
                lte[line]=[b[starn][line][-3]]
                nlte[line] = [b[starn][line][-2]]
    for line in lte:
        print(line, np.mean(lte[line]), len(lte[line]))
        print(line, np.mean(nlte[line]), len(nlte[line]), "\n")
        if len(lte[line]) > 5:
            mean[line] = [np.mean(lte[line]), np.mean(nlte[line])]

solar_twin_lines("Ti2")