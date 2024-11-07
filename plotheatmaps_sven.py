import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl

import scipy.interpolate
import seaborn as sns

from scipy import stats
import os
import pickle
import sys
import pandas as pd
from scipy.interpolate import interpn
from scipy.interpolate import interp1d
import warnings

from astropy.io import fits

atomname = r"/home/jama2357/Documents/TiFeAtoms/DepartureGrid/GALAH_DR3_main_allspec_v2.fits"
atom = fits.open(atomname)
atom1 = atom[1].data

# Basic NLTE correction plot to show the difference between ti 1 and 2
def plot_basic_corrections(element):
    b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
    a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))

    print("Number of stars", len(a[0]), len(a[1]))
    plt.scatter(a[0], a[1], s=0.5)
    print(min(a[1]))


    #plt.title("Mean NLTE corrections to Galah for "+element+ " averaged over all lines")
    plt.xlabel("Metallicity [Fe/H]")
    if "1" in element:
        plt.ylabel("[Ti/H]$_\mathrm{Ti\ I,\ NLTE}$ - [Ti/H]$_\mathrm{Ti\ I,\ LTE}$")
    elif "2" in element:
        plt.ylabel("[Ti/H]$_\mathrm{Ti\ II,\ NLTE}$ - [Ti/H]$_{Ti II\ LTE}$")


    #plt.savefig("graphs/"+element+"NLTE_Corrections.jpg")
    plt.show()

# Basic corrections for NLTE on the same plot.
def plot_both():
    for element in ["Ti1", "Ti2"]:

        b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
        a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))

        print("Number of stars for ", element, ": ", len(a[1]))
        if element == "Ti1":
            plt.scatter(a[0], a[1], s=0.5, label = "Ti I")
        elif element == "Ti2":
            plt.scatter(a[0], a[1], s=0.5, label = "Ti II")

        print(min(a[1]))
        plt.legend()
        #plt.title("Mean NLTE corrections to Galah averaged over all lines")
        plt.xlabel("Metallicity [Fe/H]")
        if "1" in element:
            plt.ylabel("[Ti/H]$_\mathrm{Ti\ I,\ NLTE}$ - [Ti/H]$_\mathrm{Ti\ I,\ LTE}$")
        elif "2" in element:
            plt.ylabel("[Ti/H]$_\mathrm{Ti\ II,\ NLTE}$ - [Ti/H]$_{Ti II\ LTE}$")

    #plt.savefig("graphs/Both_NLTE_Corrections.jpg")
    plt.show()
import mpl_scatter_density # adds projection='scatter_density'
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Colormap

# Heat map for LTE and NLTE for Ti 2 and 1 with giants and dwarfs
def color_single(element):
    b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
    a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
    c = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
    sns.set_theme(style="whitegrid")

    hr = ((np.where(np.isnan(a[1])))[0])
    print(len(a[0]))
    print(len(b))
    fehgiants = []
    fehdwarfs = []
    nltegiants = []
    nltedwarfs = []
    dwarfcount = 0
    giantcount = 0
    index = 0

    trendnltedwarfs = []
    trendfehdwarfs = []
    trendnltegiants = []
    trendfehgiants = []


    logg = []
    # Open up the dict to find the logg. Works only because dictionary is in the same order as the list. Weak and
    # an issue as dicts arent MEANT to have order.
    for starn in list(b.keys()):

        lg = ([i[1] for i in b[starn].values()])
        lgfix = lg[0]
        logg.append(lgfix)

    # Checking every 0.2 metallicity bin so we have a good averaged trend line
    for y in np.arange(min(a[0]), max(a[0]) , 0.2):
        windwarf = np.where(np.logical_and(np.logical_and(y -0.1 <= a[0], a[0] <= y+0.1), np.asarray(logg) > 3.5))[0]
        winmeandwarf = np.mean(np.asarray(a[1])[windwarf])
        if not len(windwarf) < 30:
            trendfehdwarfs.append(y)
            trendnltedwarfs.append(winmeandwarf)

        wingiant = np.where(np.logical_and(np.logical_and(y -0.1 <= a[0], a[0] <= y+0.1), np.asarray(logg) <= 3.5))[0]
        if not len(wingiant) < 30:
            winmeangiant = np.mean(np.asarray(a[1])[wingiant])

            trendfehgiants.append(y)
            trendnltegiants.append(winmeangiant)




    lte = []
    for x in (b.keys()):
        lte.append(b[x][list(b[x].keys())[0]][4])
        g = (b[x][list(b[x].keys())[0]][1])
        feh = (b[x][list(b[x].keys())[0]][2])
        nlte = (b[x][list(b[x].keys())[0]][6])
        if g <= 3.5:
            giantcount+=1
        elif g > 3.5:
            dwarfcount += 1
    # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact
    nlte = np.delete(np.asarray(a[1]), hr)
    feh = np.delete(np.asarray(a[0]), hr)


    comparelte = False
    if comparelte:
        plt.scatter(lte, nlte, s=0.5)
        plt.xlabel("LTE preediction by GALAH")
        plt.ylabel("NLTE correction")
        plt.show()
        print("Exiting due to comparelte variable")
        exit()


    print((np.where(np.isnan(nlte))))
    print("Number of stars for ", element, ": ", len(nlte))

    #sns.kdeplot(x=feh, y=nlte, fill=True, cmap="viridis")

    ax = plt.figure().add_subplot(1,1,1, projection="scatter_density")
    den = ax.scatter_density(feh, nlte, cmap=LinearSegmentedColormap.from_list("white_viridis",  [
        (0, '#ffffff'),
        (1e-20, '#89f089'),
        (0.15, '#CDDD1A'),
        (1, '#e21616'),
    ], N=1000))
    #den = ax.scatter_density(feh, nlte, cmap=LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N).reversed())

    #den = ax.scatter_density(feh, nlte, cmap="Reds")

    plt.colorbar(den, label="Number of stars")
    if element == "Ti1":
        pass
        #plt.title("Mean NLTE corrections to Galah for Ti I averaged over all lines")
    elif element == "Ti2":
        pass
        #plt.title("Mean NLTE corrections to Galah for Ti II averaged over all lines")


    #plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(dwarfcount)+")")
    #plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(giantcount)+")")
    sns.lineplot(x=trendfehdwarfs,y= trendnltedwarfs, color="black", label = "Dwarf (N="+str(dwarfcount)+")")
    sns.lineplot(x=trendfehgiants, y=trendnltegiants, linestyle="--",  color="blue", label = "Giant (N="+str(giantcount)+")")

    plt.legend()
    plt.xlabel("Metallicity [Fe/H]")
    if "1" in element:
        plt.ylabel("[Ti/H]$_\mathrm{Ti\ I,\ NLTE}$ - [Ti/H]$_\mathrm{Ti\ I,\ LTE}$")
    elif "2" in element:
        plt.ylabel("[Ti/H]$_\mathrm{Ti\ II,\ NLTE}$ - [Ti/H]$_\mathrm{Ti\ II,\ LTE}$")
    #plt.savefig("graphs/NLTEimpact_Colormap_"+element, bbox_inches="tight")
    plt.show()

#heat map for both lte/nlte in the same one
def color_both():
    for element in ["Ti1", "Ti2"]:
        print("Not working yet, exiting")
        exit()
        b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
        a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
        hr = ((np.where(np.isnan(a[1])))[0])
        print(len(a[0]))
        print(len(b))
        fehgiants = []
        fehdwarfs = []
        nltegiants = []
        nltedwarfs = []
        dwarfcount = 0
        giantcount = 0
        index = 0
        for x in (b.keys()):

            g = (b[x][list(b[x].keys())[0]][1])
            feh = (b[x][list(b[x].keys())[0]][2])
            nlte = (b[x][list(b[x].keys())[0]][6])
            if g > 3.5:
                giantcount+=1
                fehgiants.append(feh)
                nltegiants.append(nlte)
            elif g <= 3.5:
                dwarfcount += 1
                fehdwarfs.append(feh)
                nltedwarfs.append(nlte)
        # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact
        nlte = np.delete(np.asarray(a[1]), hr)
        feh = np.delete(np.asarray(a[0]), hr)



        print((np.where(np.isnan(nlte))))
        print("Number of stars for ", element, ": ", len(nlte))
        if element == "Ti1":
            ax = plt.figure().add_subplot(1,1,1, projection="scatter_density")
        den = ax.scatter_density(feh, nlte, cmap=LinearSegmentedColormap.from_list("white_viridis",  [
            (0, '#ffffff'),
            (1e-20, '#89f089'),
            (0.15, '#CDDD1A'),
            (1, '#e21616'),
        ], N=1000))
        #den = ax.scatter_density(feh, nlte, cmap=LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N).reversed())

        #den = ax.scatter_density(feh, nlte, cmap="Reds")
        if element == "Ti2":

            plt.colorbar(den, label="Number of stars")
        z=np.polyfit(fehdwarfs, nltedwarfs, 1)
        p = np.poly1d(z)
        ax.plot(fehdwarfs, p(fehdwarfs))
    #plt.title("Mean NLTE corrections to Galah for Ti averaged over all lines")

    plt.xlabel("Metallicity [Fe/H]")
    if "1" in element:
        plt.ylabel("[Ti/H]$_\mathrm{Ti\ I,\ NLTE}$ - [Ti/H]$_\mathrm{Ti\ I,\ LTE}$")
    elif "2" in element:
        plt.ylabel("[Ti/H]$_\mathrm{Ti\ II,\ NLTE}$ - [Ti/H]$_{Ti II\ LTE}$")

    plt.show()



# Basic comparison of LTE to NLTE in one graph
def basic_ltenlte_compare(element):
    b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
    a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
    c = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
    feh = [i[2] for i in c.values()]
    lte = [i[-2] for i in c.values()]
    nlte = [i[-1]+i[-2] for i in c.values()]
    print(len(feh), len(nlte))

    plt.scatter(feh, lte, c="blue", label = "LTE abundance", s=1)
    plt.scatter(feh, nlte, c="red", label = "Non-LTE abundance", s=1)

    plt.legend()
    #plt.savefig("graphs/"+element+"_simpleLTE2.png")
    plt.show()

# Heatmap of stars now showing either LTE abundance or true NLTE abundance (not correction)
def heatmap_lte(element, nltetrue):
    b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
    a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
    c = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
    if nltetrue:
        nltestring = "NLTE"
    else:
        nltestring="LTE"
    feh = np.asarray([i[2] for i in c.values()])

    if nltetrue:
        lte = np.asarray([i[-1]+i[-2] for i in c.values()])
    else:
        lte = np.asarray([i[-2] for i in c.values()])

    logg = np.asarray([i[1] for i in c.values()])
    nltegiants = []
    nltedwarfs = []
    dwarfcount = 0
    giantcount = 0
    index = 0

    trendltedwarfs = []
    trendfehdwarfs = []
    trendltegiants = []
    trendfehgiants = []


    # Checking every 0.2 metallicity bin so we have a good averaged trend line for a mean value of nlte
    for y in np.arange(min(feh), max(feh) , 0.2):
        windwarf = np.where(np.logical_and(np.logical_and(y -0.1 <= feh, feh <= y+0.1), np.asarray(logg) > 3.5))[0]
        winmeandwarf = np.mean(np.asarray(lte)[windwarf])
        if not len(windwarf) < 30:
            trendfehdwarfs.append(y)
            trendltedwarfs.append(winmeandwarf)

        wingiant = np.where(np.logical_and(np.logical_and(y -0.1 <= feh, feh <= y+0.1), np.asarray(logg) <= 3.5))[0]
        if not len(wingiant) < 30:
            winmeangiant = np.mean(np.asarray(lte)[wingiant])

            trendfehgiants.append(y)
            trendltegiants.append(winmeangiant)




    # Now separating the feh and nlte impact for dwarfs and giants
    dwarfindex = np.where(np.asarray(logg) > 3.5)
    giantindex = np.where(np.asarray(logg) <= 3.5)






    print("Number of stars for ", element, ": ", len(lte))
    ax = plt.figure().add_subplot(1,1,1, projection="scatter_density")
    den = ax.scatter_density(feh, lte, cmap=LinearSegmentedColormap.from_list("white_viridis",  [
        (0, '#ffffff'),
        (1e-20, '#89f089'),
        (0.15, '#CDDD1A'),
        (1, '#e21616'),
    ], N=1000))
    #den = ax.scatter_density(feh, nlte, cmap=LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N).reversed())

    #den = ax.scatter_density(feh, nlte, cmap="Reds")

    plt.colorbar(den, label="Number of stars")
    if element == "Ti1":
        pass
        #plt.title("Mean "+nltestring+" stellar abundance from Galah for Ti I averaged over all lines")
    elif element == "Ti2":
        pass
        #plt.title("Mean "+nltestring+" stellar abundance from Galah for Ti II averaged over all lines")


    plt.plot(trendfehdwarfs, trendltedwarfs, c="black", label = "Dwarf (N="+str(len(dwarfindex[0]))+")")
    plt.plot(trendfehgiants, trendltegiants, "--",  c="blue", label = "Giant (N="+str(len(giantindex[0]))+")")
    plt.legend()
    plt.xlabel("Metallicity [Fe/H]")
    if "1" in element:
        plt.ylabel("[Ti/H]$_\mathrm{Ti\ I,\ " + nltestring + "}$")
    elif "2" in element:
        plt.ylabel("[Ti/H]$_\mathrm{Ti\ II,\ " + nltestring + "}$")


    #plt.savefig("graphs/"+nltestring+"Colormap_"+element+"2.png")
    plt.show()


# Finding stars like the sun
def solar_twin_lines(element):
    b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
    a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
    c = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
    feh = np.asarray([i[2] for i in c.values()])
    lte = np.asarray([i[-2] for i in c.values()])
    nlte = np.asarray([i[-1]+i[-2] for i in c.values()])
    logg = np.asarray([i[1] for i in c.values()])
    teff = np.asarray([i[0] for i in c.values()])
    # Indexes of stars close to our sun's parameters
    twindex = np.where(
        np.logical_and(
        np.logical_and(
            abs(teff-5772) <= 100,
            abs(logg-4.44) <= 0.1),
            abs(feh) <= 0.05
        )

    )[0]
    print(twindex)
    print(len(twindex))


from statistics import stdev as std

# checking line by line scatter by calculating standard deviation of each star in lte and nlte using the std
# of all lines available in both. E.g. if one line isn't available in NLTE we remove it from lte calculation too.
def scatter():
    elementlist  = ["Ti2", "Ti1"]
    for element in elementlist:
        if element == "Ti2":
            elementstr = "Ti II"
            scattercolour = "darkorange"
            scatter2 = "blue"
        else:
            elementstr = "Ti I"
            scattercolour = "darkorange"
            scatter2= "blue"
        # Lists for standard deviations
        ltestd = []
        nltestd = []
        # we don't use some stars so must make a new list to allow for the same length of x and y axis when plotting
        fehlist = []
        logglist = []
        load = True

        if load:
            print("Loading")
            fehlist = pkl.load(open(element+"STD_FeHlist", "rb"))
            logglist = pkl.load(open(element+"STD_logglist", "rb"))
            nlte_std_vals = pkl.load(open(element+"STD_NLTE", "rb"))
            lte_std_vals = pkl.load(open(element+"STD_LTE", "rb"))
            linedict = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))

            plt.scatter(fehlist, nlte_std_vals, label = elementstr + " standard deviation in NLTE", alpha = 0.7, s=2, c=scattercolour)
            plt.scatter(fehlist, lte_std_vals, label = elementstr + " standard deviation in LTE", alpha = 0.7, s=2, c=scatter2)



            for lteval in [True, False]:
                trendfehdwarfs = []
                trenddwarfs = []
                trendfehgiants = []
                trendgiants = []
                if lteval:
                    std_val = lte_std_vals
                    ltestring = "LTE"
                    linec = "red"
                else:
                    std_val = nlte_std_vals
                    ltestring = "NLTE"
                    linec="black"
                # Find the mean STD in a metallicity range (y) and plot the mean and the mertallciity (y)
                for y in np.arange(min(fehlist), max(fehlist), 0.2):
                    windwarf = np.where(np.logical_and(np.logical_and(y - 0.1 <= fehlist, fehlist <= y + 0.1), np.asarray(logglist) > 3.5))[0]

                    winmeandwarf = np.mean(np.asarray(std_val)[windwarf])
                    if not len(windwarf) < 15:
                        trendfehdwarfs.append(y)
                        trenddwarfs.append(winmeandwarf)

                    wingiant = \
                    np.where(np.logical_and(np.logical_and(y - 0.1 <= fehlist, fehlist <= y + 0.1), np.asarray(logglist) <= 3.5))[0]
                    if not len(wingiant) < 15:
                        winmeangiant = np.mean(np.asarray(std_val)[wingiant])

                        trendfehgiants.append(y)
                        trendgiants.append(winmeangiant)
                plt.plot(trendfehgiants, trendgiants, label = elementstr + " giants in "+ltestring, color=linec, linestyle = "--", linewidth=2)
                plt.plot(trendfehdwarfs, trenddwarfs, label = elementstr+" dwarfs in " + ltestring, color=linec, linewidth=2)
        else:
            linedict = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
            a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
            c = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
            feh = np.asarray([i[2] for i in c.values()])
            lte = np.asarray([i[-2] for i in c.values()])
            nlte = np.asarray([i[-1] + i[-2] for i in c.values()])
            logg = np.asarray([i[1] for i in c.values()])
            teff = np.asarray([i[0] for i in c.values()])

            for starn in linedict.keys():

                lineplot = []
                nltelineplot = []
                for line in linedict[starn]:
                    lte = linedict[starn][line][-3]
                    nlte = linedict[starn][line][-2]
                    impact = linedict[starn][line][-1]
                    if np.isnan(lte) or np.isnan(nlte):
                        continue
                    # Won't work with float32
                    lineplot.append(np.float64(lte))
                    nltelineplot.append(np.float64(nlte))
                if len(lineplot) < 4 or len(nltelineplot) < 4:
                    continue
                fehlist.append(linedict[starn][line][2])
                ltestd.append(std(lineplot))
                nltestd.append(std(nltelineplot))
                logglist.append(linedict[starn][line][1])
            #plt.scatter(fehlist,ltestd, c="blue", label = "LTE", s=1)
            #plt.scatter(fehlist, nltestd, c="red", label = "NLTE", s=1)
            pkl.dump(np.asarray(nltestd), open(element+"STD_NLTE", "wb"))
            pkl.dump(np.asarray(ltestd), open(element+"STD_LTE", "wb"))
            pkl.dump(logglist, open(element+"STD_logglist", "wb"))
            plt.scatter(fehlist, np.asarray(nltestd), label = element + " standard deviation", s=1, alpha = 0.3)
            plt.scatter(fehlist, np.asarray(ltestd), label = element + " standard deviation", s=1, alpha = 0.3)
            pkl.dump(fehlist, open(element+"STD_FeHlist", "wb"))
        plt.ylabel("Standard Deviation of [Ti/H]")
        plt.xlabel("[Fe/H]")
        plt.legend()
        #plt.hlines(0, min(fehlist), max(fehlist), colors="purple")
        #plt.savefig("graphs/"+element+"_STD2.png")
        plt.show()

# Plotion imbalance for dwarfs and giants
def ion_imbalance():
    verbose = False

    for element in ["Ti1"]:
        # Teff logg feh vmic lteabund nlteIMPACT
        impact = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        if element == "Ti1":
            alt_impact = pkl.load(open("Ti2_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        else:
            alt_impact = pkl.load(open("Ti1_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        nlte_imbalance_list = []
        lte_imbalance_list = []
        logg_list = []
        feh_list = []
        """load = False
        if load:
            lists = pkl.load(open("ion_imbalance_feh_nlte_lte2_svenlines.pkl", "rb"))
            nlte_imbalance_list = lists[1]
            feh_list = lists[0]
            lte_imbalance_list=lists[2]
        else:"""
        # Make a list of imbalance in non lte and lte, including the feh and logg as we limit to stars with both
        # ti2 and ti1 lines (which is not done previously)
        for starn in impact:
            if starn in alt_impact:
                ti1lte = impact[starn][-2]
                ti2lte = alt_impact[starn][-2]
                ti1nlteimpact = impact[starn][-1]
                ti2nlteimpact = alt_impact[starn][-1]
                ti1nlte = ti1lte + ti1nlteimpact
                ti2nlte = ti2lte + ti2nlteimpact
                feh_list.append(impact[starn][2])
                nlte_imbalance_list.append(ti2nlte-ti1nlte)
                lte_imbalance_list.append(ti2lte-ti1lte)
                logg_list.append(impact[starn][1])
                if verbose:
                    if impact[starn][2] < -2:
                        # [teff, logg, feh, xi, ltemean, nlte mmean impact]
                        print("LTE, NLTE impact, NLTE abund of Ti I:   ", impact[starn][-2],  impact[starn][-1], impact[starn][-2]+ impact[starn][-1])
                        print("LTE, NLTE impact, NLTE abund of Ti II:  ", alt_impact[starn][-2], alt_impact[starn][-1], alt_impact[starn][-2]+  alt_impact[starn][-1])
                        print("NLTE and LTE ion imbalance, NLTE effect:", ti2nlte-ti1nlte, ti2lte-ti1lte, abs(ti2nlte-ti1nlte)-abs(ti2lte-ti1lte))
                        print("\n")

        # As usual, plot the mean values in a bin range of metallicity, separating dwarfs and giants based on logg
        for lteval in [True, False]:
            trendfehdwarfs = []
            trenddwarfs = []
            trendfehgiants = []
            trendgiants = []

            if lteval:
                imbalance_list = lte_imbalance_list
                ltestring = "LTE"
                clr = "red"
            else:
                imbalance_list = nlte_imbalance_list
                ltestring = "Non-LTE"
                clr = "black"
            for y in np.arange(min(feh_list), max(feh_list), 0.3):
                windwarf = np.where(np.logical_and(np.logical_and(y - 0.15 <= feh_list, feh_list <= y + 0.15), np.asarray(logg_list) > 3.5))[0]

                winmeandwarf = np.mean(np.asarray(imbalance_list)[windwarf])
                if not len(windwarf) < 10:
                    trendfehdwarfs.append(y)
                    trenddwarfs.append(winmeandwarf)

                wingiant = \
                np.where(np.logical_and(np.logical_and(y - 0.15 <= feh_list, feh_list <= y + 0.15), np.asarray(logg_list) <= 3.5))[0]
                if not len(wingiant) < 10:
                    winmeangiant = np.nanmean(np.asarray(imbalance_list)[wingiant])

                    trendfehgiants.append(y)
                    trendgiants.append(winmeangiant)
                    if verbose:
                        print(lteval, y, winmeangiant, len(wingiant))
            plt.scatter(feh_list, imbalance_list, label=ltestring+" ionisation imbalance", s=4)
            plt.plot(trendfehgiants, trendgiants, label = ltestring+" trend line of Giant imbalance", linestyle="--", color = clr)
            plt.plot(trendfehdwarfs, trenddwarfs, label = ltestring+" trend line of Dwarf imbalance", linestyle = "-", color = clr)

        #pkl.dump([feh_list, nlte_imbalance_list, lte_imbalance_list, logg_list], open("ion_imbalance_feh_nlte_lte2_svenlines.pkl", "wb"))


        #plt.plot(trendfehgiants, trendgiants, label = "Trend line of Giant STD", color="blue")
        #plt.plot(trendfehdwarfs, trenddwarfs, label = "Trend line of Dwarf STD", color="black")

        plt.xlabel("[Fe/H]")
        plt.ylabel("[Ti/H]$_\mathrm{Ti\ II}$ - [Ti/H]$_\mathrm{Ti\ I}$")
        plt.legend()
        #plt.savefig("graphs/ionisation_imbalance_both2.png", bbox_inches="tight")
        plt.show()


def excitation_balance():
    for lteval in [False, True]:
        if lteval:
            ltestring = "LTE"
        else:
            ltestring = "NLTE"
        element = "Ti1"
        elementstr = "Ti I"
        linedict = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
        lineless = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        """for x in [329822, 390082, 165167, 64889, 302568, 285544, 595404, 155698, 500335, 211862, 542323, 588259, 488379, 233073, 614187]:
            print(lineless[x])"""

        gradientlist = []
        feh_list = []
        feh = [i[2] for i in lineless.values()]
        counter = 0
        length = len(lineless)
        for starn in linedict:
            ltelist = []
            linelist = []

            counter+=1
            if counter % 5000 == 0:
                print(str(counter)+"/"+str(length))
            """teff = linedict[starn][linedict[starn].keys()[0]][0]
            feh = linedict[starn][4758][2]
            logg = linedict[starn][4758][1]"""
            #feh = linedict[x][list(b[x].keys())[0]][2]
            if not (len(linedict[starn])) > 6:
                continue
            for line in linedict[starn]:
                if lteval:
                    lte = linedict[starn][line][4]
                else:
                    lte = linedict[starn][line][5]
                nlte_impact = linedict[starn][line][6]
                if (np.isnan(lte)):
                    continue
                ltelist.append(lte)
                linelist.append(line)
            feh_list.append(linedict[starn][line][2])

            # calc the gradient of the excitation imbalance. 0 is best
            z = np.polyfit(linelist, ltelist, 1)
            gradientlist.append(z[0])
        print("Ready to plot")
        plt.scatter(feh_list, gradientlist, s=8)
        """ax = plt.figure().add_subplot(1,1,1, projection="scatter_density")
        den = ax.scatter_density(feh, gradientlist, cmap=LinearSegmentedColormap.from_list("white_viridis",  [
            (0, '#ffffff'),
            (1e-20, '#89f089'),
            (0.15, '#CDDD1A'),
            (1, '#e21616'),
        ], N=1000))"""
        #den = ax.scatter_density(feh, nlte, cmap=LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N).reversed())

        #den = ax.scatter_density(feh, nlte, cmap="Reds")
        #plt.colorbar(den, label="Number of stars")
        #plt.title("Gradient of excitation imbalances in "+ltestring)
        #plt.savefig("graphs/Exctitation_imbalance_"+ltestring+"2.png")
        plt.xlabel("[Fe/H]")
        plt.ylabel("Gradient of excitation imbalance")
        plt.show()

# Checking actual LTE vs NLTE abundances to comapre to GCE trends
def galactic_evolution():
    for element in ["Ti1", "Ti2"]:
        for lteval in [True, False]:
            b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
            a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
            lineless = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
            sns.set_theme(style="whitegrid")
            y = 0
            hr = ((np.where(np.isnan(a[1])))[0])
            print(len(a[0]))
            print(len(b))
            fehgiants = []
            fehdwarfs = []
            nltegiants = []
            nltedwarfs = []
            dwarfcount = 0
            giantcount = 0
            index = 0

            trendnltedwarfs = []
            trendfehdwarfs = []
            trendnltegiants = []
            trendfehgiants = []


            # Checking numbero f dwarfs and giants by counting those with more or less than 3.5 logg ([1] in lineless)
            for g in [i[1] for i in lineless.values()]:
                if g <= 3.5:
                    giantcount+=1
                elif g > 3.5:
                    dwarfcount += 1
            # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact





            feh = [i[2] for i in lineless.values()]
            logg = [i[1] for i in lineless.values()]
            if lteval:
                abundance = [i[-2] for i in lineless.values()]-np.asarray(feh)
            else:
                abundance = [i[-1]+i[-2] for i in lineless.values()]-np.asarray(feh)
            #sns.kdeplot(x=feh, y=nlte, fill=True, cmap="viridis")
            # Checking every 0.2 metallicity bin so we have a good averaged trend line

            for y in np.arange(min(feh), max(feh) , 0.2):
                windwarf = np.where(np.logical_and(np.logical_and(y -0.1 <= feh, feh <= y+0.1), np.asarray(logg) > 3.5))[0]
                winmeandwarf = np.mean(np.asarray(abundance)[windwarf])
                if not len(windwarf) < 30:
                    trendfehdwarfs.append(y)
                    trendnltedwarfs.append(winmeandwarf)

                wingiant = np.where(np.logical_and(np.logical_and(y -0.1 <= feh, feh <= y+0.1), np.asarray(logg) <= 3.5))[0]
                if not len(wingiant) < 30:
                    winmeangiant = np.mean(np.asarray(abundance)[wingiant])

                    trendfehgiants.append(y)
                    trendnltegiants.append(winmeangiant)



            if element == "Ti1":
                if lteval:
                    pos = 1
                else:
                    pos = 3
            else:
                if lteval:
                    pos = 2
                else:
                    pos = 4
            ax = plt.figure().add_subplot(1,1,1, projection="scatter_density")
            den = ax.scatter_density(feh, abundance, cmap=LinearSegmentedColormap.from_list("white_viridis",  [
                (0, '#ffffff'),
                (1e-20, '#89f089'),
                (0.15, '#CDDD1A'),
                (1, '#e21616'),
            ], N=1000))
            #den = ax.scatter_density(feh, nlte, cmap=LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N).reversed())

            #den = ax.scatter_density(feh, nlte, cmap="Reds")

            plt.colorbar(den, label="Number of stars")
            if lteval:
                ltestring = "LTE"
            else:
                ltestring = "NLTE"
            if element == "Ti1":
                pass
                #plt.title("Mean "+ltestring+" abundance to Galah for Ti I averaged over all lines")
            elif element == "Ti2":
                pass
                #plt.title("Mean "+ltestring+" abundance to Galah for Ti II averaged over all lines")


            #plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(dwarfcount)+")")
            #plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(giantcount)+")")
            sns.lineplot(x=trendfehdwarfs,y= trendnltedwarfs, color="black", label = "Dwarf (N="+str(dwarfcount)+")")
            sns.lineplot(x=trendfehgiants, y=trendnltegiants, linestyle="--",  color="blue", label = "Giant (N="+str(giantcount)+")")
            if lteval:

                plt.legend()
            else:
                plt.legend([], [], frameon=False)
            plt.xlabel("Metallicity [Fe/H]")
            if "1" in element:
                plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I,\ " + ltestring + "}$")
            elif "2" in element:
                plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ II,\ " + ltestring + "}$")

            #plt.savefig("graphs/GCE_abundances_"+element+ltestring+"2.png")
            plt.show()

# The trendline to match up for GCE. NLTE and LTE in Ti 1 and 2.
def galactic_evolution_trendlines():
    trenddict = {}
    for element in ["Ti1", "Ti2"]:
        trenddict[element] = {}
        for lteval in [False, True]:
            b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
            a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
            lineless = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
            sns.set_theme(style="whitegrid")

            hr = ((np.where(np.isnan(a[1])))[0])
            print(len(a[0]))
            print(len(b))
            fehgiants = []
            fehdwarfs = []
            nltegiants = []
            nltedwarfs = []
            dwarfcount = 0
            giantcount = 0
            index = 0

            trendnltedwarfs = []
            trendfehdwarfs = []
            trendnltegiants = []
            trendfehgiants = []

            # Checking numbero f dwarfs and giants by counting those with more or less than 3.5 logg ([1] in lineless)
            for g in [i[1] for i in lineless.values()]:
                if g <= 3.5:
                    giantcount += 1
                elif g > 3.5:
                    dwarfcount += 1
            # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact

            feh = [i[2] for i in lineless.values()]
            logg = [i[1] for i in lineless.values()]
            if lteval:
                abundance = [i[-2] for i in lineless.values()]
            else:
                abundance = [i[-1] + i[-2] for i in lineless.values()]
            # sns.kdeplot(x=feh, y=nlte, fill=True, cmap="viridis")
            # Checking every 0.2 metallicity bin so we have a good averaged trend line

            for y in np.arange(min(feh), max(feh), 0.2):
                windwarf = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) > 3.5))[0]
                winmeandwarf = np.mean(np.asarray(abundance)[windwarf])
                if not len(windwarf) < 30:
                    trendfehdwarfs.append(y)
                    trendnltedwarfs.append(winmeandwarf)

                wingiant = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) <= 3.5))[0]
                if not len(wingiant) < 30:
                    winmeangiant = np.mean(np.asarray(abundance)[wingiant])

                    trendfehgiants.append(y)
                    trendnltegiants.append(winmeangiant)


            # plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(dwarfcount)+")")
            # plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(giantcount)+")")

            trenddict[element][lteval] = [trendnltedwarfs, trendnltegiants, trendfehdwarfs, trendfehgiants]





    for element in ["Ti1", "Ti2"]:
        for lteval in [True, False]:

            if lteval:
                ltestring = "LTE"
                if element == "Ti1":
                    elementstring = "Ti I"
                    col = "blue"

                else:
                    col = "red"
                    elementstring = "Ti II"

            else:
                ltestring = "NLTE"
                if element == "Ti1":
                    col = "black"
                    elementstring = "Ti I"

                else:
                    col = "darkorange"
                    elementstring = "Ti II"
            print(lteval, trenddict[element][lteval][2], trenddict[element][lteval][0], "\n", trenddict[element][lteval][3], trenddict[element][lteval][1])
            sns.lineplot(x=trenddict[element][lteval][2], y=trenddict[element][lteval][0], label=elementstring+" "+ ltestring+" dwarfs", color = col)
            sns.lineplot(x=trenddict[element][lteval][3], y=trenddict[element][lteval][1], linestyle="--", color=col,
                         label=elementstring+" " + ltestring+" giants")
        #plt.title("Titanium abundance trend lines")
        plt.legend()
        plt.xlabel("Metallicity [Fe/H]")
        if "1" in element:
            plt.ylabel("[Ti/H]$_\mathrm{Ti\ I}$")
        elif "2" in element:
            plt.ylabel("[Ti/H]$_\mathrm{Ti\ II}$")

        #plt.savefig("graphs/"+element+"GCE_Trendlines2.png", bbox_inches="tight")
        plt.show()

def hr_diag():
    for element in ["Ti1"]:
        # Teff logg feh vmic lteabund nlteIMPACT
        impact = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        if element == "Ti1":
            alt_impact = pkl.load(open("Ti2_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        else:
            alt_impact = pkl.load(open("Ti1_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        nlteimpact = []
        lte_imbalance_list = []
        logg_list = []
        feh_list = []
        temp_list = []
        """load = False
        if load:
            lists = pkl.load(open("ion_imbalance_feh_nlte_lte2_svenlines.pkl", "rb"))
            nlte_imbalance_list = lists[1]
            feh_list = lists[0]
            lte_imbalance_list=lists[2]
        else:"""
        # Make a list of imbalance in non lte and lte, including the feh and logg as we limit to stars with both
        # ti2 and ti1 lines (which is not done previously)
        for starn in impact:
            ti1lte = impact[starn][-2]
            ti1nlteimpact = impact[starn][-1]
            feh_list.append(impact[starn][2])
            nlteimpact.append(ti1nlteimpact)
            logg_list.append(impact[starn][1])
            temp_list.append(impact[starn][0])
        plt.scatter(temp_list, logg_list, c=nlteimpact, s= 2, cmap=LinearSegmentedColormap.from_list("white_viridis", [
            (0, '#ffffff'),
            (1e-20, '#ffffff'),
            (0.25, '#1f32db'),
            (0.5, '#740699'),
            (1, '#e21616'),
        ], N=1000))
        """ax = plt.figure().add_subplot(1, 1, 1, projection="scatter_density")
        den = ax.scatter_density(temp_list, logg_list, cmap=LinearSegmentedColormap.from_list("white_viridis", [
            (0, '#ffffff'),
            (1e-20, '#89f089'),
            (0.15, '#CDDD1A'),
            (1, '#e21616'),
        ], N=1000))
        """
        plt.colorbar( label="$\Delta$[Ti/Fe]$_\mathrm{NLTE}$")
        plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()
        # den = ax.scatter_density(feh, nlte, cmap=LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N).reversed())

        # den = ax.scatter_density(feh, nlte, cmap="Reds")


        #pkl.dump([feh_list, nlte_imbalance_list, lte_imbalance_list, logg_list], open("ion_imbalance_feh_nlte_lte2_svenlines.pkl", "wb"))


        #plt.plot(trendfehgiants, trendgiants, label = "Trend line of Giant STD", color="blue")
        #plt.plot(trendfehdwarfs, trenddwarfs, label = "Trend line of Dwarf STD", color="black")

        plt.xlabel("T$_{eff}$")
        plt.ylabel("Log($g$)")
        #plt.savefig("graphs/HRdiagNLTE2.png")
        plt.show()

from operator import itemgetter
def exc_energy(element):
    b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
    a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
    sns.set_theme(style="whitegrid")
    hr = ((np.where(np.isnan(a[1])))[0])
    fehgiants = []
    fehdwarfs = []
    nltegiants = []
    nltedwarfs = []
    dwarfcount = 0
    giantcount = 0
    index = 0

    if element == "Ti1":
        lines = ["4758", "4759", "4781", "4797", "4820", "5689", "5716", "5720", "5739", "5866"
                 ]
        exc_energy_list = [2.2492, 2.2556,  0.8484,  2.3341, 1.5025, 2.2968, 2.2968, 2.2916, 2.2492,  1.0666]
    elif element == "Ti2":
        lines = ["4874", "4865", "4849", "4798", "4764", "4719"]
        exc_energy_list = [3.0948, 1.1156, 1.1305, 1.0800, 1.2369, 1.2429]
    print((list(b)))


    for metal in [-2, -1.5, -1, -0.5, 0, 0.5]:
        nlte_means = []

        metal_array = np.asarray(a[0])
        metal_indexes = np.where(np.logical_and(metal_array <= metal, metal_array > metal-0.5))
        star_numbers = np.asarray(list(b))[metal_indexes]

        for line_index in range(len(lines)):
            line = lines[line_index]
            energy = exc_energy_list[line_index]

            nlte_values = []
            for starn in star_numbers:
                if int(line) in b[starn]:
                    nlte_values.append(b[starn][int(line)][-1])
                else:
                    continue
                    nlte_values.append( np.nan)
            nlte_mean = np.nanmean(nlte_values)
            nlte_means.append(nlte_mean)
            print(metal, energy, nlte_mean)

        plot_nlte = [x for _, x in sorted(zip(exc_energy_list, nlte_means))]
        plot_energy = sorted(exc_energy_list)
        plt.plot(plot_energy, plot_nlte, label="[Fe/H] <= "+str(metal))

    plt.ylabel("$\Delta$[Ti/Fe]$_{\mathrm{Ti\ I,\ NLTE}}$")
    plt.xlabel("Excitation energy (eV)")
    plt.legend()
    #plt.savefig("graphs/"+str(element)+"excitation_energy_NLTE2.png")
    plt.show()





def color_single_fe(element):
    b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
    a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
    c = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
    sns.set_theme(style="whitegrid")

    hr = ((np.where(np.isnan(a[1])))[0])
    print(len(a[0]))
    print(len(b))
    fehgiants = []
    fehdwarfs = []
    nltegiants = []
    nltedwarfs = []
    dwarfcount = 0
    giantcount = 0
    index = 0

    trendnltedwarfs = []
    trendfehdwarfs = []
    trendnltegiants = []
    trendfehgiants = []


    logg = []
    # Open up the dict to find the logg. Works only because dictionary is in the same order as the list. Weak and
    # an issue as dicts arent MEANT to have order.
    for starn in list(b.keys()):

        lg = ([i[1] for i in b[starn].values()])
        lgfix = lg[0]
        logg.append(lgfix)

    # Checking every 0.2 metallicity bin so we have a good averaged trend line
    for y in np.arange(min(a[0]), max(a[0]) , 0.2):
        windwarf = np.where(np.logical_and(np.logical_and(y -0.1 <= a[0], a[0] <= y+0.1), np.asarray(logg) > 3.5))[0]
        winmeandwarf = np.mean(np.asarray(a[1])[windwarf])
        if not len(windwarf) < 30:
            trendfehdwarfs.append(y)
            trendnltedwarfs.append(winmeandwarf)

        wingiant = np.where(np.logical_and(np.logical_and(y -0.1 <= a[0], a[0] <= y+0.1), np.asarray(logg) <= 3.5))[0]
        if not len(wingiant) < 30:
            winmeangiant = np.mean(np.asarray(a[1])[wingiant])

            trendfehgiants.append(y)
            trendnltegiants.append(winmeangiant)




    lte = []
    for x in (b.keys()):
        lte.append(b[x][list(b[x].keys())[0]][4])
        g = (b[x][list(b[x].keys())[0]][1])
        feh = (b[x][list(b[x].keys())[0]][2])
        nlte = (b[x][list(b[x].keys())[0]][6])
        if g <= 3.5:
            giantcount+=1
        elif g > 3.5:
            dwarfcount += 1
    # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact
    nlte = np.delete(np.asarray(a[1]), hr)
    feh = np.delete(np.asarray(a[0]), hr)


    comparelte = False
    if comparelte:
        plt.scatter(lte, nlte, s=0.5)
        plt.xlabel("LTE preediction by GALAH")
        plt.ylabel("NLTE correction")
        plt.show()
        print("Exiting due to comparelte variable")
        exit()


    print((np.where(np.isnan(nlte))))
    print("Number of stars for ", element, ": ", len(nlte))

    #sns.kdeplot(x=feh, y=nlte, fill=True, cmap="viridis")

    ax = plt.figure().add_subplot(1,1,1, projection="scatter_density")
    den = ax.scatter_density(feh, nlte, cmap=LinearSegmentedColormap.from_list("white_viridis",  [
        (0, '#ffffff'),
        (1e-20, '#89f089'),
        (0.15, '#CDDD1A'),
        (1, '#e21616'),
    ], N=1000))
    #den = ax.scatter_density(feh, nlte, cmap=LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N).reversed())

    #den = ax.scatter_density(feh, nlte, cmap="Reds")

    plt.colorbar(den, label="Number of stars")
    if element == "Ti1":
        pass
        #plt.title("Mean NLTE corrections to Galah for Ti I averaged over all lines")
    elif element == "Ti2":
        pass
        #plt.title("Mean NLTE corrections to Galah for Ti II averaged over all lines")


    #plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(dwarfcount)+")")
    #plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(giantcount)+")")
    sns.lineplot(x=trendfehdwarfs,y= trendnltedwarfs, color="black", label = "Dwarf (N="+str(dwarfcount)+")")
    sns.lineplot(x=trendfehgiants, y=trendnltegiants, linestyle="--",  color="blue", label = "Giant (N="+str(giantcount)+")")

    plt.legend()
    plt.xlabel("Metallicity [Fe/H]")
    if "1" in element:

        plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I,\ NLTE}$ - [Ti/Fe]$_\mathrm{Ti\ I,\ LTE}$")
        plt.ylabel("$\Delta\mathrm{[Ti/Fe]_{Ti\ I,\ NLTE}}$")
    elif "2" in element:
        plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ II,\ NLTE}$ - [Ti/Fe]$_\mathrm{Ti\ II,\ LTE}$")
        plt.ylabel("$\Delta\mathrm{[Ti/Fe]_{Ti\ II,\ NLTE}}$")
    plt.ylim(-0.05, 0.25)
    plt.xlim(-2.5, 0.5)
    #plt.savefig("graphs/NLTEimpact_Colormap_Fe_"+element+"_2.png", bbox_inches="tight")
    plt.show()
def scatter_fe():
    elementlist  = ["Ti2", "Ti1"]
    for element in elementlist:
        if element == "Ti2":
            elementstr = "Ti II"
            scattercolour = "darkorange"
            scatter2 = "blue"
        else:
            elementstr = "Ti I"
            scattercolour = "darkorange"
            scatter2= "blue"
        # Lists for standard deviations
        ltestd = []
        nltestd = []
        # we don't use some stars so must make a new list to allow for the same length of x and y axis when plotting
        fehlist = []
        logglist = []
        load = True

        if load:
            print("Loading")
            fehlist = pkl.load(open(element+"STD_FeHlist", "rb"))
            logglist = pkl.load(open(element+"STD_logglist", "rb"))
            nlte_std_vals = pkl.load(open(element+"STD_NLTE", "rb"))
            lte_std_vals = pkl.load(open(element+"STD_LTE", "rb"))
            linedict = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))

            sns.scatterplot(x=fehlist, y=nlte_std_vals, label = elementstr + " standard deviation in NLTE",  s=8)
            sns.scatterplot(x=fehlist, y=lte_std_vals, label = elementstr + " standard deviation in LTE",  s=8)



            for lteval in [True, False]:
                trendfehdwarfs = []
                trenddwarfs = []
                trendfehgiants = []
                trendgiants = []
                if lteval:
                    std_val = lte_std_vals
                    ltestring = "LTE"
                    linec = "red"
                else:
                    std_val = nlte_std_vals
                    ltestring = "NLTE"
                    linec="black"
                # Find the mean STD in a metallicity range (y) and plot the mean and the mertallciity (y)
                for y in np.arange(min(fehlist), max(fehlist), 0.2):
                    windwarf = np.where(np.logical_and(np.logical_and(y - 0.1 <= fehlist, fehlist <= y + 0.1), np.asarray(logglist) > 3.5))[0]

                    winmeandwarf = np.mean(np.asarray(std_val)[windwarf])
                    if not len(windwarf) < 15:
                        trendfehdwarfs.append(y)
                        trenddwarfs.append(winmeandwarf)

                    wingiant = \
                    np.where(np.logical_and(np.logical_and(y - 0.1 <= fehlist, fehlist <= y + 0.1), np.asarray(logglist) <= 3.5))[0]
                    if not len(wingiant) < 15:
                        winmeangiant = np.mean(np.asarray(std_val)[wingiant])

                        trendfehgiants.append(y)
                        trendgiants.append(winmeangiant)
                plt.plot(trendfehgiants, trendgiants, label = elementstr + " giants in "+ltestring, color=linec, linestyle = "--", linewidth=2)
                plt.plot(trendfehdwarfs, trenddwarfs, label = elementstr+" dwarfs in " + ltestring, color=linec, linewidth=2)
        else:
            linedict = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
            a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
            c = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
            feh = np.asarray([i[2] for i in c.values()])
            lte = np.asarray([i[-2] for i in c.values()])
            nlte = np.asarray([i[-1] + i[-2] for i in c.values()])
            logg = np.asarray([i[1] for i in c.values()])
            teff = np.asarray([i[0] for i in c.values()])

            for starn in linedict.keys():

                lineplot = []
                nltelineplot = []
                for line in linedict[starn]:
                    lte = linedict[starn][line][-3]-linedict[starn][line][2]
                    nlte = linedict[starn][line][-2]-linedict[starn][line][2]
                    impact = linedict[starn][line][-1]
                    if np.isnan(lte) or np.isnan(nlte):
                        continue
                    # Won't work with float32
                    lineplot.append(np.float64(lte))
                    nltelineplot.append(np.float64(nlte))
                if len(lineplot) < 4 or len(nltelineplot) < 4:
                    continue
                fehlist.append(linedict[starn][line][2])
                ltestd.append(std(lineplot))
                nltestd.append(std(nltelineplot))
                logglist.append(linedict[starn][line][1])
            #plt.scatter(fehlist,ltestd, c="blue", label = "LTE", s=1)
            #plt.scatter(fehlist, nltestd, c="red", label = "NLTE", s=1)
            pkl.dump(np.asarray(nltestd), open(element+"STD_NLTE", "wb"))
            pkl.dump(np.asarray(ltestd), open(element+"STD_LTE", "wb"))
            pkl.dump(logglist, open(element+"STD_logglist", "wb"))
            plt.scatter(fehlist, np.asarray(nltestd), label = element + " standard deviation", s=1, alpha = 0.3)
            plt.scatter(fehlist, np.asarray(ltestd), label = element + " standard deviation", s=1, alpha = 0.3)
            pkl.dump(fehlist, open(element+"STD_FeHlist", "wb"))
        plt.ylabel("Standard Deviation of [Ti/Fe]")
        plt.xlabel("[Fe/H]")
        plt.legend()
        #plt.hlines(0, min(fehlist), max(fehlist), colors="purple")
        #plt.savefig("graphs/"+element+"_STD_fe2.png")
        plt.show()
import pandas as pd
def ion_imbalance_fe():
    verbose = True
    heatmap = True
    if heatmap:
        ax = plt.figure().add_subplot(1, 1, 1, projection="scatter_density")

    for element in ["Ti1"]:
        # Teff logg feh vmic lteabund nlteIMPACT
        impact = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        # just roll with it, it's intentional
        if element == "Ti1":
            alt_impact = pkl.load(open("Ti2_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        else:
            alt_impact = pkl.load(open("Ti1_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        nlte_imbalance_list = []
        lte_imbalance_list = []
        logg_list = []
        feh_list = []
        """load = False
        if load:
            lists = pkl.load(open("ion_imbalance_feh_nlte_lte2_svenlines.pkl", "rb"))
            nlte_imbalance_list = lists[1]
            feh_list = lists[0]
            lte_imbalance_list=lists[2]
        else:"""
        # Make a list of imbalance in non lte and lte, including the feh and logg as we limit to stars with both
        # ti2 and ti1 lines (which is not done previously)
        for starn in impact:
            if starn in alt_impact:
                # Removing Fe/H of the star from the LTE value to get Ti/Fe as it's saved as Ti/H. No imapct elsewhere
                # as the nlte different is the same for both Ti/H and /Fe for maths reasons
                ti1lte = impact[starn][-2]-impact[starn][2]
                ti2lte = alt_impact[starn][-2]-impact[starn][2]
                ti1nlteimpact = impact[starn][-1]
                ti2nlteimpact = alt_impact[starn][-1]
                ti1nlte = ti1lte + ti1nlteimpact
                ti2nlte = ti2lte + ti2nlteimpact
                feh_list.append(impact[starn][2])
                nlte_imbalance_list.append(ti1nlte-ti2nlte)
                lte_imbalance_list.append(ti1lte-ti2lte)
                logg_list.append(impact[starn][1])
                if verbose:
                    if impact[starn][2] < -2:
                        # [teff, logg, feh, xi, ltemean, nlte mmean impact]
                        print("LTE, NLTE impact, NLTE abund of Ti I:   ", impact[starn][-2],  impact[starn][-1], impact[starn][-2]+ impact[starn][-1])
                        print("LTE, NLTE impact, NLTE abund of Ti II:  ", alt_impact[starn][-2], alt_impact[starn][-1], alt_impact[starn][-2]+  alt_impact[starn][-1])
                        print("NLTE and LTE ion imbalance, NLTE effect:", ti2nlte-ti1nlte, ti2lte-ti1lte, abs(ti2nlte-ti1nlte)-abs(ti2lte-ti1lte))
                        print("\n")

        # As usual, plot the mean values in a bin range of metallicity, separating dwarfs and giants based on logg
        for lteval in [True, False]:
            trendfehdwarfs = []
            trenddwarfs = []
            trendfehgiants = []
            trendgiants = []

            if lteval:
                imbalance_list = lte_imbalance_list
                ltestring = "LTE"
                clr = "red"
                scatclr = "tab:orange"
            else:
                imbalance_list = nlte_imbalance_list
                ltestring = "Non-LTE"
                clr = "black"
                scatclr = "blue"
            for y in np.arange(min(feh_list), max(feh_list), 0.2):
                windwarf = np.where(np.logical_and(np.logical_and(y - 0.1 <= feh_list, feh_list <= y + 0.1), np.asarray(logg_list) > 3.5))[0]

                winmeandwarf = np.mean(np.asarray(imbalance_list)[windwarf])
                if not len(windwarf) < 30:
                    trendfehdwarfs.append(y)
                    trenddwarfs.append(winmeandwarf)
                    if verbose:
                        print("dwarf, lteval?, metallicity, mean imbalance, number of stars", lteval, y, winmeandwarf, len(windwarf))

                wingiant = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh_list, feh_list <= y + 0.1), np.asarray(logg_list) <= 3.5))[0]
                if not len(wingiant) < 30:
                    winmeangiant = np.nanmean(np.asarray(imbalance_list)[wingiant])

                    trendfehgiants.append(y)
                    trendgiants.append(winmeangiant)
                    if verbose:
                        print("giant, lteval?, metallicity, mean imbalance, number of stars", lteval, y, winmeangiant, len(wingiant))

            if not heatmap:
                sns.scatterplot(x=feh_list,y= imbalance_list, label=ltestring+" ionisation imbalance", s=8)
            #plt.plot(trendfehgiants, trendgiants, label = ltestring+" trend line of Giant imbalance", linestyle="--", color = clr)
            #plt.plot(trendfehdwarfs, trenddwarfs, label = ltestring+" trend line of Dwarf imbalance", linestyle = "-", color = clr)
            sns.lineplot(x=trendfehgiants, y=trendgiants, linestyle="--", color = clr)
            sns.lineplot(x=trendfehdwarfs, y=trenddwarfs, linestyle = "-", color = clr)

        #pkl.dump([feh_list, nlte_imbalance_list, lte_imbalance_list, logg_list], open("ion_imbalance_feh_nlte_lte2_svenlines.pkl", "wb"))


        #plt.plot(trendfehgiants, trendgiants, label = "Trend line of Giant STD", color="blue")
        #plt.plot(trendfehdwarfs, trenddwarfs, label = "Trend line of Dwarf STD", color="black")
        if heatmap:
            den = ax.scatter_density(feh_list, imbalance_list, cmap=LinearSegmentedColormap.from_list("white_viridis", [
                (0, '#ffffff'),
                (1e-20, '#89f089'),
                (0.15, '#CDDD1A'),
                (1, '#e21616'),
            ], N=1000))

            plt.colorbar(den, label="Number of stars")

        plt.xlabel("[Fe/H]")
        plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I}$ - [Ti/Fe]$_\mathrm{Ti\ II}$")
        plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I-Ti\ II}$")
        #plt.legend()
        #plt.savefig("graphs/ionisation_imbalance_both_Fe2.png", bbox_inches="tight")
        plt.show()

# trendlines for GCE but for ti1 and ti2, ignoring giants and dwarf differentials
def galactic_evolution_trendlines_ionmean():
    trenddict = {}
    for element in ["Ti1", "Ti2"]:
        trenddict[element] = {}
        for lteval in [False, True]:
            b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
            a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
            lineless = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
            sns.set_theme(style="whitegrid")

            hr = ((np.where(np.isnan(a[1])))[0])
            print(len(a[0]))
            print(len(b))
            fehgiants = []
            fehdwarfs = []
            nltegiants = []
            nltedwarfs = []
            dwarfcount = 0
            giantcount = 0
            index = 0

            trendnltedwarfs = []
            trendfehdwarfs = []
            trendnltegiants = []
            trendfehgiants = []

            # Checking numbero f dwarfs and giants by counting those with more or less than 3.5 logg ([1] in lineless)
            for g in [i[1] for i in lineless.values()]:
                if g <= 3.5:
                    giantcount += 1
                elif g > 3.5:
                    dwarfcount += 1
            # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact

            feh = [i[2] for i in lineless.values()]
            logg = [i[1] for i in lineless.values()]

            if lteval:
                abundance = np.asarray([i[-2] for i in lineless.values()])-np.asarray(feh)
            else:
                abundance = np.asarray([i[-1] + i[-2] for i in lineless.values()])-np.asarray(feh)
            # sns.kdeplot(x=feh, y=nlte, fill=True, cmap="viridis")
            # Checking every 0.2 metallicity bin so we have a good averaged trend line

            for y in np.arange(min(feh), max(feh), 0.1):
                windwarf = \
                np.where(np.logical_and(y - 0.05 <= feh, feh <= y + 0.05))[0]
                winmeandwarf = np.mean(np.asarray(abundance)[windwarf])
                if not len(windwarf) < 30:
                    trendfehdwarfs.append(y)
                    trendnltedwarfs.append(winmeandwarf)

            # plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(dwarfcount)+")")
            # plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(giantcount)+")")

            trenddict[element][lteval] = [trendnltedwarfs, trendfehdwarfs]





    for element in ["Ti1", "Ti2"]:
        for lteval in [True, False]:

            if lteval:
                ltestring = "LTE"
                if element == "Ti1":
                    elementstring = "Ti I"
                    col = "blue"

                else:
                    col = "red"
                    elementstring = "Ti II"

            else:
                ltestring = "NLTE"
                if element == "Ti1":
                    col = "black"
                    elementstring = "Ti I"

                else:
                    col = "darkorange"
                    elementstring = "Ti II"
            sns.lineplot(x=trenddict[element][lteval][1], y=trenddict[element][lteval][0], label=elementstring+" "+ ltestring, color = col)
        #plt.title("Titanium abundance trend lines")
        plt.legend()
        plt.xlabel("Metallicity [Fe/H]")
        plt.ylabel("[Ti/Fe]")
    print(len(lineless))
    #plt.savefig("graphs/GCE_Trendlines_overall2.png", bbox_inches="tight")
    plt.show()

sns.set_theme(style="whitegrid")
sns.color_palette("dark")
def galactic_evolution_trendlines_overallmean():
    trenddict = {}
    for lteval in [False, True]:
        fehgiants = []
        fehdwarfs = []
        nltegiants = []
        nltedwarfs = []
        dwarfcount = 0
        giantcount = 0
        index = 0

        trendnltedwarfs = []
        trendfehdwarfs = []
        trendnltegiants = []
        trendfehgiants = []

        for element in ["Ti1", "Ti2"]:

            b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
            a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
            lineless = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
            sns.set_theme(style="whitegrid")

            hr = ((np.where(np.isnan(a[1])))[0])
            print(len(a[0]))
            print(len(b))

            # Checking numbero f dwarfs and giants by counting those with more or less than 3.5 logg ([1] in lineless)
            for g in [i[1] for i in lineless.values()]:
                if g <= 3.5:
                    giantcount += 1
                elif g > 3.5:
                    dwarfcount += 1
            # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact

            feh = [i[2] for i in lineless.values()]
            logg = [i[1] for i in lineless.values()]

            if lteval:
                abundance = np.asarray([i[-2] for i in lineless.values()])-np.asarray(feh)
            else:
                abundance = np.asarray([i[-1] + i[-2] for i in lineless.values()])-np.asarray(feh)
            # sns.kdeplot(x=feh, y=nlte, fill=True, cmap="viridis")
            # Checking every 0.2 metallicity bin so we have a good averaged trend line

            for y in np.arange(min(feh), max(feh), 0.2):
                windwarf = \
                np.where(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1))[0]
                winmeandwarf = np.mean(np.asarray(abundance)[windwarf])
                if not len(windwarf) < 20:
                    trendfehdwarfs.append(y)
                    trendnltedwarfs.append(winmeandwarf)

            # plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(dwarfcount)+")")
            # plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(giantcount)+")")
            print(lteval)
            trenddict[lteval] = [trendnltedwarfs, trendfehdwarfs]





    for lteval in [True, False]:

        if lteval:
            ltestring = "LTE"
            col = "blue"

            col = "red"

        else:
            ltestring = "NLTE"
            col = "black"

        sns.lineplot(x=trenddict[lteval][1], y=trenddict[lteval][0], label=ltestring, color = col)
    #plt.title("Titanium abundance trend lines")
    plt.legend()
    plt.xlabel("Metallicity [Fe/H]")
    plt.ylabel("[Ti/Fe]")

    ##plt.savefig("graphs/GCE_Trendlines_overall2.png", bbox_inches="tight")
    plt.show()
def galactic_evolution_trendlines_fe():
    trenddict = {}
    for element in ["Ti1", "Ti2"]:
        ax = plt.figure().add_subplot(1, 1, 1, projection="scatter_density")

        trenddict[element] = {}
        for lteval in [False, True]:
            b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
            a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
            lineless = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
            sns.set_theme(style="whitegrid")

            hr = ((np.where(np.isnan(a[1])))[0])
            print(len(a[0]))
            print(len(b))
            fehgiants = []
            fehdwarfs = []
            nltegiants = []
            nltedwarfs = []
            dwarfcount = 0
            giantcount = 0
            index = 0

            trendnltedwarfs = []
            trendfehdwarfs = []
            trendnltegiants = []
            trendfehgiants = []

            # Checking numbero f dwarfs and giants by counting those with more or less than 3.5 logg ([1] in lineless)
            for g in [i[1] for i in lineless.values()]:
                if g <= 3.5:
                    giantcount += 1
                elif g > 3.5:
                    dwarfcount += 1
            # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact

            feh = [i[2] for i in lineless.values()]
            logg = [i[1] for i in lineless.values()]
            # turn from /H to /Fe by minusing metallicity fe/h
            if lteval:
                abundance = np.asarray([i[-2] for i in lineless.values()])-np.asarray(feh)
                if element == "Ti1":
                    abundanceuse = np.asarray([i[-2] for i in lineless.values()])-np.asarray(feh)
                    fehuse = feh
            else:
                abundance = np.asarray([i[-1] + i[-2] for i in lineless.values()])-np.asarray(feh)
            # sns.kdeplot(x=feh, y=nlte, fill=True, cmap="viridis")
            # Checking every 0.2 metallicity bin so we have a good averaged trend line

            for y in np.arange(min(feh), max(feh), 0.2):
                windwarf = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) > 3.5))[0]
                winmeandwarf = np.mean(np.asarray(abundance)[windwarf])
                if not len(windwarf) < 30:
                    trendfehdwarfs.append(y)
                    trendnltedwarfs.append(winmeandwarf)

                wingiant = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) <= 3.5))[0]
                if not len(wingiant) < 30:
                    winmeangiant = np.mean(np.asarray(abundance)[wingiant])

                    trendfehgiants.append(y)
                    trendnltegiants.append(winmeangiant)


            # plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(dwarfcount)+")")
            # plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(giantcount)+")")

            trenddict[element][lteval] = [trendnltedwarfs, trendnltegiants, trendfehdwarfs, trendfehgiants]





    for element in ["Ti1", "Ti2"]:
        for lteval in [True, False]:

            if lteval:
                ltestring = "LTE"
                if element == "Ti1":
                    elementstring = "Ti I"
                    col = "blue"

                else:
                    col = "red"
                    elementstring = "Ti II"

            else:
                ltestring = "NLTE"
                if element == "Ti1":
                    col = "black"
                    elementstring = "Ti I"

                else:
                    col = "darkorange"
                    elementstring = "Ti II"
            print(lteval, trenddict[element][lteval][2], trenddict[element][lteval][0], "\n", trenddict[element][lteval][3], trenddict[element][lteval][1])
            print(trenddict[element][lteval][0][-3] - trenddict[element][lteval][1][-3], "\n")

            sns.lineplot(x=trenddict[element][lteval][2], y=trenddict[element][lteval][0], label=elementstring+" "+ ltestring+" dwarfs", color = col)
            sns.lineplot(x=trenddict[element][lteval][3], y=trenddict[element][lteval][1], linestyle="--", color=col,
                         label=elementstring+" " + ltestring+" giants")
        den = ax.scatter_density(fehuse, abundanceuse, cmap=LinearSegmentedColormap.from_list("white_viridis", [
            (0, '#ffffff'),
            (1e-20, '#89f089'),
            (0.15, '#CDDD1A'),
            (1, '#e21616'),
        ], N=1000))

        plt.colorbar(den, label="Number of stars")

        #plt.title("Titanium abundance trend lines")
        plt.legend()
        plt.xlabel("Metallicity [Fe/H]")
        if "1" in element:
            plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I}$")
        elif "2" in element:
            plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ II}$")

        #plt.savefig("graphs/"+element+"GCE_Trendlines_fe2.png", bbox_inches="tight")
        plt.show()

def galactic_evolution_trendlines_fe_line(wl):
    trenddict = {}
    if wl in [4874, 4865, 4849, 4798, 4764, 4719]:
        element = "Ti2"
    else:
        element = "Ti1"
    print("WL", wl)
    trenddict[element] = {}
    abundancelte = []
    abundancenlte = []

    for lteval in [False, True]:
        lined = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
        a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
        lineless = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        sns.set_theme(style="whitegrid")

        hr = ((np.where(np.isnan(a[1])))[0])
        fehgiants = []
        fehdwarfs = []
        nltegiants = []
        nltedwarfs = []
        dwarfcount = 0
        giantcount = 0
        index = 0

        trendnltedwarfs = []
        trendfehdwarfs = []
        trendnltegiants = []
        trendfehgiants = []

        # Checking numbero f dwarfs and giants by counting those with more or less than 3.5 logg ([1] in lineless)
        for g in [i[1] for i in lineless.values()]:
            if g <= 3.5:
                giantcount += 1
            elif g > 3.5:
                dwarfcount += 1
        # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact

        feh = [i[2] for i in lineless.values()]
        logg = [i[1] for i in lineless.values()]
        feh = []
        logg = []
        abundance = []
        # turn from /H to /Fe by minusing metallicity fe/h
        for i in lined.values():
            try:
                if lteval:
                    abundance.append(i[wl][-3]-i[wl][2])
                    abundancelte.append(i[wl][-3]-i[wl][2])
                    feh.append(i[wl][2])
                    logg.append(i[wl][1])

                else:
                    abundance.append(i[wl][-2]-i[wl][2])
                    abundancenlte.append(i[wl][-2]-i[wl][2])
                    feh.append(i[wl][2])
                    logg.append(i[wl][1])

            except KeyError:
                pass
        print(feh)
        for y in np.arange(min(feh), max(feh), 0.2):
            windwarf = \
            np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) > 3.5))[0]
            winmeandwarf = np.mean(np.asarray(abundance)[windwarf])
            if not len(windwarf) < 30:
                trendfehdwarfs.append(y)
                trendnltedwarfs.append(winmeandwarf)

            wingiant = \
            np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) <= 3.5))[0]
            if not len(wingiant) < 30:
                winmeangiant = np.mean(np.asarray(abundance)[wingiant])

                trendfehgiants.append(y)
                trendnltegiants.append(winmeangiant)


        # plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(dwarfcount)+")")
        # plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(giantcount)+")")

        trenddict[element][lteval] = [trendnltedwarfs, trendnltegiants, trendfehdwarfs, trendfehgiants]





    for lteval in [True, False]:

        if lteval:
            abundance_use = abundancelte
            ltestring = "LTE"
            if element == "Ti1":
                elementstring = "Ti I"
                col = "blue"

            else:
                col = "blue"
                elementstring = "Ti II"

        else:
            abundance_use = abundancenlte
            ltestring = "NLTE"
            if element == "Ti1":
                col = "black"
                elementstring = "Ti I"

            else:
                col = "black"
                elementstring = "Ti II"
        print(lteval, trenddict[element][lteval][2], trenddict[element][lteval][0], "\n", trenddict[element][lteval][3], trenddict[element][lteval][1])
        print(trenddict[element][lteval][0][-3] - trenddict[element][lteval][1][-3], "\n")
        sns.lineplot(x=trenddict[element][lteval][2], y=trenddict[element][lteval][0], label=elementstring+" "+ ltestring+" dwarfs", color = col)
        sns.lineplot(x=trenddict[element][lteval][3], y=trenddict[element][lteval][1], linestyle="--", color=col,
                     label=elementstring+" " + ltestring+" giants")
        sns.scatterplot(x=feh, y =abundance_use )
    #plt.title("Titanium abundance trend lines")
    plt.legend()
    plt.xlabel("Metallicity [Fe/H]")
    if "1" in element:
        plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I}$")
    elif "2" in element:
        plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ II}$")

    #plt.savefig("graphs/wlscatter/"+element+"GCE_Trendlines_fe_"+str(wl)+"2.png", bbox_inches="tight")
    plt.show(block=False)
    plt.close()
def ion_imbalance_fe_line():
    count4778 = 0
    countothers = 0
    counterval = 0
    verbose = True
    heatmap = True
    if heatmap:
        ax = plt.figure().add_subplot(1, 1, 1, projection="scatter_density", sharey='row')

    for element in ["Ti1"]:

        # Teff logg feh vmic lteabund nlteIMPACT
        impact_orig = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
        impact = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
        # just roll with it, it's intentional
        if element == "Ti1":
            alt_impact = pkl.load(open("Ti2_Galah_impacts2_svenlines.pkl", "rb"))


        nlte_imbalance_list = []
        lte_imbalance_list = []
        galah_imbalance_list = []
        logg_list = []
        feh_list = []
        ti1meanlist = []
        ti2meanlist = []
        """load = False
        if load:
            lists = pkl.load(open("ion_imbalance_feh_nlte_lte2_svenlines.pkl", "rb"))
            nlte_imbalance_list = lists[1]
            feh_list = lists[0]
            lte_imbalance_list=lists[2]
        else:"""
        # Make a list of imbalance in non lte and lte, including the feh and logg as we limit to stars with both
        # ti2 and ti1 lines (which is not done previously)
        counter = 0
        for starn in impact:
            if counter % 500 == 0:
                print(counter,"/", len(impact))
            counter+=1
            ti1count = 0
            ti2count = 0
            if starn in alt_impact:
                # Removing Fe/H of the star from the LTE value to get Ti/Fe as it's saved as Ti/H. No imapct elsewhere
                # as the nlte different is the same for both Ti/H and /Fe for maths reasons

                ti2ltelist = []
                ti2nltelist = []
                ti1ltelist = []
                ti1nltelist = []
                fehhere = ([i[2] for i in alt_impact[starn].values()][0])

                if fehhere >= -2:
                    pass

                for wl in alt_impact[starn]:
                    # applied in apply nlte galah
                    """try:
                        flag = atom1[starn]['ind_flag_Ti' + str(wl)]
                    except KeyError:
                        flag = atom1[starn]['ind_flag_Ti' + str(wl+1)]
                    if flag == 1:
                        continue"""
                    if wl == 4849:
                        pass
                    ti2count += 1

                    # LTE value by taking LTE and - feh to get ti/fe instead of ti/h
                    ti2ltelist.append(alt_impact[starn][wl][-3] - alt_impact[starn][wl][2])

                    ti2nltelist.append(alt_impact[starn][wl][-2] - alt_impact[starn][wl][2])



                for wl in impact[starn]:

                    if wl == 4778:
                        count4778+=1
                        continue
                    # we apply it in apply)_nlte_galah now
                    """try:
                        flag = atom1[starn]['ind_flag_Ti' + str(wl)]
                    except KeyError:
                        flag = atom1[starn]['ind_flag_Ti' + str(wl+1)]
                    if flag == 1:
                        continue"""

                    ti1count+=1

                    ti1ltelist.append(impact[starn][wl][-3] - impact[starn][wl][2])
                    countothers+=1

                    ti1nltelist.append(impact[starn][wl][-2] - impact[starn][wl][2])
                if len(ti1ltelist) == 0 or len(ti1nltelist) == 0 or ti1count <=2 or ti2count <=2:
                    continue
                ti2lte = np.nanmean(ti2ltelist)
                ti2nlte = np.nanmean(ti2nltelist)

                ti1lte = np.nanmean(ti1ltelist)
                ti1nlte = np.nanmean(ti1nltelist)

                feh_list.append(impact[starn][wl][2])
                nlte_imbalance_list.append(ti1nlte-ti2nlte)
                lte_imbalance_list.append(ti1lte-ti2lte)
                logg_list.append(impact[starn][wl][1])
                #galah_imbalance_list.append(atom1[starn]['Ti_fe']-atom1[starn]['Ti2_fe'])
                ti1meanlist.append(ti1count)
                ti2meanlist.append(ti2count)
                counterval += 1

                if verbose:
                    if impact[starn][wl][2] < -2:
                        # [teff, logg, feh, xi, ltemean, nlte mmean impact]
                        print("NLTE and LTE ion imbalance, NLTE effect:", ti2nlte-ti1nlte, ti2lte-ti1lte, abs(ti2nlte-ti1nlte)-abs(ti2lte-ti1lte))
                        print("\n")

        # As usual, plot the mean values in a bin range of metallicity, separating dwarfs and giants based on logg
        for lteval in [True, False]:
            trendfehdwarfs = []
            trenddwarfs = []
            trendfehgiants = []
            trendgiants = []
            trendgalah = []
            if lteval:
                imbalance_list = lte_imbalance_list
                ltestring = "LTE"
                clr = "red"
                scatclr = "tab:orange"
            else:
                imbalance_list = nlte_imbalance_list
                ltestring = "Non-LTE"
                clr = "black"
                scatclr = "blue"
            for y in np.arange(min(feh_list), max(feh_list), 0.2):
                windwarf = np.where(np.logical_and(np.logical_and(y - 0.1 <= feh_list, feh_list <= y + 0.1), np.asarray(logg_list) > 3.5))[0]

                winmeandwarf = np.mean(np.asarray(imbalance_list)[windwarf])
                if not len(windwarf) < 30:
                    trendfehdwarfs.append(y)
                    trenddwarfs.append(winmeandwarf)
                    if verbose:
                        print("dwarf, lteval?, metallicity, mean imbalance, number of stars", lteval, y, winmeandwarf, len(windwarf))

                wingiant = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh_list, feh_list <= y + 0.1), np.asarray(logg_list) <= 3.5))[0]
                if not len(wingiant) < 30:
                    winmeangiant = np.nanmean(np.asarray(imbalance_list)[wingiant])
                    #meangalah = np.nanmean(np.asarray(galah_imbalance_list)[wingiant])
                    trendfehgiants.append(y)
                    trendgiants.append(winmeangiant)
                    #trendgalah.append(meangalah)
                    if verbose:
                        print("giant, lteval?, metallicity, mean imbalance, number of stars", lteval, y, winmeangiant, len(wingiant))

            if not heatmap:
                sns.scatterplot(x=feh_list,y= imbalance_list, label=ltestring+" ionisation imbalance", s=8)
            #plt.plot(trendfehgiants, trendgiants, label = ltestring+" trend line of Giant imbalance", linestyle="--", color = clr)
            #plt.plot(trendfehdwarfs, trenddwarfs, label = ltestring+" trend line of Dwarf imbalance", linestyle = "-", color = clr)
            sns.lineplot(x=trendfehgiants, y=trendgiants, linestyle="--", color = clr, label=ltestring + " giants")
            sns.lineplot(x=trendfehdwarfs, y=trenddwarfs, linestyle = "-", color = clr, label = ltestring + " dwarfs")
            #sns.lineplot(x=trendfehgiants, y=trendgalah, linestyle="--", color = "blue", label=ltestring + " GALAH")

        #pkl.dump([feh_list, nlte_imbalance_list, lte_imbalance_list, logg_list], open("ion_imbalance_feh_nlte_lte2_svenlines.pkl", "wb"))


        #plt.plot(trendfehgiants, trendgiants, label = "Trend line of Giant STD", color="blue")
        #plt.plot(trendfehdwarfs, trenddwarfs, label = "Trend line of Dwarf STD", color="black")
        if heatmap:
            den = ax.scatter_density(feh_list, lte_imbalance_list, cmap=LinearSegmentedColormap.from_list("white_viridis", [
                (0, '#ffffff'),
                (1e-20, '#89f089'),
                (0.15, '#CDDD1A'),
                (1, '#e21616'),
            ], N=1000))

            plt.colorbar(den, label="Number of stars")
            ax.set_ylim([-0.6, 0.6])

        plt.xlabel("[Fe/H]")
        plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I}$ - [Ti/Fe]$_\mathrm{Ti\ II}$")
        plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I-Ti\ II}$")
        plt.legend()
        #plt.savefig("graphs/ionimb_lines/ionisation_imbalance_both_Fe_linespecific2.png", bbox_inches="tight")
        print(count4778, countothers, counterval)
        print(np.mean(ti1meanlist), np.mean(ti2meanlist))
        plt.show()



def galactic_evolution_trendlines_fe_line_removal():
    trenddict = {}

    for element in ["Ti2", "Ti1"]:
        ax = plt.figure().add_subplot(1, 1, 1, projection="scatter_density")
        print("Element:", element)
        abundancelte = []
        abundancenlte = []

        trenddict[element] = {}

        for lteval in [False, True]:
            print("LTE?", lteval)
            lined = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
            a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
            lineless = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
            sns.set_theme(style="whitegrid")


            lineless_original = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
            lineless = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))

            counter = 0
            print(len(lineless))
            for star in lineless_original:
                counter+=1
                if counter % 10000 == 0:
                    print(counter, "/", len(lineless_original))
                if atom1[star]['e_logg'] > 0.4 or np.isnan(atom1[star]['e_logg']):
                    lineless.pop(star)
                    lined.pop(star)

                elif atom1[star]['e_teff'] > 120 or np.isnan(atom1[star]['e_teff']):
                    lineless.pop(star)
                    lined.pop(star)

                elif atom1[star]['e_fe_h'] > 0.35 or np.isnan(atom1[star]['e_fe_h']):
                    lineless.pop(star)
                    lined.pop(star)

            hr = ((np.where(np.isnan(a[1])))[0])
            fehgiants = []
            fehdwarfs = []
            nltegiants = []
            nltedwarfs = []
            dwarfcount = 0
            giantcount = 0
            index = 0

            trendnltedwarfs = []
            trendfehdwarfs = []
            trendnltegiants = []
            trendfehgiants = []

            # Checking numbero f dwarfs and giants by counting those with more or less than 3.5 logg ([1] in lineless)
            for g in [i[1] for i in lineless.values()]:
                if g <= 3.5:
                    giantcount += 1
                elif g > 3.5:
                    dwarfcount += 1
            # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact

            feh = []
            logg = []
            abundance = []
            # turn from /H to /Fe by minusing metallicity fe/h
            starcount = 0
            for starn in lined.keys():
                i = lined[starn]
                tilinecount = 0
                if starcount % 500 == 0:
                    print(starcount, "/", len(lined))
                getnltemean = []
                getltemean = []
                for wl in i.keys():
                    # for ti1
                    if wl in [4778, 5866]:
                        continue
                    # for ti2
                    #        lines = ["4874", "4798", "4764"]

                    if wl in [4865, 4849, 4719]:
                        continue
                    try:
                        flag = atom1[starn]['ind_flag_Ti' + str(wl)]
                    except KeyError:
                        flag = atom1[starn]['ind_flag_Ti' + str(wl+1)]
                    if flag == 1:
                        continue

                    lte = i[wl][-3] - i[wl][2]
                    nlte = i[wl][-2]-i[wl][2]
                    getltemean.append(lte)
                    getnltemean.append(nlte)
                    tilinecount+=1

                if len(getnltemean) == 0 or len(getltemean) == 0 or tilinecount<2:
                    continue

                ltemean_lines = np.nanmean(getltemean)
                nltemean_lines = np.nanmean(getnltemean)



                try:
                    if lteval:
                        abundance.append(ltemean_lines)
                        abundancelte.append(ltemean_lines)
                        feh.append(i[wl][2])
                        logg.append(i[wl][1])

                    else:
                        abundance.append(nltemean_lines)
                        abundancenlte.append(nltemean_lines)
                        feh.append(i[wl][2])
                        logg.append(i[wl][1])

                except KeyError:
                    pass
                starcount+=1

            print("Star count for", element, ": ", starcount)
            for y in np.arange(min(feh), max(feh), 0.2):
                windwarf = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) > 3.5))[0]
                winmeandwarf = np.mean(np.asarray(abundance)[windwarf])
                if not len(windwarf) < 30:
                    trendfehdwarfs.append(y)
                    trendnltedwarfs.append(winmeandwarf)

                wingiant = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) <= 3.5))[0]
                if not len(wingiant) < 30:
                    winmeangiant = np.mean(np.asarray(abundance)[wingiant])

                    trendfehgiants.append(y)
                    trendnltegiants.append(winmeangiant)


            # plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(dwarfcount)+")")
            # plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(giantcount)+")")

            trenddict[element][lteval] = [trendnltedwarfs, trendnltegiants, trendfehdwarfs, trendfehgiants]





        for lteval in [True, False]:

            if lteval:
                abundance_use = abundancelte
                ltestring = "LTE"
                if element == "Ti1":
                    elementstring = "Ti I"
                    col = "blue"

                else:
                    col = "blue"
                    elementstring = "Ti II"

            else:
                abundance_use = abundancenlte
                ltestring = "NLTE"
                if element == "Ti1":
                    col = "black"
                    elementstring = "Ti I"

                else:
                    col = "black"
                    elementstring = "Ti II"
            sns.lineplot(x=trenddict[element][lteval][2], y=trenddict[element][lteval][0], label=elementstring+" "+ ltestring+" dwarfs", color = col)
            sns.lineplot(x=trenddict[element][lteval][3], y=trenddict[element][lteval][1], linestyle="--", color=col,
                         label=elementstring+" " + ltestring+" giants")

        den = ax.scatter_density(feh, abundancelte, cmap=LinearSegmentedColormap.from_list("white_viridis", [
            (0, '#ffffff'),
            (1e-20, '#89f089'),
            (0.15, '#CDDD1A'),
            (1, '#e21616'),
        ], N=1000))

        plt.colorbar(den, label="Number of stars")

        #sns.scatterplot(x=feh, y =abundance_use )

        #plt.title("Titanium abundance trend lines")
        plt.legend()
        plt.xlabel("Metallicity [Fe/H]")
        if "1" in element:
            plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I}$")
        elif "2" in element:
            plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ II}$")

        #plt.savefig("graphs/wlscatter/"+element+"GCE_Trendlines_fe_removedlines2.png", bbox_inches="tight")
        plt.show()



def galactic_evolution_trendlines_fe_splitup():
    trenddict = {}
    fig = plt.figure()
    fig.set_figwidth(15)
    fig.set_figheight(10)
    for element in ["Ti1", "Ti2"]:

        trenddict[element] = {}
        for lteval in [False, True]:
            if element == "Ti1":
                if lteval:
                    ax = fig.add_subplot(2, 2, 1, projection="scatter_density")
                    plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I}$")

                else:
                    ax = fig.add_subplot(2, 2, 2, projection="scatter_density")

            else:
                if lteval:
                    ax = fig.add_subplot(2, 2, 3, projection="scatter_density")
                    plt.xlabel("Metallicity [Fe/H]")
                    plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ II}$")

                else:
                    ax = fig.add_subplot(2, 2, 4, projection="scatter_density")
                    plt.xlabel("Metallicity [Fe/H]")

            b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
            a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
            lineless = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
            sns.set_theme(style="whitegrid")

            hr = ((np.where(np.isnan(a[1])))[0])
            print(len(a[0]))
            print(len(b))
            fehgiants = []
            fehdwarfs = []
            nltegiants = []
            nltedwarfs = []
            dwarfcount = 0
            giantcount = 0
            index = 0

            trendnltedwarfs = []
            trendfehdwarfs = []
            trendnltegiants = []
            trendfehgiants = []

            # Checking numbero f dwarfs and giants by counting those with more or less than 3.5 logg ([1] in lineless)
            for g in [i[1] for i in lineless.values()]:
                if g <= 3.5:
                    giantcount += 1
                elif g > 3.5:
                    dwarfcount += 1
            # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact

            feh = [i[2] for i in lineless.values()]
            logg = [i[1] for i in lineless.values()]
            # turn from /H to /Fe by minusing metallicity fe/h
            if lteval:
                abundance = np.asarray([i[-2] for i in lineless.values()])-np.asarray(feh)
                if element == "Ti1":
                    abundanceuse = np.asarray([i[-2] for i in lineless.values()])-np.asarray(feh)
                    fehuse = feh
            else:
                abundance = np.asarray([i[-1] + i[-2] for i in lineless.values()])-np.asarray(feh)
            # sns.kdeplot(x=feh, y=nlte, fill=True, cmap="viridis")
            # Checking every 0.2 metallicity bin so we have a good averaged trend line

            for y in np.arange(min(feh), max(feh), 0.2):
                windwarf = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) > 3.5))[0]
                winmeandwarf = np.mean(np.asarray(abundance)[windwarf])
                if not len(windwarf) < 30:
                    trendfehdwarfs.append(y)
                    trendnltedwarfs.append(winmeandwarf)

                wingiant = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) <= 3.5))[0]
                if not len(wingiant) < 30:
                    winmeangiant = np.mean(np.asarray(abundance)[wingiant])

                    trendfehgiants.append(y)
                    trendnltegiants.append(winmeangiant)


            # plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(dwarfcount)+")")
            # plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(giantcount)+")")

            trenddict[element][lteval] = [trendnltedwarfs, trendnltegiants, trendfehdwarfs, trendfehgiants]

            if lteval:
                ltestring = "LTE"
                if element == "Ti1":
                    elementstring = "Ti I"
                    col = "blue"

                else:
                    col = "blue"
                    elementstring = "Ti II"

            else:
                ltestring = "NLTE"
                if element == "Ti1":
                    col = "black"
                    elementstring = "Ti I"

                else:
                    col = "black"
                    elementstring = "Ti II"
            print(lteval, trenddict[element][lteval][2], trenddict[element][lteval][0], "\n", trenddict[element][lteval][3], trenddict[element][lteval][1])
            print(trenddict[element][lteval][0][-3] - trenddict[element][lteval][1][-3], "\n")

            sns.lineplot(x=trenddict[element][lteval][2], y=trenddict[element][lteval][0], label=elementstring+" "+ ltestring+" dwarfs", color = col)
            sns.lineplot(x=trenddict[element][lteval][3], y=trenddict[element][lteval][1], linestyle="--", color=col,
                         label=elementstring+" " + ltestring+" giants")
            print(feh, abundance)
            den = ax.scatter_density(feh, abundance, cmap=LinearSegmentedColormap.from_list("white_viridis", [
                                                 (0, '#ffffff'),
                                                 (1e-20, '#89f089'),
                                                 (0.15, '#CDDD1A'),
                                                 (1, '#e21616'),
                                             ], N=1000))

            ax.set_ylim([-0.6, 0.6])

            #plt.title("Titanium abundance trend lines")
            plt.legend(loc="lower left")
    plt.colorbar(den, label="Number of stars")

    plt.savefig("graphs/Split_GCE_Trendlines_fe2_sven.png", bbox_inches="tight")
    plt.show()



def ion_imbalance_fe_split():
    verbose = True
    heatmap = True
    fig = plt.figure()
    fig.set_figwidth(15)
    fig.set_figheight(10)

    counter = 0

    for element in ["Ti1"]:

        # Teff logg feh vmic lteabund nlteIMPACT
        impact_orig = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        impact = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        # just roll with it, it's intentional
        if element == "Ti1":
            alt_impact = pkl.load(open("Ti2_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        else:
            alt_impact = pkl.load(open("Ti1_Galah_impacts_lineless2_svenlines.pkl", "rb"))
        print(impact_orig)
        print(len(impact))
        for star in impact_orig:
            counter+=1
            if counter % 10000 == 0:
                print(counter, "/", len(impact_orig))
            if atom1[star]['e_logg'] > 0.4 or np.isnan(atom1[star]['e_logg']):
                impact.pop(star)
            elif atom1[star]['e_teff'] > 100 or np.isnan(atom1[star]['e_teff']):
                impact.pop(star)
            elif atom1[star]['e_fe_h'] > 0.1 or np.isnan(atom1[star]['e_fe_h']):
                impact.pop(star)

        print(len(impact))

        nlte_imbalance_list = []
        lte_imbalance_list = []
        logg_list = []
        feh_list = []
        """load = False
        if load:
            lists = pkl.load(open("ion_imbalance_feh_nlte_lte2_svenlines.pkl", "rb"))
            nlte_imbalance_list = lists[1]
            feh_list = lists[0]
            lte_imbalance_list=lists[2]
        else:"""
        # Make a list of imbalance in non lte and lte, including the feh and logg as we limit to stars with both
        # ti2 and ti1 lines (which is not done previously)
        for starn in impact:
            if starn in alt_impact:
                # Removing Fe/H of the star from the LTE value to get Ti/Fe as it's saved as Ti/H. No imapct elsewhere
                # as the nlte different is the same for both Ti/H and /Fe for maths reasons
                ti1lte = impact[starn][-2]-impact[starn][2]
                ti2lte = alt_impact[starn][-2]-impact[starn][2]
                ti1nlteimpact = impact[starn][-1]
                ti2nlteimpact = alt_impact[starn][-1]
                ti1nlte = ti1lte + ti1nlteimpact
                ti2nlte = ti2lte + ti2nlteimpact
                feh_list.append(impact[starn][2])
                nlte_imbalance_list.append(ti1nlte-ti2nlte)
                lte_imbalance_list.append(ti1lte-ti2lte)
                logg_list.append(impact[starn][1])
                if verbose:
                    if impact[starn][2] < -2:
                        # [teff, logg, feh, xi, ltemean, nlte mmean impact]
                        print("LTE, NLTE impact, NLTE abund of Ti I:   ", impact[starn][-2],  impact[starn][-1], impact[starn][-2]+ impact[starn][-1])
                        print("LTE, NLTE impact, NLTE abund of Ti II:  ", alt_impact[starn][-2], alt_impact[starn][-1], alt_impact[starn][-2]+  alt_impact[starn][-1])
                        print("NLTE and LTE ion imbalance, NLTE effect:", ti2nlte-ti1nlte, ti2lte-ti1lte, abs(ti2nlte-ti1nlte)-abs(ti2lte-ti1lte))
                        print("\n")

        # As usual, plot the mean values in a bin range of metallicity, separating dwarfs and giants based on logg
        for lteval in [True, False]:

            trendfehdwarfs = []
            trenddwarfs = []
            trendfehgiants = []
            trendgiants = []

            if lteval:
                imbalance_list = lte_imbalance_list
                ltestring = "LTE"
                clr = "red"
                scatclr = "tab:orange"
            else:
                imbalance_list = nlte_imbalance_list
                ltestring = "Non-LTE"
                clr = "black"
                scatclr = "blue"

            dwarffeh = []
            dwarfion = []
            giantfeh = []
            giantion = []
            dwarfs = np.where(np.asarray(logg_list)>3.5)
            giants = np.where(np.asarray(logg_list)<=3.5)
            dwarffeh=(np.asarray(feh_list)[dwarfs])
            dwarfion=(np.asarray(imbalance_list)[dwarfs])
            giantfeh=(np.asarray(feh_list)[giants])

            giantion=(np.asarray(imbalance_list)[giants])

            for y in np.arange(min(feh_list), max(feh_list), 0.2):
                windwarf = np.where(np.logical_and(np.logical_and(y - 0.1 <= feh_list, feh_list <= y + 0.1), np.asarray(logg_list) > 3.5))[0]

                winmeandwarf = np.mean(np.asarray(imbalance_list)[windwarf])
                if not len(windwarf) < 30:
                    trendfehdwarfs.append(y)
                    trenddwarfs.append(winmeandwarf)
                    if verbose:
                        print("dwarf, lteval?, metallicity, mean imbalance, number of stars", lteval, y, winmeandwarf, len(windwarf))

                wingiant = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh_list, feh_list <= y + 0.1), np.asarray(logg_list) <= 3.5))[0]

                if not len(wingiant) < 30:
                    winmeangiant = np.nanmean(np.asarray(imbalance_list)[wingiant])
                    trendfehgiants.append(y)
                    trendgiants.append(winmeangiant)
                    if verbose:
                        print("giant, lteval?, metallicity, mean imbalance, number of stars", lteval, y, winmeangiant, len(wingiant))
            print(dwarfion, "\n", dwarffeh)
            if not heatmap:
                sns.scatterplot(x=feh_list,y= imbalance_list, label=ltestring+" ionisation imbalance", s=8)
            #plt.plot(trendfehgiants, trendgiants, label = ltestring+" trend line of Giant imbalance", linestyle="--", color = clr)
            #plt.plot(trendfehdwarfs, trenddwarfs, label = ltestring+" trend line of Dwarf imbalance", linestyle = "-", color = clr)

        #pkl.dump([feh_list, nlte_imbalance_list, lte_imbalance_list, logg_list], open("ion_imbalance_feh_nlte_lte2_svenlines.pkl", "wb"))


                #plt.plot(trendfehgiants, trendgiants, label = "Trend line of Giant STD", color="blue")
                #plt.plot(trendfehdwarfs, trenddwarfs, label = "Trend line of Dwarf STD", color="black")
            for star in ['giant', 'dwarf']:
                if star == "giant":

                    if lteval:

                        ax = fig.add_subplot(2, 2, 1, projection="scatter_density")
                        plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I-Ti\ II}$")

                    else:

                        ax = fig.add_subplot(2, 2, 2, projection="scatter_density")

                    den = ax.scatter_density(giantfeh, giantion,
                                             cmap=LinearSegmentedColormap.from_list("white_viridis", [
                                                 (0, '#ffffff'),
                                                 (1e-20, '#89f089'),
                                                 (0.15, '#CDDD1A'),
                                                 (1, '#e21616'),
                                             ], N=1000))

                    plt.colorbar(den, label="Number of stars")
                    sns.lineplot(x=trendfehgiants, y=trendgiants, linestyle="--", color=clr, label = ltestring+" giants ("+str(len(giantfeh))+")")

                else:
                    if lteval:
                        ax = fig.add_subplot(2, 2, 3, projection="scatter_density")
                        plt.xlabel("Metallicity [Fe/H]")
                        plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I-Ti\ II}$")

                    else:
                        ax = fig.add_subplot(2, 2, 4, projection="scatter_density")

                        plt.xlabel("Metallicity [Fe/H]")
                    den = ax.scatter_density(dwarffeh, dwarfion,
                                             cmap=LinearSegmentedColormap.from_list("white_viridis", [
                                                 (0, '#ffffff'),
                                                 (1e-20, '#89f089'),
                                                 (0.15, '#CDDD1A'),
                                                 (1, '#e21616'),
                                             ], N=1000))
                    sns.lineplot(x=trendfehdwarfs, y=trenddwarfs, linestyle="-", color=clr, ax=ax, label = ltestring+ " dwarfs ("+ str(len(dwarffeh))+")")

                    plt.colorbar(den, label="Number of stars")
                ax.set_ylim([-0.6, 0.8])
                ax.set_xlim([-2.5, 0.5])

        #plt.legend()
        plt.savefig("graphs/Comb_ionimb2_sven.png", bbox_inches="tight")
        plt.show()


def galactic_evolution_trendlines_fe_splitup_testlogg():
    trenddict = {}
    fig = plt.figure()
    fig.set_figwidth(15)
    fig.set_figheight(10)
    for element in ["Ti2", "Ti1"]:

        trenddict[element] = {}
        for lteval in [False, True]:
            if element == "Ti1":
                if lteval:
                    ax = fig.add_subplot(2, 2, 1, projection="scatter_density")
                    plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I}$")

                else:
                    ax = fig.add_subplot(2, 2, 2, projection="scatter_density")

            else:
                if lteval:
                    ax = fig.add_subplot(2, 2, 3, projection="scatter_density")
                    plt.xlabel("Metallicity [Fe/H]")
                    plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ II}$")

                else:
                    ax = fig.add_subplot(2, 2, 4, projection="scatter_density")
                    plt.xlabel("Metallicity [Fe/H]")

            b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
            a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
            lineless_original = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
            lineless = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))

            counter = 0
            print(len(lineless))
            for star in lineless_original:
                counter+=1
                if counter % 10000 == 0:
                    print(counter, "/", len(lineless_original))
                if atom1[star]['e_logg'] > 0.4 or np.isnan(atom1[star]['e_logg']):
                    lineless.pop(star)
                elif atom1[star]['e_teff'] > 100 or np.isnan(atom1[star]['e_teff']):
                    lineless.pop(star)
                elif atom1[star]['e_fe_h'] > 0.1 or np.isnan(atom1[star]['e_fe_h']):
                    lineless.pop(star)

            print(len(lineless))

            sns.set_theme(style="whitegrid")


            dwarfcount = 0
            giantcount = 0
            index = 0

            trendnltedwarfs = []
            trendfehdwarfs = []
            trendnltegiants = []
            trendfehgiants = []

            # Checking numbero f dwarfs and giants by counting those with more or less than 3.5 logg ([1] in lineless)
            for g in [i[1] for i in lineless.values()]:
                if g <= 3.5:
                    giantcount += 1
                elif g > 3.5:
                    dwarfcount += 1
            # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact

            feh = [i[2] for i in lineless.values()]
            logg = [i[1] for i in lineless.values()]
            # turn from /H to /Fe by minusing metallicity fe/h
            if lteval:
                abundance = np.asarray([i[-2] for i in lineless.values()])-np.asarray(feh)
                if element == "Ti1":
                    abundanceuse = np.asarray([i[-2] for i in lineless.values()])-np.asarray(feh)
                    fehuse = feh
            else:
                abundance = np.asarray([i[-1] + i[-2] for i in lineless.values()])-np.asarray(feh)
            # sns.kdeplot(x=feh, y=nlte, fill=True, cmap="viridis")
            # Checking every 0.2 metallicity bin so we have a good averaged trend line

            for y in np.arange(min(feh), max(feh), 0.2):
                windwarf = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) > 3.5))[0]
                winmeandwarf = np.mean(np.asarray(abundance)[windwarf])
                if not len(windwarf) < 30:
                    trendfehdwarfs.append(y)
                    trendnltedwarfs.append(winmeandwarf)

                wingiant = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) <= 3.5))[0]
                if not len(wingiant) < 30:
                    winmeangiant = np.mean(np.asarray(abundance)[wingiant])

                    trendfehgiants.append(y)
                    trendnltegiants.append(winmeangiant)


            # plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(dwarfcount)+")")
            # plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(giantcount)+")")

            trenddict[element][lteval] = [trendnltedwarfs, trendnltegiants, trendfehdwarfs, trendfehgiants]

            if lteval:
                ltestring = "LTE"
                if element == "Ti1":
                    elementstring = "Ti I"
                    col = "blue"

                else:
                    col = "blue"
                    elementstring = "Ti II"

            else:
                ltestring = "NLTE"
                if element == "Ti1":
                    col = "black"
                    elementstring = "Ti I"

                else:
                    col = "black"
                    elementstring = "Ti II"

            sns.lineplot(x=trenddict[element][lteval][2], y=trenddict[element][lteval][0], label=elementstring+" "+ ltestring+" dwarfs", color = col)
            sns.lineplot(x=trenddict[element][lteval][3], y=trenddict[element][lteval][1], linestyle="--", color=col,
                         label=elementstring+" " + ltestring+" giants")
            den = ax.scatter_density(feh, abundance, cmap=LinearSegmentedColormap.from_list("white_viridis", [
                                                 (0, '#ffffff'),
                                                 (1e-20, '#89f089'),
                                                 (0.15, '#CDDD1A'),
                                                 (1, '#e21616'),
                                             ], N=1000))

            plt.colorbar(den, label="Number of stars")
            ax.set_ylim([-0.6, 0.6])
            ax.set_xlim([-3, 0.5])

            #plt.title("Titanium abundance trend lines")
            plt.legend(loc="lower left")
    #plt.savefig("graphs/logg_Split_GCE_Trendlines_fe2.png", bbox_inches="tight")
    plt.show()


def galactic_evolution_trendlines_fe_splitup_separatedwarfs():
    for startype in ['dwarf', 'giant']:
        trenddict = {}

        fig = plt.figure()
        fig.set_figwidth(15)
        fig.set_figheight(10)

        for element in ["Ti1", "Ti2"]:

            trenddict[element] = {}
            for lteval in [False, True]:
                if element == "Ti1":
                    if lteval:
                        ax = fig.add_subplot(2, 2, 1, projection="scatter_density")
                        plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I}$")

                    else:
                        ax = fig.add_subplot(2, 2, 2, projection="scatter_density")

                else:
                    if lteval:
                        ax = fig.add_subplot(2, 2, 3, projection="scatter_density")
                        plt.xlabel("Metallicity [Fe/H]")
                        plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ II}$")

                    else:
                        ax = fig.add_subplot(2, 2, 4, projection="scatter_density")
                        plt.xlabel("Metallicity [Fe/H]")

                b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
                a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
                lineless = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
                sns.set_theme(style="whitegrid")

                hr = ((np.where(np.isnan(a[1])))[0])
                print(len(a[0]))
                print(len(b))
                fehgiants = []
                fehdwarfs = []
                nltegiants = []
                nltedwarfs = []
                dwarfcount = 0
                giantcount = 0
                index = 0

                trendnltedwarfs = []
                trendfehdwarfs = []
                trendnltegiants = []
                trendfehgiants = []

                # Checking numbero f dwarfs and giants by counting those with more or less than 3.5 logg ([1] in lineless)
                for g in [i[1] for i in lineless.values()]:
                    if g <= 3.5:
                        giantcount += 1
                    elif g > 3.5:
                        dwarfcount += 1
                # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact

                feh = [i[2] for i in lineless.values()]
                logg = [i[1] for i in lineless.values()]
                if startype == "dwarf":
                    dwarfs = np.where(np.asarray(logg) > 3.5)
                else:
                    dwarfs = np.where(np.asarray(logg) <= 3.5)

                # turn from /H to /Fe by minusing metallicity fe/h
                if lteval:
                    abundance = np.asarray([i[-2] for i in lineless.values()])-np.asarray(feh)
                    if element == "Ti1":
                        abundanceuse = np.asarray([i[-2] for i in lineless.values()])-np.asarray(feh)
                        fehuse = feh
                else:
                    abundance = np.asarray([i[-1] + i[-2] for i in lineless.values()])-np.asarray(feh)
                # sns.kdeplot(x=feh, y=nlte, fill=True, cmap="viridis")
                # Checking every 0.2 metallicity bin so we have a good averaged trend line

                for y in np.arange(min(feh), max(feh), 0.2):
                    windwarf = \
                    np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) > 3.5))[0]
                    winmeandwarf = np.mean(np.asarray(abundance)[windwarf])
                    if not len(windwarf) < 30:
                        trendfehdwarfs.append(y)
                        trendnltedwarfs.append(winmeandwarf)

                    wingiant = \
                    np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) <= 3.5))[0]
                    if not len(wingiant) < 30:
                        winmeangiant = np.mean(np.asarray(abundance)[wingiant])

                        trendfehgiants.append(y)
                        trendnltegiants.append(winmeangiant)


                # plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(dwarfcount)+")")
                # plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(giantcount)+")")

                trenddict[element][lteval] = [trendnltedwarfs, trendnltegiants, trendfehdwarfs, trendfehgiants]

                if lteval:
                    ltestring = "LTE"
                    if element == "Ti1":
                        elementstring = "Ti I"
                        col = "blue"

                    else:
                        col = "blue"
                        elementstring = "Ti II"

                else:
                    ltestring = "NLTE"
                    if element == "Ti1":
                        col = "black"
                        elementstring = "Ti I"

                    else:
                        col = "black"
                        elementstring = "Ti II"
                print(lteval, trenddict[element][lteval][2], trenddict[element][lteval][0], "\n", trenddict[element][lteval][3], trenddict[element][lteval][1])
                print(trenddict[element][lteval][0][-3] - trenddict[element][lteval][1][-3], "\n")

                if startype == "dwarf":
                    sns.lineplot(x=trenddict[element][lteval][2], y=trenddict[element][lteval][0], label=elementstring+" "+ ltestring+" dwarfs ("+str(len(dwarfs[0]))+")", color = col)
                else:
                    sns.lineplot(x=trenddict[element][lteval][3], y=trenddict[element][lteval][1], linestyle="--", color=col,
                             label=elementstring+" " + ltestring+" giants ("+str(len(dwarfs[0]))+")")
                print(feh, abundance)
                den = ax.scatter_density(np.asarray(feh)[dwarfs], np.asarray(abundance)[dwarfs], cmap=LinearSegmentedColormap.from_list("white_viridis", [
                                                     (0, '#ffffff'),
                                                     (1e-20, '#89f089'),
                                                     (0.15, '#CDDD1A'),
                                                     (1, '#e21616'),
                                                 ], N=1000))

                plt.colorbar(den, label="Number of stars")
                ax.set_ylim([-0.6, 0.6])

                #plt.title("Titanium abundance trend lines")
                plt.legend(loc="lower left")
        #plt.savefig("graphs/Split_GCE_Trendlines_fe2_"+startype+".png", bbox_inches="tight")
        plt.show()

def galactic_evolution_trendlines_fe_splitup_testheat():
    trenddict = {}
    fig = plt.figure()
    fig.set_figwidth(15)
    fig.set_figheight(10)
    for element in ["Ti1", "Ti2"]:

        trenddict[element] = {}
        for lteval in [False, True]:
            if element == "Ti1":
                if lteval:
                    ax = fig.add_subplot(2, 2, 1, projection="scatter_density")
                    plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ I}$")

                else:
                    ax = fig.add_subplot(2, 2, 2, projection="scatter_density")

            else:
                if lteval:
                    ax = fig.add_subplot(2, 2, 3, projection="scatter_density")
                    plt.xlabel("Metallicity [Fe/H]")
                    plt.ylabel("[Ti/Fe]$_\mathrm{Ti\ II}$")

                else:
                    ax = fig.add_subplot(2, 2, 4, projection="scatter_density")
                    plt.xlabel("Metallicity [Fe/H]")

            b = pkl.load(open(element + "_Galah_impacts2_svenlines.pkl", "rb"))
            a = pkl.load(open(element + "_NLTE_list2_svenlines.pkl", "rb"))
            lineless = pkl.load(open(element + "_Galah_impacts_lineless2_svenlines.pkl", "rb"))
            sns.set_theme(style="whitegrid")

            hr = ((np.where(np.isnan(a[1])))[0])
            print(len(a[0]))
            print(len(b))
            fehgiants = []
            fehdwarfs = []
            nltegiants = []
            nltedwarfs = []
            dwarfcount = 0
            giantcount = 0
            index = 0

            trendnltedwarfs = []
            trendfehdwarfs = []
            trendnltegiants = []
            trendfehgiants = []

            # Checking numbero f dwarfs and giants by counting those with more or less than 3.5 logg ([1] in lineless)
            for g in [i[1] for i in lineless.values()]:
                if g <= 3.5:
                    giantcount += 1
                elif g > 3.5:
                    dwarfcount += 1
            # b is [starn][line] = teff, logg, feh, vmic, lte, nlte, nlte impact

            feh = [i[2] for i in lineless.values()]
            logg = [i[1] for i in lineless.values()]
            # turn from /H to /Fe by minusing metallicity fe/h
            if lteval:
                abundance = np.asarray([i[-2] for i in lineless.values()])-np.asarray(feh)
                if element == "Ti1":
                    abundanceuse = np.asarray([i[-2] for i in lineless.values()])-np.asarray(feh)
                    fehuse = feh
            else:
                abundance = np.asarray([i[-1] + i[-2] for i in lineless.values()])-np.asarray(feh)
            # sns.kdeplot(x=feh, y=nlte, fill=True, cmap="viridis")
            # Checking every 0.2 metallicity bin so we have a good averaged trend line

            for y in np.arange(min(feh), max(feh), 0.2):
                windwarf = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) > 3.5))[0]
                winmeandwarf = np.mean(np.asarray(abundance)[windwarf])
                if not len(windwarf) < 30:
                    trendfehdwarfs.append(y)
                    trendnltedwarfs.append(winmeandwarf)

                wingiant = \
                np.where(np.logical_and(np.logical_and(y - 0.1 <= feh, feh <= y + 0.1), np.asarray(logg) <= 3.5))[0]
                if not len(wingiant) < 30:
                    winmeangiant = np.mean(np.asarray(abundance)[wingiant])

                    trendfehgiants.append(y)
                    trendnltegiants.append(winmeangiant)


            # plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(dwarfcount)+")")
            # plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(giantcount)+")")

            trenddict[element][lteval] = [trendnltedwarfs, trendnltegiants, trendfehdwarfs, trendfehgiants]

            if lteval:
                ltestring = "LTE"
                if element == "Ti1":
                    elementstring = "Ti I"
                    col = "blue"

                else:
                    col = "blue"
                    elementstring = "Ti II"

            else:
                ltestring = "NLTE"
                if element == "Ti1":
                    col = "black"
                    elementstring = "Ti I"

                else:
                    col = "black"
                    elementstring = "Ti II"
            print(lteval, trenddict[element][lteval][2], trenddict[element][lteval][0], "\n", trenddict[element][lteval][3], trenddict[element][lteval][1])
            print(trenddict[element][lteval][0][-3] - trenddict[element][lteval][1][-3], "\n")

            sns.lineplot(x=trenddict[element][lteval][2], y=trenddict[element][lteval][0], label=elementstring+" "+ ltestring+" dwarfs", color = col)
            sns.lineplot(x=trenddict[element][lteval][3], y=trenddict[element][lteval][1], linestyle="--", color=col,
                         label=elementstring+" " + ltestring+" giants")
            print(feh, abundance)


            den = ax.scatter_density(feh, abundance, cmap=LinearSegmentedColormap.from_list("white_viridis", [
                                                 (0, '#ffffff'),
                                                 (1e-20, '#befeff'),
                                                 (0.3, '#38dddf'),
                                                    (0.7, "#1d7dcd"),
                                                 (1, '#1d27cd'),
                                             ], N=1000))

            ax.set_ylim([-0.6, 0.6])
            plt.legend(loc="lower left")

    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(den, cax=cbar_ax)
    fig.subplots_adjust(right=0.8)
    #plt.title("Titanium abundance trend lines")
    #plt.colorbar(den, label="Number of stars")

    ##plt.savefig("graphs/Split_GCE_Trendlines_fe2.png", bbox_inches="tight")
    plt.show()

ion_imbalance_fe_split()