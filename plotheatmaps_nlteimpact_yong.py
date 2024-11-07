import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl





def plot_basic_corrections(element):
    b = pkl.load(open(element + "_Yong_impacts.pkl", "rb"))
    a = pkl.load(open(element + "_Yong_NLTE_list.pkl", "rb"))

    print("Number of stars", len(a[0]), len(a[1]))
    plt.scatter(a[0], a[1], s=0.5)
    print(min(a[1]))


    plt.title("Mean NLTE corrections to Yong for "+element+ " averaged over all lines")
    plt.xlabel("Metallicity [Fe/H]")
    plt.ylabel("[Ti/Fe]$_{NLTE}$ - [Ti/Fe]$_{LTE}$")


    plt.savefig("graphs/"+element+"NLTE_Corrections.jpg")
    plt.show()

def plot_both():
    for element in ["Ti1", "Ti2"]:

        b = pkl.load(open(element + "_Yong_impacts.pkl", "rb"))
        a = pkl.load(open(element + "_Yong_NLTE_list.pkl", "rb"))

        print("Number of stars for ", element, ": ", len(a[1]))
        if element == "Ti1":
            plt.scatter(a[0], a[1], s=0.5, label = "Ti I")
        elif element == "Ti2":
            plt.scatter(a[0], a[1], s=0.5, label = "Ti II")

        print(min(a[1]))
        plt.legend()
        plt.title("Mean NLTE corrections to Yong averaged over all lines")
        plt.xlabel("Metallicity [Fe/H]")
        plt.ylabel("[Ti/Fe]$_{NLTE}$ - [Ti/Fe]$_{LTE}$")

    plt.savefig("graphs/Both_NLTE_Corrections.jpg")
    plt.show()
#plot_both()
import mpl_scatter_density # adds projection='scatter_density'
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Colormap


def color_single(element):
    b = pkl.load(open(element + "_Yong_impacts.pkl", "rb"))
    a = pkl.load(open(element + "_Yong_NLTE_list.pkl", "rb"))
    c = pkl.load(open(element + "_Yong_impacts_lineless.pkl", "rb"))

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


    comparelte = False
    if comparelte:
        plt.scatter(lte, nlte, s=0.5)
        plt.xlabel("LTE preediction by Yong")
        plt.ylabel("NLTE correction")
        plt.show()
        print("Exiting due to comparelte variable")
        exit()


    print((np.where(np.isnan(nlte))))
    print("Number of stars for ", element, ": ", len(nlte))
    ax = plt.figure().add_subplot(1,1,1, projection="scatter_density")
    den = ax.scatter_density(feh, nlte, cmap=LinearSegmentedColormap.from_list("white_viridis",  [
        (0, '#ffffff'),
        (1e-20, '#1ADD1A'),
        (0.15, '#CDDD1A'),
        (1, '#e21616'),
    ], N=1000))
    #den = ax.scatter_density(feh, nlte, cmap=LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N).reversed())

    #den = ax.scatter_density(feh, nlte, cmap="Reds")

    plt.colorbar(den, label="Number of stars")
    if element == "Ti1":
        plt.title("Mean NLTE corrections to Yong for Ti I averaged over all lines")
    elif element == "Ti2":
        plt.title("Mean NLTE corrections to Yong for Ti II averaged over all lines")


    plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarf (N="+str(len(nltedwarfs))+")")
    plt.plot(trendfehgiants, trendnltegiants, "--",  c="blue", label = "Giant (N="+str(len(nltegiants))+")")
    plt.legend()
    plt.xlabel("Metallicity [Fe/H]")
    plt.ylabel("[Ti/Fe]$_{NLTE}$ - [Ti/Fe]$_{LTE}$")
    plt.savefig("graphs/NLTEColormap_"+element)
    plt.show()

def color_both():
    for element in ["Ti1", "Ti2"]:
        print("Not working yet, exiting")
        exit()
        b = pkl.load(open(element + "_Yong_impacts.pkl", "rb"))
        a = pkl.load(open(element + "_Yong_NLTE_list.pkl", "rb"))
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
            (1e-20, '#1ADD1A'),
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
    plt.title("Mean NLTE corrections to Yong for Ti averaged over all lines")

    plt.xlabel("Metallicity [Fe/H]")
    plt.ylabel("[Ti/Fe]$_{NLTE}$ - [Ti/Fe]$_{LTE}$")
    plt.show()
color_single("Ti1")