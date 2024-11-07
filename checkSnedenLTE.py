import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns

sns.set_theme(style="whitegrid")
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Colormap


def find_lines_in_all_stars():
    ti1lines = []
    ti2lines = []
    ti1all = []
    ti1dict = {}
    ti2all = []
    ti2dict = {}
    for filename in os.listdir("/home/jama2357/Documents/TiFeAtoms/DepartureGrid/CoG/sneden_stars/"):

        ti1this = []
        ti2this = []
        with open("/home/jama2357/Downloads/"+filename) as file:
            for line in file:

                splitline  = (line.split())
                if "/" in splitline[0]:
                    print(splitline)
                    feh = float(splitline[2])
                    teff = float(splitline[0].replace("/", ""))
                    logg = float(splitline[1].replace("/", ""))
                    vmic = float(splitline[-1])
                    if feh > -1.5:
                        break
                if splitline[0] =="Ti" and splitline[1] == "I":
                    line = float(splitline[2])
                    abund = float(splitline[7])
                    reduced = float(splitline[6])
                    enh = float(splitline[-1])
                    # in milliangstroms.
                    ew = float(splitline[5])

                    print("made red." , np.log10((ew/1000)/line))
                    print("      reduced, EW, teff, logg, feh, enhancement,    line")
                    print("reduced", reduced, ew, teff, logg, feh, enh, line, filename)
                    if reduced > -4.9:
                        continue

                    if line not in ti1dict:
                        ti1dict[line] = {}

                        ti1dict[line]["feh"] = [feh]
                        ti1dict[line]["abund"] = [abund]
                    else:
                        ti1dict[line]["feh"].append(feh)
                        ti1dict[line]["abund"].append(abund)
                elif splitline[0] =="Ti" and splitline[1] == "II":
                    line = float(splitline[2])
                    abund = float(splitline[7])

                    if line not in ti2dict:
                        ti2dict[line] = {}

                        ti2dict[line]["feh"] = [feh]
                        ti2dict[line]["abund"] = [abund]
                    else:
                        ti2dict[line]["feh"].append(feh)
                        ti2dict[line]["abund"].append(abund)
    return ti1dict, ti2dict
def plot_wls(tidict, tistring):

    import matplotlib.pyplot as plt
    maxnum = -100
    minnum = 100
    maxfeh = -100
    minfeh=100
    for line in tidict.keys():
        numb = len(tidict[line]['feh'])
        feh = tidict[line]['feh']
        if numb > maxnum:
            maxnum = numb
        if numb < minnum:
            minnum = numb
        if max(feh) > maxfeh:
            maxfeh = max(feh)
        if min(feh) < minfeh:
            minfeh = min(feh)

    for line in tidict.keys():
        fehs = np.asarray(tidict[line]['feh'])
        abunds = np.asarray(tidict[line]['abund'])
        abund = np.mean(abunds)
        numb = len(tidict[line]['feh'])
        meanfeh = np.mean(tidict[line]['feh'])
        if numb < 10:
            pass
        """if meanfeh >-2:
            maxfeh = -2
            continue"""
        cmap = LinearSegmentedColormap.from_list("white_viridis", [
            (0, '#1f32db'),
            (1, '#e21616'),
        ], N=1000)
        #plt.scatter(line, abund, c=numb, vmin=minnum, vmax= maxnum)
        plt.scatter(line, abund, c=numb, vmin=minnum, vmax= maxnum, cmap=cmap)

        #plt.scatter(line, abund, c=meanfeh, vmin=minfeh, vmax=maxfeh)
    print(tidict)
    print(len(tidict))
    plt.ylim(2.3, 3.4)
    plt.xlim(3000, 5500)
    plt.colorbar(label="Number of stars")
    plt.legend()
    plt.xlabel("Line ($\mathrm{\AA}$)")
    plt.ylabel("A(Ti)")
    plt.title(tistring)

    plt.savefig("graphs/sneden/wl_number_"+tistring+".png")
    plt.show()


def plot_wl_all():
    ti1, ti2 = find_lines_in_all_stars()
    ti1tr = True

    for ti in [ti1, ti2]:
        if ti1tr:
            tistring = "Ti I"
        else:
            tistring = "Ti II"
        plot_wls(ti, tistring)
        ti1tr = False


def plotlte():
    fig, ax = plt.subplots(nrows=2, ncols=2)


    ti1lines = []
    ti2lines = []
    ti1all = []
    ti1dict = {}
    ti2all = []
    ti2dict = {}
    for filename in os.listdir("/home/jama2357/Documents/TiFeAtoms/DepartureGrid/CoG/sneden_stars/"):

        ti1this = []
        ti2this = []
        ti1dict[filename]= {}
        ti2dict[filename] = {}
        with open("/home/jama2357/Downloads/" + filename) as file:
            abunds1 = []
            abunds2 = []

            for line in file:
                splitline = (line.split())
                if "/" in splitline[0]:
                    feh = float(splitline[2])
                    ti1dict[filename]["feh"] = [feh]
                    ti2dict[filename]["feh"] = [feh]

                    if feh > -1.5:
                        pass
                if splitline[0] == "Ti" and splitline[1] == "I":
                    line = float(splitline[2])
                    abund = float(splitline[7])

                    ti1dict[filename][line] = abund
                    abunds1.append(abund)
                if splitline[0] == "Ti" and splitline[1] == "II":
                    line = float(splitline[2])
                    abund = float(splitline[7])

                    ti2dict[filename][line] = abund
                    abunds2.append(abund)
            ti1dict[filename]['abund'] = np.mean(abunds1)
            ti2dict[filename]['abund'] = np.mean(abunds2)
    for filename in ti1dict:

        feh = ti1dict[filename]['feh']
        feh2 = ti2dict[filename]['feh']

        abund = ti1dict[filename]['abund']-feh-4.9

        abund2 = ti2dict[filename]['abund']-feh2-4.9
        print(abund, abund2)
        ax[0,0].scatter(feh, abund)
        ax[0,1].scatter(feh2, abund2)
        ax[1,0].scatter(feh, abund-abund2)

    plt.show()
    return ti1dict, ti2dict
find_lines_in_all_stars()