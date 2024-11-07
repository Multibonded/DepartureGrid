import pickle as pkl

element = "Ti1"
reliable     = pkl.load(open(element+"Accurate_param_indexes_galah.pkl", "rb"))
element = "Ti2"
reliable2     = pkl.load(open(element+"Accurate_param_indexes_galah.pkl", "rb"))


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

logg = atom1['logg']
c2 = atom1['snr_c2_iraf']
lines = ["4758", "4759", "4782", "4820", "5739"]
lines2 = ["4720", "4765", "4799", "4866"]
lines = ["4758", "4759", "4782", "4798", "4820", "5689", "5716", "5720", "5739"]
lines2 = ["4874", "4798", "4765"]

lines = ["4758", "4759", "4782", "5739", "4798"]
lines2 = ["4798", "4765", "4720"]



def imbalance_and_abunds():
    stardict = {}
    stardict2 = {}
    ti1giants = []
    ti2giants = []
    ti1dwarfs = []
    ti2dwarfs = []
    stars = []
    stardwarf = []
    ti_galah = []
    ti_galah_dwarfs = []

    temp = []
    count = 0
    errdict = {}
    errdict2 = {}
    for line in lines:
        errdict[line] = atom1['ind_cov_e_Ti' + line]
        stardict[line] = atom1['ind_Ti' + line + '_fe']
    for line in lines2:
        errdict2[line] = atom1['ind_cov_e_Ti' + line]
        stardict2[line] = atom1['ind_Ti' + line + '_fe']


    for starn in reliable:
        count+=1

        if count % 500 == 0:
            print(count, "/", len(reliable))
        if starn not in reliable2:
            continue
        starl=[]
        stare = []
        starl2=[]
        stare2 = []
        for line in lines:
            if not np.isnan(stardict[line][starn]):
                starl.append(stardict[line][starn])
                stare.append(errdict[line][starn])
        for line in lines2:
            if not np.isnan(stardict2[line][starn]):

                starl2.append(stardict2[line][starn])
                stare2.append(errdict2[line][starn])
        if len(starl2) <2 :
            continue
        if len(starl) < 2:
            continue
        if starn in reliable2:
            if logg[starn] <= 3.5:
                stars.append(starn)
                avg1 = np.average(starl, weights = stare)
                avg2 = np.average(starl2, weights = stare2)
                ti_galah.append(avg1-avg2)
                ti1giants.append(avg1)
                ti2giants.append(avg2)
            else:
                stardwarf.append(starn)
                avg1 = np.average(starl, weights = stare)
                avg2 = np.average(starl2, weights = stare2)
                ti_galah_dwarfs.append(avg1-avg2)
                ti1dwarfs.append(avg1)
                ti2dwarfs.append(avg2)

    feh_galah = atom1[stars]['fe_h']
    feh_dwarfs = atom1[stardwarf]['fe_h']
    feh = atom1['fe_h']
    #ti_galah = atom1[stars]['ti_fe']-atom1[stars]['ti2_fe']


    feh = []
    ti1= []
    ti2=[]
    ti = []

    trendfehgiants = []
    trendnltegiants = []

    ti1giantplot = []
    ti2giantplot = []
    print(len(feh_galah))
    for y in np.arange(min(feh_galah), max(feh_galah), 0.2):

        wingiant = \
            np.where((np.logical_and(y - 0.1 <= feh_galah, feh_galah <= y + 0.1)))[0]
        print("len", len(wingiant))
        if not len(wingiant) < 30:
            winmeangiant = np.nanmean(np.asarray(ti_galah)[wingiant])
            ti1giantmean = np.nanmean(np.asarray(ti1giants)[wingiant])
            ti2giantmean = np.nanmean(np.asarray(ti2giants)[wingiant])
            trendfehgiants.append(y)
            trendnltegiants.append(winmeangiant)
            ti1giantplot.append(ti1giantmean)
            ti2giantplot.append(ti2giantmean)

    trendfehdwarfs = []
    trendnltedwarfs = []

    ti1dwarfplot = []
    ti2dwarfplot = []

    for y in np.arange(min(feh_dwarfs), max(feh_dwarfs), 0.2):

        windwarf = \
            np.where((np.logical_and(y - 0.1 <= feh_dwarfs, feh_dwarfs <= y + 0.1)))[0]
        print("len", len(windwarf))
        if not len(windwarf) < 30:
            winmeandwarf= np.nanmean(np.asarray(ti_galah_dwarfs)[windwarf])
            ti1dwarfmean = np.nanmean(np.asarray(ti1dwarfs)[windwarf])
            ti2dwarfmean = np.nanmean(np.asarray(ti2dwarfs)[windwarf])

            trendfehdwarfs.append(y)
            trendnltedwarfs.append(winmeandwarf)
            ti1dwarfplot.append(ti1dwarfmean)
            ti2dwarfplot.append(ti2dwarfmean)

    print(trendfehgiants, trendnltegiants)
    print(len(stars), len(ti_galah))
    plt.scatter(feh_galah, ti_galah, s=2)
    plt.scatter(feh_dwarfs, ti_galah_dwarfs, s=2)

    plt.plot(trendfehgiants, trendnltegiants, c="red", label = "Giants")
    plt.plot(trendfehdwarfs, trendnltedwarfs, c="black", label = "Dwarfs")
    plt.legend()
    plt.ylabel("Ti I - Ti II")
    plt.savefig("graphs/svencheck/ionimb_ours_test.png")
    plt.show()

    plt.plot(trendfehgiants, ti1giantplot, c = "blue", label = "Ti1 Giants")
    plt.plot(trendfehdwarfs, ti1dwarfplot, c = "red", label = "Ti1 Dwarfs")

    plt.plot(trendfehgiants, ti2giantplot, c = "black", label="Ti2 Giants")
    plt.plot(trendfehdwarfs, ti2dwarfplot, c = "green", label = "Ti2 Dwarfs")
    plt.legend()
    plt.ylabel("Ti/Fe")
    plt.savefig("graphs/svencheck/abunds_ours_test.png")
    plt.show()

def linescatter():
    for line in lines:
        print(line)
        stardict = {}
        ti1giants = []
        ti1dwarfs = []
        stars = []
        stardwarf = []
        stargiants=[]
        temp = []
        count = 0
        errdict = {}
        errdict2 = {}
        errdict[line] = atom1['ind_cov_e_Ti' + line]
        stardict[line] = atom1['ind_Ti' + line + '_fe']
        for starn in reliable:
            count += 1

            if count % 500 == 0:
                print(count, "/", len(reliable))
            if starn not in reliable2:
                continue
            starl = []
            stare = []
            if not np.isnan(stardict[line][starn]):
                starl.append(stardict[line][starn])
                stare.append(errdict[line][starn])
            if len(starl)<1:
                continue
            if starn in reliable2:
                if logg[starn] <= 3.5:
                    stars.append(starn)
                    avg1 = np.average(starl, weights=stare)
                    ti1giants.append(avg1)
                else:
                    stardwarf.append(starn)
                    avg1 = np.average(starl, weights=stare)
                    ti1dwarfs.append(avg1)

        feh_galah = atom1[stars]['fe_h']
        feh_dwarfs = atom1[stardwarf]['fe_h']
        feh = atom1['fe_h']
        # ti_galah = atom1[stars]['ti_fe']-atom1[stars]['ti2_fe']

        feh = []
        ti1 = []
        ti2 = []
        ti = []

        trendfehgiants = []

        ti1giantplot = []
        for y in np.arange(min(feh_galah), max(feh_galah), 0.2):

            wingiant = \
                np.where((np.logical_and(y - 0.1 <= feh_galah, feh_galah <= y + 0.1)))[0]
            if not len(wingiant) < 30:
                ti1giantmean = np.nanmean(np.asarray(ti1giants)[wingiant])
                trendfehgiants.append(y)
                ti1giantplot.append(ti1giantmean)

        trendfehdwarfs = []

        ti1dwarfplot = []

        for y in np.arange(min(feh_dwarfs), max(feh_dwarfs), 0.2):

            windwarf = \
                np.where((np.logical_and(y - 0.1 <= feh_dwarfs, feh_dwarfs <= y + 0.1)))[0]
            if not len(windwarf) < 30:
                ti1dwarfmean = np.nanmean(np.asarray(ti1dwarfs)[windwarf])

                trendfehdwarfs.append(y)
                ti1dwarfplot.append(ti1dwarfmean)



        plt.plot(trendfehgiants, ti1giantplot, c="blue", label="Giants mean")
        plt.plot(trendfehdwarfs, ti1dwarfplot, c="red", label="Dwarfs mean")
        plt.scatter(feh_galah, ti1giants, c="black", label="Giants")
        plt.scatter(feh_dwarfs, ti1dwarfs, c="green", label = "dwarfs")
        plt.title(line)
        plt.legend()

        plt.savefig("graphs/svencheck/"+str(line)+"line_ours_test.png")
        plt.show(block=False)
        plt.close()
imbalance_and_abunds()