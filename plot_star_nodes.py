"""Plotting Temp vs grav, with our nodes for interpolation and all reliable galah stars also plotted"""

import os
import sys
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from scipy.interpolate import interpn
from scipy.interpolate import interp1d
import warnings
import seaborn as sns
"""Here we make the nlte impact on abnundance by interpolating n/lte EWs on the atmospheric grid,"""
# file,element,teff,logg,feh,xi = sys.argv[1],sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),float(sys.argv[6])

element = "Ti1"
file = element+"_lines"
# Only use one line, as we already get indexes from the reliable file, so no need to check between different lines
just_plot_simple = True
if just_plot_simple:
    file = "Ti1_line_single"
"""
Star    Teff    Logg    Fe/H    Mic
Sun     5772    4.44    0       0.9
Arcturus 4286   1,64    -0.53   1.3
84937   6356    4.06    -2.06   1.2
140283  5792    3.65    -2.36   1.3
122563  4646    1.4     -2,5    1.8

"""
teff = 5792
logg = 3.65
feh = -2.36
xi = 1.3
# plot paramtetrs or nlte?
parameters = True
data = pd.read_csv(file, delim_whitespace=True)
index = data.index
line = data["line"]
print("line", line)
pklfile = open('grids/{}.pkl'.format(element), 'rb')
grid = pkl.load(pklfile)

if element == 'Na':
    solar = 6.17
if element == 'Mg':
    solar = 7.53
if element == 'Al':
    solar = 6.37
if element == "Ti" or element == "Ti1" or element == "Ti2":
    solar = 4.90
tg = grid["tg"]
gg = grid["gg"]
xg = np.array(grid["xg"])
fg = np.array(grid["fg"])
ag = grid["ag"]
ag = [round(x, 4) for x in ag]
wg = grid["wg"]
wl = grid["wl"]
wn = grid["wn"]
ti = np.where(tg > teff)
ti = min(ti[0]) - 1
gi = np.where(gg > logg)
gi = min(gi[0]) - 1
fi = np.where(fg > feh)
fi = min(fi[0]) - 1
print(tg, gg, fg)

print("suf", tg, gg, fg)

print("A", [ti, ti + 2])
print("B", wn.shape)
print("D",ti)

print("D", ti)
print("E", wn[np.isfinite(wn)], len(wn[np.isfinite(wn)]))
print("F", wl[np.isfinite(wl)], len(wl[np.isfinite(wl)]))
print("F2", len(wl[np.isnan(wl)]))
print("G", wl[0][0][0][0][0])
print("H", np.where(wl == -99), len(np.where(wl == -99)))
abund_dict = {}

fig, axs = plt.subplots(2, 2)

fourcount = 0
fivecount = 0
# lists to plot for different metallicity values of stars
fourlist=[]
threelist = []
twolist = []
onelist = []

from astropy.io import fits

atomname = r"/home/jama2357/Documents/TiFeAtoms/DepartureGrid/GALAH_DR3_main_allspec_v2.fits"
atom = fits.open(atomname)
print(repr(atom[1].header))
atom1 = atom[1].data
# now plot galah paramaters using the msot accurate stars
onelist = []
twolist =[]
threelist = []
fourlist = []
oneval = 0
twoval = 0
threeval = 0
fourval = 0
reliable1 = pkl.load(open("Ti1Accurate_param_indexes_galah.pkl", "rb"))
reliable2 = pkl.load(open("Ti2Accurate_param_indexes_galah.pkl", "rb"))
reliable = reliable1+reliable2
reliable=set(reliable)
reliablefinal = pkl.load(open("final_stars_used.pkl", "rb"))

stars = []
for star in reliablefinal:
    if star not in reliable:
        print("not in", star)
        stars.append(star)
print("number not in", len(stars))


nancount = 0
for starn in (reliablefinal):
    metval = atom1[starn]['fe_h']
    tempval = atom1[starn]['teff']
    graval = atom1[starn]['logg']
    if np.isnan(atom1[starn]['ti_fe']) and np.isnan(atom1[starn]['ti2_fe']):
        nancount+=1
        print("What", starn, atom1[starn]['ti_fe'])
    if 0 < metval <= 1:
        onelist.append([tempval, graval])
        oneval += 1
    elif -0.5 < metval <= 0:
        twolist.append([tempval, graval])
        twoval += 1

    elif -1.5 < metval <= -0.5:
        threelist.append([tempval, graval])
        threeval += 1

    elif -3.6 < metval <= -1.5:
        fourlist.append([tempval, graval])
        fourval += 1

print("nan count", nancount)
axs[0, 0].scatter(([x[0] for x in onelist]), ([y[1] for y in onelist]), label="0 < [Fe/H] <= 1.0 ("+f"{oneval:,d}"+" stars)", c="black", s=1.5, alpha=0.2)
axs[0, 1].scatter(([x[0] for x in twolist]), ([y[1] for y in twolist]), label="-0.5 < [Fe/H] <= 0 ("+f"{twoval:,d}"+" stars)", c="black", s=1.5, alpha=0.2)
axs[1, 0].scatter(([x[0] for x in threelist]), ([y[1] for y in threelist]), label="-1.5 < [Fe/H] <= -0.5 ("+f"{threeval:,d}"+" stars)", c="black", s=1.5, alpha=0.2)
axs[1, 1].scatter(([x[0] for x in fourlist]), ([y[1] for y in fourlist]), label="-3.6 < [Fe/H] <= -1.5 ("+f"{fourval:,d}"+" stars)", c="black", s=1.5, alpha=0.2)
onelist = []
twolist =[]
threelist = []
fourlist = []


realnodes = pkl.load(open("atmosphere_nodes.pkl", "rb"))
"""for windex in index:
    alte = -9
    anlte = -9
    elte = 9
    enlte = 9


    wli = wl[:, :, :, :, :, windex]
    wni = wn[:, :, :, :, :, windex]
    wla = []
    wna = []
    print(windex)
    for temp in range(len(wn)):
        tempval = tg[temp]
        if tempval not in abund_dict:
            abund_dict[tempval] = {}

        for grav in range(len(wn[0])):
            graval = gg[grav]
            if graval not in abund_dict[tempval]:
                abund_dict[tempval][graval] = {}
            for met in range(len(wn[0][0])):
                metval = fg[met]
                if metval not in abund_dict[tempval][graval]:
                    abund_dict[tempval][graval][metval] = {}
                if [tempval, graval, metval] not in realnodes:
                    continue
                for vmic in range(len(wn[0][0][0])):
                    vmicval = xg[vmic]
                    if vmicval not in abund_dict[tempval][graval][metval]:
                        abund_dict[tempval][graval][metval][vmicval] = {}

                    wavelength = line[windex]

                    nlte_ew = wn[temp, grav, met, vmic, :, windex]
                    lte_ew = wl[temp, grav, met, vmic, :, windex]
                    if -98 in nlte_ew or np.isnan(nlte_ew[0]):
                        delta = np.zeros(len(ag))+np.nan
                        if metval == -4:
                            fourcount +=1
                        if metval == -5:
                            fivecount+=1

                        #print(delta, metval, 10**wl[temp, grav, met, vmic, :, windex] ,graval, tempval, 10**wn[temp, grav, met, vmic, :, windex])


                    if max(nlte_ew)-min(nlte_ew) < 0.00001:
                        delta = np.zeros(len(ag))+np.nan
                        if metval == -4:
                            fourcount +=1
                        if metval == -5:
                            fivecount+=1

                    else:
                        inter = interp1d(wn[temp, grav, met, vmic, :, windex], ag,fill_value="extrapolate")
                        interl = inter(wl[temp, grav, met, vmic, :, windex])
                        delta = interl-ag
                    # Some parts messed up and have identical values at all abundances, clearly wrong but didn't catch
                    # a nan or -98. Seem to ave 7 and 10 values at -5 adnd -4 respectively

                    if np.inf in (delta) or -np.inf in delta:
                        delta = np.zeros(len(ag))+np.nan

                        if metval == -4:
                            fourcount +=1
                        if metval == -5:
                            fivecount+=1


                    # in two places we see at 1.79 abund we have a random value of -6 delta fir no reason apparently
                    # at met -5 grav 3 temp5500
                    #
                    if min(delta) < -6:
                        delta = np.zeros(len(ag)) + np.nan

                    for enhance in range(len(ag)):
                        if ag[enhance] not in abund_dict[tempval][graval][metval][vmicval]:

                            abund_dict[tempval][graval][metval][vmicval][ag[enhance]] = {}

                        try:
                            abund_dict[tempval][graval][metval][vmicval][ag[enhance]][line[windex]] =delta[enhance]
                        except KeyError:
                            #print("line2", line[windex], ag[enhance])
                            abund_dict[tempval][graval][metval][vmicval][ag[enhance]][line[windex]] = {}
                            abund_dict[tempval][graval][metval][vmicval][ag[enhance]][line[windex]] =delta[enhance]
                    if not np.isnan(delta).any():
                        
                        if tempval == 5250 and metval == -5.0:
                            #print("good!", graval, delta, metval)
                            pass
                        if parameters:
                            if metval == 0.5:
                                onelist.append([tempval, graval])

                            elif metval == 0:
                                twolist.append([tempval, graval])
                            elif metval == -1:
                                threelist.append([tempval, graval])

                            elif metval == -2.5:
                                fourlist.append([tempval, graval])

                        else:
                            # plot nlte impact
                            plt.scatter(np.ones(len(ag))*metval,  (interl-ag))

                #plt.scatter(metval, (interl-ag)[-1])
"""
for atmos in realnodes:
    temp = atmos[0]
    grav = atmos[1]
    feh = atmos[2]
    if feh == -2:
        fourlist.append([temp, grav])
    elif feh == -1:
        threelist.append([temp, grav])
    elif feh== 0:
        twolist.append([temp, grav])
    elif feh == 0.5:
        onelist.append([temp, grav])
axs[0, 0].scatter(([x[0] for x in onelist]), ([y[1] for y in onelist]), label="Nodes: [Fe/H] = 0.5", c="b", marker = "s",  s=5.5)
axs[0, 1].scatter(([x[0] for x in twolist]), ([y[1] for y in twolist]), label="Nodes: [Fe/H] = 0.0", c="b", marker = "s",  s=5.5)
axs[1, 0].scatter(([x[0] for x in threelist]), ([y[1] for y in threelist]), label="Nodes: [Fe/H] = -1.0",  marker = "s", c="b", s=5.5)
axs[1, 1].scatter(([x[0] for x in fourlist]), ([y[1] for y in fourlist]), label="Nodes: [Fe/H] = -2.5", marker = "s", c="b", s=5.5)



for ax in axs.flat:
    ax.axis(ymin=-1, ymax=6)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_xlabel('T$_{eff}$ / K', fontsize=18)
    ax.set_ylabel('log(g / cm s$^{-2}$)', fontsize=18)
    ax.legend()

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
print("number of errors in atmospgeres with - 4 and -5 met", fourcount, fivecount)
fig2 = plt.gcf()
fig2.set_size_inches(18.5, 10.5)
plt.savefig("graphs/StellarParamatersWithNodes.jpg", bbox_inches="tight")
plt.show()
