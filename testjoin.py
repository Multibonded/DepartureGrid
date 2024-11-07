import os
import sys
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from scipy.interpolate import interpn
from scipy.interpolate import interp1d
import warnings
"""Here we make the nlte impact on abnundance by interpolating n/lte EWs on the atmospheric grid,"""
# file,element,teff,logg,feh,xi = sys.argv[1],sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),float(sys.argv[6])

star = "Sun"
file = star + "_Ti1_ew"
element = "Ti1"

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
ews = data["ew"]
ers = data["err"]
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
ag = grid["ag"] + solar + feh

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

fourcount = 0
fivecount = 0

for windex in index:
    alte = -9
    anlte = -9
    elte = 9
    enlte = 9

    ew = ews[windex]
    print("EW?", ew, len(ews), np.log10(ew))
    er = ers[windex]
    print("fin", wl[np.isfinite(wl)])

    wli = wl[:, :, :, :, :, windex]
    wni = wn[:, :, :, :, :, windex]
    print("fin", wli[np.isfinite(wli)])
    wla = []
    wna = []
    print("ag", ag)
    for temp in range(len(wn)):
        tempval = tg[temp]
        abund_dict[tempval] = {}

        for grav in range(len(wn[0])):
            graval = gg[grav]
            abund_dict[tempval][graval] = {}
            for met in range(len(wn[0][0])):
                metval = fg[met]
                abund_dict[tempval][graval][metval] = {}

                for vmic in range(len(wn[0][0][0])):
                    vmicval = xg[vmic]
                    abund_dict[tempval][graval][metval][vmicval] = {}

                    wavelength = line[windex]
                    abund_dict[tempval][graval][metval][vmicval][wavelength] = {}

                    nlte_ew = wn[temp, grav, met, vmic, :, windex]
                    lte_ew = wl[temp, grav, met, vmic, :, windex]
                    if -98 in nlte_ew or np.isnan(nlte_ew[0]):
                        delta = np.array(len(ag))+np.nan
                        if metval == -4:
                            fourcount +=1
                        if metval == -5:
                            fivecount+=1
                        print(delta, metval, 10**wl[temp, grav, met, vmic, :, windex] ,graval, tempval, 10**wn[temp, grav, met, vmic, :, windex])

                        continue
                    if max(nlte_ew)-min(nlte_ew) < 0.00001:
                        delta = np.array(len(ag))+np.nan
                        if metval == -4:
                            fourcount +=1
                        if metval == -5:
                            fivecount+=1
                        print(delta, metval, 10**wl[temp, grav, met, vmic, :, windex] ,graval, tempval, 10**wn[temp, grav, met, vmic, :, windex])

                        continue

                    inter = interp1d(wn[temp, grav, met, vmic, :, windex], ag,fill_value="extrapolate")
                    interl = inter(wl[temp, grav, met, vmic, :, windex])
                    delta = interl-ag
                    # Some parts messed up and have identical values at all abundances, clearly wrong but didn't catch
                    # a nan or -98. Seem to ave 7 and 10 values at -5 adnd -4 respectively
                    if np.inf in (delta) or -np.inf in delta:
                        delta = np.array(len(ag))+np.nan
                        print(delta, metval, 10**wl[temp, grav, met, vmic, :, windex] ,graval, tempval, 10**wn[temp, grav, met, vmic, :, windex])

                        if metval == -4:
                            fourcount +=1
                        if metval == -5:
                            fivecount+=1

                        continue
                    # in two places we see at 1.79 abund we have a random value of -6 delta fir no reason apparently
                    # at met -5 grav 3 temp5500
                    #
                    if min(delta) < -6:
                        print(delta, metval, 10**wl[temp, grav, met, vmic, :, windex] ,graval, tempval, 10**wn[temp, grav, met, vmic, :, windex])

                    for enhance in range(len(ag)):
                        abund_dict[tempval][graval][metval][vmicval][ag[enhance]][line[windex]] =interl[enhance]-ag[enhance]
                        #abund_dict[tempval][graval][metval][vmicval][line[windex]] =interl-ag

                    if parameters:
                        if metval == -5:
                            plt.scatter(tempval, graval, marker="x", c="red")
                        elif metval == -4:
                            plt.scatter(tempval, graval, marker="o", alpha=0.1, c="blue")

                    else:
                        # plot nlte impact
                        plt.scatter(np.ones(len(ag))*metval,  (interl-ag))#

                #plt.scatter(metval, (interl-ag)[-1])

if parameters:
    plt.xlabel("Teff")
    plt.ylabel("Log(g)")
else:
    plt.xlabel("Metallicity")
    plt.ylabel("NLTE impact on abundace")
print("number of errors in atmospgeres with - 4 and -5 met", fourcount, fivecount)

plt.show()
pkl.dump(abund_dict, open("nlte_impact_dict.pkl", "wb"))
