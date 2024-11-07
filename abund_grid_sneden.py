"""Makes a dictionary with all NLTE impacts for the grid nodes we have for all lines. Effects found by interpolating
The equivalent widths in NLTE to the NLTE abundances, and then use that function on the LTE equivalent widths and
taking the difference."""

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


reduced = False

element = "Ti2"
file = element+"_lines_yongsneden"
just_plot_simple = False
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
"""if reduced:
    pklfile = open('grids/{}_sneden_reduced.pkl'.format(element), 'rb')
else:
    pklfile = open('grids/{}_sneden.pkl'.format(element), 'rb')
"""
pklfile = open('grids/{}_sneden.pkl'.format(element), 'rb')
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

print(tg, gg, fg)
print("suf", tg, gg, fg)

print("B", wn.shape)
print("E", wn[np.isfinite(wn)], len(wn[np.isfinite(wn)]))
print("F", wl[np.isfinite(wl)], len(wl[np.isfinite(wl)]))
print("F2", len(wl[np.isnan(wl)]))
print("G", wl[0][0][0][0][0])
print("H", np.where(wl == -99), len(np.where(wl == -99)))
abund_dict = {}
ew_dict = {}

fourcount = 0
fivecount = 0
fourlist=[]
fivelist = []
for windex in index:
    alte = -9
    anlte = -9
    elte = 9
    enlte = 9
    print(line[windex])

    wli = wl[:, :, :, :, :, windex]
    wni = wn[:, :, :, :, :, windex]
    wla = []
    wna = []

    for temp in range(len(wn)):
        tempval = tg[temp]
        if tempval not in abund_dict:
            abund_dict[tempval] = {}
            ew_dict[tempval] = {}
        for grav in range(len(wn[0])):
            graval = gg[grav]
            if graval not in abund_dict[tempval]:
                abund_dict[tempval][graval] = {}
                ew_dict[tempval][graval] = {}

            for met in range(len(wn[0][0])):
                metval = fg[met]

                log_eps = np.round(ag + metval + 4.9, 4)
                if metval not in abund_dict[tempval][graval]:
                    abund_dict[tempval][graval][metval] = {}
                    ew_dict[tempval][graval][metval] = {}

                for vmic in range(len(wn[0][0][0])):
                    vmicval = xg[vmic]
                    if vmicval not in abund_dict[tempval][graval][metval]:
                        abund_dict[tempval][graval][metval][vmicval] = {}
                        ew_dict[tempval][graval][metval][vmicval] = {}

                    wavelength = line[windex]

                    nlte_ew = wn[temp, grav, met, vmic, :, windex]
                    lte_ew = wl[temp, grav, met, vmic, :, windex]
                    if -98 in nlte_ew or np.isnan(nlte_ew[0]):
                        delta = np.zeros(len(log_eps))+np.nan
                        if metval == -4:
                            fourcount +=1
                        if metval == -5:
                            fivecount+=1
                        #print(delta, metval, 10**wl[temp, grav, met, vmic, :, windex] ,graval, tempval, 10**wn[temp, grav, met, vmic, :, windex])

                    # There should be an EW change when abundance/enhancement changes. I can't see why there
                    # would be no change throughout, except some weird numerical error.
                    elif max(nlte_ew)-min(nlte_ew) < 0.00001:
                        delta = np.zeros(len(log_eps))+np.nan
                        if metval == -4:
                            fourcount +=1
                        if metval == -5:
                            fivecount+=1

                    else:
                        """print("ttt", temp, grav, met, vmic, log_eps, windex)
                        inter = interp1d(wn[temp, grav, met, vmic, :, windex], log_eps,fill_value="extrapolate", kind="linear")
                        interl = inter(wl[temp, grav, met, vmic, :, windex])
                        delta = interl-log_eps
                        print("delta!", delta)"""
                        try:
                            inter = interp1d(wn[temp, grav, met, vmic, :, windex], log_eps,fill_value="extrapolate", kind="linear")
                            interl = inter(wl[temp, grav, met, vmic, :, windex])
                            delta = interl - log_eps

                        except ValueError:

                            u, indices = np.unique(wn[temp, grav, met, vmic, :, windex], return_index=True)
                            interl = inter(wl[temp, grav, met, vmic, :, windex])
                            delta = interl - log_eps

                            try:
                                inter = interp1d(wn[temp, grav, met, vmic, :, windex][indices], log_eps[indices],fill_value="extrapolate", kind="cubic")
                            # If there aren't enough unique values we can't interpolate cubicly and I don't think
                            # it has a scientific result anyway, otherwise we could just use linear interpolation.
                            except ValueError:
                                delta = np.zeros(len(log_eps)) + np.nan

                    # Some parts messed up and have identical values at all abundances, clearly wrong but didn't catch
                    # a nan or -98. Seem to ave 7 and 10 values at -5 adnd -4 respectively

                    if np.inf in (delta) or -np.inf in delta:
                        delta = np.zeros(len(log_eps))+np.nan

                        if metval == -4:
                            fourcount +=1
                        if metval == -5:
                            fivecount+=1


                    # in two places we see at 1.79 abund we have a random value of -6 delta fir no reason apparently
                    # at met -5 grav 3 temp5500
                    # Some values are MESSED up, one goes up to 1E5!!! Insanity and not sure why.
                    if min(delta) < -6:
                        print("weird", delta, metval, 10**wl[temp, grav, met, vmic, :, windex] ,graval, tempval, 10**wn[temp, grav, met, vmic, :, windex])
                        delta = np.zeros(len(log_eps)) + np.nan

                    for enhance in range(len(log_eps)):
                        if log_eps[enhance] not in abund_dict[tempval][graval][metval][vmicval]:

                            abund_dict[tempval][graval][metval][vmicval][log_eps[enhance]] = {}
                            ew_dict[tempval][graval][metval][vmicval][log_eps[enhance]] = {}

                        try:
                            abund_dict[tempval][graval][metval][vmicval][log_eps[enhance]][line[windex]] =delta[enhance]
                            ew_dict[tempval][graval][metval][vmicval][log_eps[enhance]][line[windex]] =lte_ew[enhance]
                        except KeyError:
                            #print("line2", line[windex], ag[enhance])
                            abund_dict[tempval][graval][metval][vmicval][log_eps[enhance]][line[windex]] = {}
                            abund_dict[tempval][graval][metval][vmicval][log_eps[enhance]][line[windex]] =delta[enhance]
                            ew_dict[tempval][graval][metval][vmicval][log_eps[enhance]][line[windex]] = {}
                            ew_dict[tempval][graval][metval][vmicval][log_eps[enhance]][line[windex]] = lte_ew[enhance]
                    if not np.isnan(delta).any():

                        if parameters:
                            if metval == -5:
                                plt.scatter(tempval, graval, marker="x", c="red")
                            elif metval == -4:

                                plt.scatter(tempval, graval, marker="o", alpha=0.1, c="blue")

                        else:
                            # plot nlte impact
                            plt.scatter(np.ones(len(log_eps))*metval,  (interl-log_eps))#

                #plt.scatter(metval, (interl-ag)[-1])

if parameters:
    plt.xlabel("Teff")
    plt.ylabel("Log(g)")
    plt.title("Low metallicities of -4 and -5")
    plt.savefig(element+r"Teff and Log(g) grid available at low metallicities.png")

else:
    plt.xlabel("Metallicity")
    plt.ylabel("NLTE impact on abundance")
    plt.savefig(element+r"NLTE at different metallicities.png")
print("number of errors in atmospgeres with - 4 and -5 met", fourcount, fivecount)

"""if reduced:
    pkl.dump(abund_dict, open(element+"nlte_impact_dict_sneden_reduced.pkl", "wb"))
else:
    pkl.dump(abund_dict, open(element + "nlte_impact_dict_sneden.pkl", "wb"))
"""
pkl.dump(abund_dict, open(element + "nlte_impact_dict_sneden.pkl", "wb"))
pkl.dump(ew_dict, open(element + "ew_dict_sneden.pkl", "wb"))
plt.show()