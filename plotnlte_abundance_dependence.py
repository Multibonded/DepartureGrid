

import os
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
full_galah_nlte = {}
element = "Ti1"

stars = {}
stars['Sun'] = [5772, 4.44, 0.00, 1.0]

stars["Arcturus"] = [4286, 1.64, -0.53, 1.3]
#stars["arcturus"] = [4000.0, 2.0, -0.5, 1]

#stars[122563] = [4500.0, 1.5, -2.5, 2]
#stars['sun'] = [5000.0, 4.5, 0, 1]
stars[84937] = [6356, 4.06, -2.06, 1.2]

stars[140283] = [5792, 3.65, -2.36, 1.3]
#stars[140283] = [6000.0, 3.5, -2.5, 1]
#stars[84937] = [6500.0, 4.0, -2.0, 1]
stars[122563] = [4636, 1.40, -2.50, 1.8]



def specific_stars(element):
    pklfile = open('grids/{}_all.pkl'.format(element), 'rb')

    grid = pkl.load(pklfile)
    tg = grid["tg"]

    gg = grid["gg"]
    xg = np.array(grid["xg"])
    fg = np.array(grid["fg"])
    # This does not affect ionisation imbalance just definitive abundance predictions in nlte_interp. Unsure how
    # but as we're applying dep. coff. ontothe abund of GALAH I don't hink it matters here for us.
    # log_eps = grid['ag']
    # log_eps = [round(x, 4) for x in log_eps]

    log_eps = grid["ag"]  # + solar + feh

    log_eps = [round(x, 4) for x in log_eps]
    star_mean_dict = {}

    for abund in log_eps:

        redew_dict = {}
        counter = 1
        feh_list = []
        nlte_list = []

        abund_dict = pkl.load(open(element + "nlte_impact_dict_all.pkl", "rb"))
        ew_dict = pkl.load(open(element + "ew_dict_all.pkl", "rb"))
        star_count = pkl.load(open(element + "Accurate_param_indexes_galah.pkl", "rb"))

        # This one we only want star info, including mean abundance nlte impact, but not line info as that clutters
        # it up . For that we use the full Galah dict.
        # Finds the mean A(Ti) of titanium for solar twins so we can more accurately find [Ti/H] which is better than Ti/fe
        # for reasons. Given as dict[line] = [lte, nlte]
        if "1" in element:
            lines = ["4758", "4759", "4778", "4781", "4797", "4801", "4820", "5689", "5716", "5720", "5739", "5866"]
            lines = ["4758", "4759", "4781", "4797", "4820", "5689", "5716", "5720", "5739"]
            lines = ["4758", "4759", "5689", "5716", "5739"]


        elif "2" in element:
            lines = ["4874", "4865", "4849", "4798", "4764", "4719"]
            lines = ["4719", "4798", "4764", "4866"]
            lines = ["4798", "4719", "4874"]




        for starn in stars:
            if starn not in star_mean_dict:
                star_mean_dict[starn] = {}
            if starn not in redew_dict:
                redew_dict[starn] = {}
            counter += 1
            # number of lines (after we restrict them) per star, remove is less than 3
            linecount = 0
            lines_mean = []
            errors = []
            for line_index in range(len(lines)):
                line = int(lines[line_index])


                teff = (stars[starn][0])
                logg = (stars[starn][1])
                xi = (stars[starn][3])
                feh = (stars[starn][2])

                # paramS 6324.3784 3.9686449 1.4542778 -0.38920593
                pklfile = open('grids/{}_reduced.pkl'.format(element), 'rb')

                grid = pkl.load(pklfile)
                tg = grid["tg"]

                gg = grid["gg"]
                xg = np.array(grid["xg"])
                fg = np.array(grid["fg"])
                # This does not affect ionisation imbalance just definitive abundance predictions in nlte_interp. Unsure how
                # but as we're applying dep. coff. ontothe abund of GALAH I don't hink it matters here for us.
                #log_eps = grid['ag']
                #log_eps = [round(x, 4) for x in log_eps]

                log_eps = grid["ag"]# + solar + feh

                log_eps = [round(x, 4) for x in log_eps]

                ti = np.where(tg > teff)
                ti = min(ti[0]) - 1
                gi = np.where(gg > logg)
                gi = min(gi[0]) - 1
                fi = np.where(fg > feh)
                fi = min(fi[0]) - 1

                log_eps_i = np.where(log_eps > abund)
                if len(log_eps_i[0]) == 0:
                    log_eps_i = np.where(log_eps == abund)
                log_eps_i = min(log_eps_i[0]) - 1

                """tg = tg[ti:ti + 2]
                gg = gg[gi:gi + 2]
                fg = fg[fi:fi + 2]
                log_eps = log_eps[log_eps_i:log_eps_i + 2]"""
                tg = tg[ti:ti + 2]
                gg = gg[gi:gi + 2]
                fg = fg[fi:fi + 2]
                log_eps = log_eps[log_eps_i:log_eps_i + 2]

                abundance_corrections = np.zeros([len(tg), len(gg), len(fg), len(xg), len(log_eps)]) + np.nan
                ew_corrections = np.zeros([len(tg), len(gg), len(fg), len(xg), len(log_eps)]) + np.nan

                # MAkes a np array with all parameters in the correct slots to allow easy itnterpolation as the nlte corrections
                # are int he sane shape as the parameter list.
                for t in range(len(tg)):
                    for g in range(len(gg)):
                        for f in range(len(fg)):
                            for x in range(len(xg)):
                                for a in range(len(log_eps)):
                                    # Abund dict has the NLTE impact. Log eps is still just [Ti/Fe] at this point, change it to A(Ti) with +fg[f] + 4.9 (Star metalcitiy and solar abundance)
                                    abundance_corrections[t, g, f, x, a] = ([abund_dict[tg[t]][gg[g]][fg[f]][xg[x]][np.round(log_eps[a]+fg[f]+4.9, 4)][line]][0])
                                    ew=(([ew_dict[tg[t]][gg[g]][fg[f]][xg[x]][np.round(log_eps[a]+fg[f]+4.9, 4)][line]][0]))
                                    ew_corrections[t, g, f, x, a] = ew
                                    #print(line, t, g, f, x, a, [abund_dict[tg[t]][gg[g]][fg[f]][xg[x]][log_eps[a]][line]][0])
                                    """if starn == 59329:
                                        print("\nLine      Teff       Logg      Feh       Vmic       Enhance      NLTE effect")
                                        print(line, tg[t], "    ", gg[g], "    ", fg[f], "      ", xg[x], "    ",
                                              log_eps[a], "            ", round(abund_dict[tg[t]][gg[g]][fg[f]][xg[x]][log_eps[a]][line], 4))
                                    """
                # Log eps is still just [Ti/Fe] at this point, change it to A(Ti) with +fg[f] + 4.9
                nlte_impact = interpn([tg, gg, fg, xg, log_eps], abundance_corrections, [teff, logg, feh, xi, abund], method="linear")

                if starn == "Sun" and abund == 0.0:
                    print("nlte_impact and line", nlte_impact, line)

                interp_ew = interpn([tg, gg, fg, xg, log_eps], ew_corrections, [teff, logg, feh, xi, abund], method="linear")
                reduced_ew  = (np.log10(((10 ** interp_ew) / 1000) / line))
                if line in redew_dict[starn]:
                    redew_dict[starn][line].append(reduced_ew)
                else:
                    redew_dict[starn][line] = []
                    redew_dict[starn][line].append(reduced_ew)

                if reduced_ew > -4.9:# or reduced_ew < -6:
                    continue
                if np.isnan(nlte_impact):
                    continue

                if starn not in full_galah_nlte:
                    full_galah_nlte[starn] = {}
                # Turn from Ti/Fe to A(Ti) and then later to Ti/H (which is what we want to use to calc mean, not ti/fe. We turn to ti/fe in another script when applying it)
                lte_abundance = round(abund, 4)
                nlte_abundance = round(nlte_impact[0] + abund, 4)

                # We input actual A(ti), and then remove the A(ti) of the solar twins in the other function to get Ti/H and then remove feh to get Ti/fe later in another script

                full_galah_nlte[starn][line] = [teff, logg, feh, xi, 0, lte_abundance ,nlte_abundance, round(nlte_abundance-lte_abundance, 4)]
                linecount+=1
                lines_mean.append(line)
            # If the key is never made it means we don't have any abund at that star for whatever reason
            # (likely interpolation range), or if we have too few lines in the star (under 3?)
            if linecount == 0:
                print("too few lines")
                continue

            # mean is the nlte impact
            dev = np.std([i[-1] for i in full_galah_nlte[starn].values()])
            mean = np.average([i[-1] for i in full_galah_nlte[starn].values()])
            # ltemean is the ACTUAL lte abundance for the star
            ltemean = np.average([i[-3] for i in full_galah_nlte[starn].values()])
            star_mean_dict[starn][abund] = [teff, logg, feh, xi, ltemean, dev, mean]
            feh_list.append(feh)

            nlte_list.append(mean)

    print(star_mean_dict)
    shifter = 1
    for starn in star_mean_dict:
        shifter+=1
        enhance = star_mean_dict[starn].keys()

        """if shifter %2 == 0:
            enhance = star_mean_dict[starn].keys()
        else:
            enhance = np.asarray(list(star_mean_dict[starn].keys()))

        if starn == 122563:
            enhance = np.asarray(list(star_mean_dict[starn].keys())) - 0.015
        if starn == 84937:
            enhance = np.asarray(list(star_mean_dict[starn].keys())) + 0.015"""
        nlte_impact = [x[-1] for x in star_mean_dict[starn].values()]
        dev  = [x[-2] for x in star_mean_dict[starn].values()]
        teff = star_mean_dict[starn][-1.0][0]
        logg = star_mean_dict[starn][-1.0][1]
        feh = star_mean_dict[starn][-1.0][2]
        vmic = star_mean_dict[starn][-1.0][3]
        #plt.plot(star_mean_dict[starn].keys(), nlte_impact, label=str(starn))
        #eb = plt.errorbar(enhance, nlte_impact, dev, label=str(starn), markersize=5, capsize=5, marker="o")
        plt.plot(enhance, nlte_impact, label=str(starn), markersize=5, marker="o")
        #eb[-1][0].set_linestyle('--')
        print("STARN", starn, nlte_impact, dev)

    plt.legend()
    plt.xlabel("[Ti/Fe]$_{\mathrm{LTE}}$")
    if "2" in element:
        plt.ylabel("$\Delta$ A(Ti)$_{\mathrm{Ti\, II, \, NLTE}}$")
    else:
        plt.ylabel("$\Delta$ A(Ti)$_{\mathrm{Ti\, I, \, NLTE}}$")
    fig = plt.gcf()
    if "1" in element:

        plt.xlim(-1.05,0.8)
        plt.ylim(-0.05, 0.4)
    else:
        plt.xlim(-1.05,0.8)
        plt.ylim(-0.05, 0.4)

    """fig.set_figwidth(8)
    fig.set_figheight(7)"""
    plt.savefig("graphs/"+str(element)+"nlte_abundancedependence.png", bbox_inches="tight")

    plt.show()
specific_stars(element)

def grid_points():
    import pickle as pkl
    import numpy as np
    import matplotlib.pyplot as plt
    element = "Ti1"
    abund_dict = pkl.load(open(element + "nlte_impact_dict.pkl", "rb"))
    ew_dict = pkl.load(open(element + "ew_dict.pkl", "rb"))

    stars = [[5000.0, 4.5, 0, 1], [4000.0, 2.0, -0.5, 1], [6500.0, 4.0, -2.0, 1], [6000.0, 3.5, -2.5, 1],
             [4500.0, 1.5, -2.5, 2]]
    for teff in abund_dict:
        for logg in abund_dict[teff]:
            for feh in abund_dict[teff][logg]:

                for vmic in abund_dict[teff][logg][feh]:

                    if [teff, logg, feh, vmic] not in stars:
                        continue
                    print(teff, logg, feh, vmic)
                    impact_plot = []
                    plotab = []

                    for eps in abund_dict[teff][logg][feh][vmic]:
                        epslist = [i for i in (abund_dict[teff][logg][feh][vmic].keys())]
                        abslist = epslist - feh - 4.9
                        impactlist = []
                        for line in abund_dict[teff][logg][feh][vmic][eps]:
                            impact = abund_dict[teff][logg][feh][vmic][eps][line]
                            ew = ew_dict[teff][logg][feh][vmic][eps][line]

                            red_ew = (np.log10(((10 ** ew) / 1000) / line))
                            if feh == -0.5:
                                print("here", eps, red_ew)
                            if red_ew > -4.9:
                                continue

                            impactlist.append(impact)
                        plotab.append(eps - feh - 4.9)

                        meanimpact = np.nanmean(impactlist)
                        impact_plot.append(meanimpact)
                    print(plotab, impact_plot)
                    plt.plot(plotab, impact_plot, label=(", ").join((str(teff), str(logg), str(feh), str(vmic))))
    plt.legend()
    plt.xlabel("$\epsilon$")
    plt.ylabel("$\Delta$ A(Ti)$_{\mathrm{Ti\, I}}$")
    plt.show()
grid_points()
