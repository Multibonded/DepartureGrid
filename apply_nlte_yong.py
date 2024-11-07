"""Applies NLTE effects to GALAH stars via interpolation. Produces a few files that has the info for all
galah stars and their nlte impact """
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
full_yong_nlte = {}



"""
Star    Teff    Logg    Fe/H    Mic
Sun     5772    4.44    0       0.9
Arcturus 4286   1,64    -0.53   1.3
84937   6356    4.06    -2.06   1.2
140283  5792    3.65    -2.36   1.3
122563  4646    1.4     -2,5    1.8

"""

def solar_twin_lines(element):
    if "1" in element:
        lines =['3354' ,'3653', '3741', '4533', '5193']

    elif "2" in element:
        lines = ["4874", "4865", "4849", "4798", "4764", "4719"]
        lines = ["4719", "4798", "4764", "4866"]
        lines = ["4798", "4719", "4874"]

    b = pkl.load(open(element + "_yong_impacts2_Ati.pkl", "rb"))
    c = pkl.load(open(element + "_yong_impacts_lineless2_Ati.pkl", "rb"))
    feh = np.asarray([i[2] for i in c.values()])
    logg = np.asarray([i[1] for i in c.values()])
    teff = np.asarray([i[0] for i in c.values()])
    print("teff", teff)
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
    errors = {}
    for starn in index_list:
        print(11,starn)

        for line in b[starn] :
            if line == 4801:
                mean[line] = [4.9, 4.9]
            if line in lte:
                lte[line].append(b[starn][line][-3])
                nlte[line].append(b[starn][line][-2])
                errors[line].append(b[starn][line][-4])
            else:
                lte[line]=[b[starn][line][-3]]
                nlte[line] = [b[starn][line][-2]]
                errors[line] = [b[starn][line][-4]]

    for line in lines:
        try:
            if len(lte[line]) > 5:
                mean[line] = [np.average(lte[line], weights=errors[line]), np.average(nlte[line], weights=errors[line])]
        except KeyError:
            mean[int(line)] = [4.9, 4.9]
    print(mean)
    return mean

def run_all(element):
    star_mean_dict = {}
    atom1 = pickle.load(open(element + "yong_abundances.pkl", "rb"))
    counter = 1
    feh_list = []
    nlte_list = []
    star_count = atom1.keys()
    abund_dict = pkl.load(open(element + "nlte_impact_dict_yong.pkl", "rb"))
    # This one we only want star info, including mean abundance nlte impact, but not line info as that clutters
    # it up . For that we use the full Galah dict.
    # Finds the mean A(Ti) of titanium for solar twins so we can more accurately find [Ti/H] which is better than Ti/fe
    # for reasons. Given as dict[line] = [lte, nlte]
    if "1" in element:
        lines = ["4533", "5193"]

    elif "2" in element:
        lines = ["4443", "4468", "4501", "4571"]

    for starn in star_count:
        counter += 1
        if counter % 500 == 0:
            print(counter, "/", len(star_count))
        # number of lines (after we restrict them) per star, remove is less than 3
        linecount = 0
        lines_mean = []
        errors = []
        for line_index in range(len(lines)):
            line = int(lines[line_index])

            # first 2 are yong, 4778 is sus
            try:
                teff = float(atom1[starn][line][0])
            except KeyError:
                try:
                    print("plus 1")
                    line = line+1
                    teff = float(atom1[starn][line][0])
                except KeyError:
                    print("No line", line, "in star", starn)
                    continue

            logg = float(atom1[starn][line][1])
            xi = float(atom1[starn][line][3])
            feh = float(atom1[starn][line][2])
            pklfile = open('grids/{}_yong.pkl'.format(element), 'rb')
            grid = pkl.load(pklfile)
            tg = grid["tg"]

            gg = grid["gg"]
            xg = np.array(grid["xg"])
            fg = np.array(grid["fg"])
            # This does not affect ionisation imbalance just definitive abundance predictions in nlte_interp. Unsure how
            # but as we're applying dep. coff. ontothe abund of GALAH I don't hink it matters here for us.
            #log_eps = grid['ag']
            #log_eps = [round(x, 4) for x in log_eps]

            log_eps = grid["ag"]# + 4.9 + feh

            log_eps = [round(x, 4) for x in log_eps]
            galah_abund_value_eps = atom1[starn][line][-1]#+feh+4.9
            galah_abund_value=galah_abund_value_eps-4.9-feh
            if not (min(tg) <= teff <= max(tg)) or not (min(gg) <= logg <= max(gg)) or not (
                    min(fg) <= feh <= max(fg)) or not (min(xg) <= xi <= max(xg)) or not (
                    min(log_eps) <= galah_abund_value  <= max(log_eps)):
                print(line, galah_abund_value, feh)
                continue

            ti = np.where(tg > teff)
            ti = min(ti[0]) - 1
            gi = np.where(gg > logg)
            gi = min(gi[0]) - 1
            fi = np.where(fg > feh)
            fi = min(fi[0]) - 1
            log_eps_i = np.where(np.asarray(log_eps) > galah_abund_value)

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
            # MAkes a np array with all parameters in the correct slots to allow easy itnterpolation as the nlte corrections
            # are int he sane shape as the parameter list.
            for t in range(len(tg)):
                for g in range(len(gg)):
                    for f in range(len(fg)):
                        for x in range(len(xg)):
                            for a in range(len(log_eps)):
                                # Abund dict has the NLTE impact. Log eps is still just [Ti/Fe] at this point, change it to A(Ti) with +fg[f] + 4.9 (Star metalcitiy and solar abundance)
                                abundance_corrections[t, g, f, x, a] = ([abund_dict[tg[t]][gg[g]][fg[f]][xg[x]][np.round(log_eps[a]+fg[f]+4.9, 4)][line]][0])
                                """if starn == "1099549853":
                                    print(line, abund_dict[tg[t]][gg[g]][fg[f]][xg[x]][np.round(log_eps[a]+fg[f]+4.9, 4)][line], tg[t], gg[g], fg[f], xg[x], np.round(log_eps[a], 4))"""

            # yong starts in A(Ti)
            nlte_impact = interpn([tg, gg, fg, xg, np.round(log_eps, 4)], abundance_corrections, [teff, logg, feh, xi, galah_abund_value], method="linear")
            galah_abund_value = galah_abund_value+feh+4.9
            """if starn == "1099549853":
                print(nlte_impact)
                print(galah_abund_value)
                exit()"""

            if np.isnan(nlte_impact):
                continue

            if starn not in full_yong_nlte:
                full_yong_nlte[starn] = {}

            # yong comes with abundances in A(Ti) so no need for changing from Ti/fe
            lte_abundance = round(galah_abund_value, 4)
            nlte_abundance = round(nlte_impact[0] + galah_abund_value, 4)

            if counter % 200 == 0:
                print(counter, "/", len(star_count))

                print("\nLine and index:", line, starn, "\n Teff, Logg, Feh, Xi, Ab:",
                      [teff, logg, feh, xi, galah_abund_value])
                print("LTE Abund      NLTE Abund       Nlte Imbalance A(Ti)")
                print(round(galah_abund_value, 4), "         ", round(nlte_impact[0] + galah_abund_value, 4),
                      "          ", round(nlte_impact[0], 4))
                print("LTE Abund      NLTE Abund       Nlte Imbalance [Ti/H]")

                print(lte_abundance, "         ", nlte_abundance, "         ", nlte_impact)
                print("LTE Abund      NLTE Abund       Nlte Imbalance [X/Fe]")
                print(lte_abundance-feh, "         ", nlte_abundance-feh, "         ", nlte_impact)

            # We input actual A(ti), and then remove the A(ti) of the solar twins in the other function to get Ti/H and then remove feh to get Ti/fe later in another script

            try:
                error = 1
            except KeyError:
                error = 1
            full_yong_nlte[starn][line] = [teff, logg, feh, xi, error, lte_abundance ,nlte_abundance, round(nlte_abundance-lte_abundance, 4)]
            linecount+=1
            lines_mean.append(line)
        # If the key is never made it means we don't have any abund at that star for whatever reason
        # (likely interpolation range), or if we have too few lines in the star (under 3?)
        if linecount == 0:
            continue
        elif linecount < 2:
            full_yong_nlte.pop(starn)
            continue
        try:
            # mean is the nlte impact
            errors = [i[-4] for i in full_yong_nlte[starn].values()]
            mean = np.average([i[-1] for i in full_yong_nlte[starn].values()], weights=errors)
            # ltemean is the ACTUAL lte abundance for the star
            ltemean = np.average([i[-3] for i in full_yong_nlte[starn].values()], weights = errors)
            if np.isnan(mean):
                full_yong_nlte.pop(starn)
                continue
            star_mean_dict[starn] = [teff, logg, feh, xi, ltemean, mean]
            feh_list.append(feh)

            nlte_list.append(mean)
        except KeyError:
            pass

    pickle.dump(full_yong_nlte, open(element + "_yong_impacts2_Ati.pkl", "wb"))
    pickle.dump(star_mean_dict, open(element + "_yong_impacts_lineless2_Ati.pkl", "wb"))

# Turn the A(ti) into Ti/H by using the new solar mean of A(Ti)_sun for each line
def turn_ti_tih(element):

    star_mean_dict = {}
    solar_values = [4.9, 4.9]
    print("Solar", solar_values)
    feh_list = []
    nlte_list = []
    for starn in full_yong_nlte:
        for line in full_yong_nlte[starn]:
            # Same value in each line. Could just do it once but eh
            teff = full_yong_nlte[starn][line][0]
            logg = full_yong_nlte[starn][line][1]
            feh = full_yong_nlte[starn][line][2]
            xi = full_yong_nlte[starn][line][3]

            solar_line = [4.9, 4.9]
            lte_abundace = full_yong_nlte[starn][line][-3]-solar_line[0]
            nlte_abundace = full_yong_nlte[starn][line][-2]-solar_line[1]
            nlte_difference = round(nlte_abundace-lte_abundace, 4)

            full_yong_nlte[starn][line] = [teff, logg, feh, xi, lte_abundace, nlte_abundace, nlte_difference]

        # mean is the nlte impact NOT mean nlte abundance

        mean = np.nanmean([i[-1] for i in full_yong_nlte[starn].values()])




        # ltemean is the ACTUAL lte abundance for the star
        ltemean = np.nanmean([i[-3] for i in full_yong_nlte[starn].values()])
        if np.isnan(mean):
            full_yong_nlte.pop(starn)
            continue
        star_mean_dict[starn] = [teff, logg, feh, xi, ltemean, mean]
        feh_list.append(feh)

        nlte_list.append(mean)

    #full_yong_nlte[starn]['average'] = np.mean(full_yong_nlte[starn])
    save = True
    if save:

        pickle.dump(full_yong_nlte, open(element+"_yong_impacts2.pkl", "wb"))
        pickle.dump(star_mean_dict, open(element+"_yong_impacts_lineless2.pkl", "wb"))

        pickle.dump([feh_list, nlte_list], open(element+"_yong_NLTE_list2.pkl", "wb"))

    print("DONE")
    plt.scatter(feh_list, nlte_list, s=0.5, c="b")

    plt.show()





element = "Ti1"



# Run once to make a dict of nlte and lte values for each star and line in A(Ti) form
run_all(element)
# Then turn to Ti/H by taking the lte/nlte mean of solar twinms. We later turn it into Ti/Fe in other scripts
turn_ti_tih(element)

