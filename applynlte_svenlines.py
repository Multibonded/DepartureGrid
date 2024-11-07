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
full_galah_nlte = {}

"""
Star    Teff    Logg    Fe/H    Mic
Sun     5772    4.44    0       0.9
Arcturus 4286   1,64    -0.53   1.3
84937   6356    4.06    -2.06   1.2
140283  5792    3.65    -2.36   1.3
122563  4646    1.4     -2,5    1.8

"""

def solar_twin_lines(element):
    b = pkl.load(open(element + "_Galah_impacts2_Ati.pkl", "rb"))
    c = pkl.load(open(element + "_Galah_impacts_lineless2_Ati.pkl", "rb"))
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
        print(11,starn)

        for line in b[starn] :
            if line == 4801:
                mean[line] = [4.9, 4.9]
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
        else:
            mean[line] = [4.9, 4.9]
    return mean

def run_all(element):
    star_mean_dict = {}

    counter = 1
    feh_list = []
    nlte_list = []
    star_count = pkl.load(open(element+"Accurate_param_indexes_galah.pkl", "rb"))
    abund_dict = pkl.load(open(element + "nlte_impact_dict.pkl", "rb"))
    # This one we only want star info, including mean abundance nlte impact, but not line info as that clutters
    # it up . For that we use the full Galah dict.
    # Finds the mean A(Ti) of titanium for solar twins so we can more accurately find [Ti/H] which is better than Ti/fe
    # for reasons. Given as dict[line] = [lte, nlte]
    if "1" in element:
        lines = ["4758", "4759", "4778", "4781", "4797", "4801", "4820", "5689", "5716", "5720", "5739", "5866"]
        lines = ["4758", "4759", "4781", "4797", "4820", "5689", "5716", "5720", "5739"]
        lines = ["4758", "4759", "4781", "4820", "5739"]
        ti_ab = {}

        for line_ab in lines:
            if line_ab == "4781":
                ti_ab[int(line_ab)] = atom1["ind_Ti4782_fe"]
            elif line_ab == "4797":
                ti_ab[int(line_ab)] = atom1["ind_Ti4798_fe"]
            elif line_ab == "4801":
                ti_ab[int(line_ab)] = atom1["ind_Ti4802_fe"]
            elif line_ab == "6716":
                ti_ab[int(line_ab)] = atom1["ind_Ti6717_fe"]
            elif line_ab == "7852":
                ti_ab[int(line_ab)] = atom1["ind_Ti7853_fe"]
            else:
                ti_ab[int(line_ab)] = atom1['ind_Ti' + str(line_ab) + "_fe"]

    elif "2" in element:
        lines = ["4874", "4865", "4849", "4798", "4764", "4719"]
        lines = ["4874", "4798", "4764"]
        lines = ["4719", "4798", "4764", "4865"]
        ti_ab = {}

        for line_ab in lines:
            if line_ab == "4865":
                ti_ab[int(line_ab)] = atom1["ind_Ti4866_fe"]
            elif line_ab == "4798":
                ti_ab[int(line_ab)] = atom1["ind_Ti4799_fe"]
            elif line_ab == "4764":
                ti_ab[int(line_ab)] = atom1["ind_Ti4765_fe"]
            elif line_ab == "4719":
                ti_ab[int(line_ab)] = atom1["ind_Ti4720_fe"]
            else:
                ti_ab[int(line_ab)] = atom1['ind_Ti' + str(line_ab) + "_fe"]



    for starn in star_count:
        counter += 1
        if counter % 500 == 0:
            print(counter, "/", len(star_count))
        # number of lines (after we restrict them) per star, remove is less than 3
        linecount = 0
        lines_mean = []
        for line_index in range(len(lines)):
            line = int(lines[line_index])

            # first 2 are sneden, 4778 is sus
            if line == 5193 or line == 4533 or line == 4778:
                continue
            teff = (atom1[starn]['teff'])
            logg = (atom1[starn]['logg'])
            xi = (atom1[starn]['vmic'])
            feh = (atom1[starn]['fe_h'])
            # check line flags
            try:
                flag = atom1[starn]['ind_flag_Ti' + str(line)]
            except KeyError:
                flag = atom1[starn]['ind_flag_Ti' + str(line + 1)]
            if flag == 1:
                continue

            # paramS 6324.3784 3.9686449 1.4542778 -0.38920593

            pklfile = open('grids/{}.pkl'.format(element), 'rb')
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
            galah_abund_value = ti_ab[line][starn]#+feh+4.9
            if not (min(tg) <= teff <= max(tg)) or not (min(gg) <= logg <= max(gg)) or not (
                    min(fg) <= feh <= max(fg)) or not (min(xg) <= xi <= max(xg)) or not (
                    min(log_eps) <= galah_abund_value  <= max(log_eps)):
                if line == 5866:
                    print("                 Teff,      logg,     feh,       xi,      enh")
                    print("Not in range:", teff, logg, feh, xi, ti_ab[line][starn], starn)
                continue

            ti = np.where(tg > teff)
            ti = min(ti[0]) - 1
            gi = np.where(gg > logg)
            gi = min(gi[0]) - 1
            fi = np.where(fg > feh)
            fi = min(fi[0]) - 1
            log_eps_i = np.where(log_eps > galah_abund_value)

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
                                #print(line, t, g, f, x, a, [abund_dict[tg[t]][gg[g]][fg[f]][xg[x]][log_eps[a]][line]][0])
                                """if starn == 59329:
                                    print("\nLine      Teff       Logg      Feh       Vmic       Enhance      NLTE effect")
                                    print(line, tg[t], "    ", gg[g], "    ", fg[f], "      ", xg[x], "    ",
                                          log_eps[a], "            ", round(abund_dict[tg[t]][gg[g]][fg[f]][xg[x]][log_eps[a]][line], 4))
                                """
            # Log eps is still just [Ti/Fe] at this point, change it to A(Ti) with +fg[f] + 4.9
            nlte_impact = interpn([tg, gg, fg, xg, np.round(log_eps+fg+4.9, 4)], abundance_corrections, [teff, logg, feh, xi, ti_ab[line][starn]+feh+4.9], method="linear")

            if starn not in full_galah_nlte:
                full_galah_nlte[starn] = {}
            # Turn from Ti/Fe to A(Ti) and then later to Ti/H (which is what we want to use to calc mean, not ti/fe. We turn to ti/fe in another script when applying it)
            lte_abundance = round(ti_ab[line][starn]+4.9+feh, 4)
            nlte_abundance = round(nlte_impact[0] + ti_ab[line][starn]+feh+4.9, 4)

            if counter % 200 == 0:
                print(counter, "/", len(star_count))

                print("\nLine and index:", line, starn, "\n Teff, Logg, Feh, Xi, Ab:",
                      [teff, logg, feh, xi, ti_ab[line][starn]])
                print("LTE Abund      NLTE Abund       Nlte Imbalance A(Ti)")
                print(round(ti_ab[line][starn], 4), "         ", round(nlte_impact[0] + ti_ab[line][starn], 4),
                      "          ", round(nlte_impact[0], 4))
                print("LTE Abund      NLTE Abund       Nlte Imbalance [Ti/H]")

                print(lte_abundance, "         ", nlte_abundance, "         ", nlte_impact)
                print("LTE Abund      NLTE Abund       Nlte Imbalance [X/Fe]")
                print(lte_abundance-feh, "         ", nlte_abundance-feh, "         ", nlte_impact)

            # We input actual A(ti), and then remove the A(ti) of the solar twins in the other function to get Ti/H and then remove feh to get Ti/fe later in another script
            full_galah_nlte[starn][line] = [teff, logg, feh, xi, lte_abundance ,nlte_abundance, round(nlte_abundance-lte_abundance, 4)]
            linecount+=1
            lines_mean.append(line)
        # If the key is never made it means we don't have any abund at that star for whatever reason
        # (likely interpolation range), or if we have too few lines in the star (under 3?)
        if linecount == 0:
            continue
        elif linecount < 2:
            full_galah_nlte.pop(starn)
            continue

        try:
            # mean is the nlte impact
            mean = np.nanmean([i[-1] for i in full_galah_nlte[starn].values()])
            # ltemean is the ACTUAL lte abundance for the star
            ltemean = np.nanmean([i[-3] for i in full_galah_nlte[starn].values()])
            if np.isnan(mean):
                full_galah_nlte.pop(starn)
                continue
            star_mean_dict[starn] = [teff, logg, feh, xi, ltemean, mean]
            feh_list.append(feh)

            nlte_list.append(mean)
        except KeyError:
            pass

    pickle.dump(full_galah_nlte, open(element + "_Galah_impacts2_Ati.pkl", "wb"))
    pickle.dump(star_mean_dict, open(element + "_Galah_impacts_lineless2_Ati.pkl", "wb"))

# Turn the A(ti) into Ti/H by using the new solar mean of A(Ti)_sun for each line
def turn_ti_tih(element):

    star_mean_dict = {}
    solar_values = solar_twin_lines(element)
    feh_list = []
    nlte_list = []
    print(full_galah_nlte)
    for starn in full_galah_nlte:
        for line in full_galah_nlte[starn]:
            # Same value in each line. Could just do it once but eh
            teff = full_galah_nlte[starn][line][0]
            logg = full_galah_nlte[starn][line][1]
            feh = full_galah_nlte[starn][line][2]
            xi = full_galah_nlte[starn][line][3]

            solar_line = solar_values[line]
            print("ya", line, solar_line)
            lte_abundace = full_galah_nlte[starn][line][-3]-solar_line[0]
            nlte_abundace = full_galah_nlte[starn][line][-2]-solar_line[1]
            nlte_difference = round(nlte_abundace-lte_abundace, 4)

            full_galah_nlte[starn][line] = [teff, logg, feh, xi, lte_abundace, nlte_abundace, nlte_difference]

        # mean is the nlte impact NOT mean nlte abundance

        mean = np.nanmean([i[-1] for i in full_galah_nlte[starn].values()])




        # ltemean is the ACTUAL lte abundance for the star
        ltemean = np.nanmean([i[-3] for i in full_galah_nlte[starn].values()])
        if np.isnan(mean):
            full_galah_nlte.pop(starn)
            continue
        star_mean_dict[starn] = [teff, logg, feh, xi, ltemean, mean]
        feh_list.append(feh)

        nlte_list.append(mean)

    #full_galah_nlte[starn]['average'] = np.mean(full_galah_nlte[starn])
    save = True
    if save:

        pickle.dump(full_galah_nlte, open(element+"_Galah_impacts2_svenlines.pkl", "wb"))
        pickle.dump(star_mean_dict, open(element+"_Galah_impacts_lineless2_svenlines.pkl", "wb"))

        pickle.dump([feh_list, nlte_list], open(element+"_NLTE_list2_svenlines.pkl", "wb"))

    print("DONE")
    plt.scatter(feh_list, nlte_list, s=0.5, c="b")

    plt.show()





element = "Ti1"


reliable     = pkl.load(open(element+"Accurate_param_indexes_galah.pkl", "rb"))

# Run once to make a dict of nlte and lte values for each star and line in A(Ti) form
run_all(element)
# Then turn to Ti/H by taking the lte/nlte mean of solar twinms. We later turn it into Ti/Fe in other scripts
turn_ti_tih(element)

