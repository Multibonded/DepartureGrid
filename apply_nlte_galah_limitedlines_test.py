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
startolook = 597618
print(atom1['ind_Ti4799_fe'][startolook], atom1['fe_h'][startolook])
print(atom1['ind_Ti4799_fe'][startolook] + atom1['fe_h'][startolook] + 4.9)
full_galah_nlte = {}
all = True
"""
Star    Teff    Logg    Fe/H    Mic
Sun     5772    4.44    0       0.9
Arcturus 4286   1,64    -0.53   1.3
84937   6356    4.06    -2.06   1.2
140283  5792    3.65    -2.36   1.3
122563  4646    1.4     -2,5    1.8

"""
reduced = True
def solar_twin_lines(element, all, linekeep):

    if reduced:
        b = pkl.load(open(element + "_Galah_impacts2_Ati_reduced.pkl", "rb"))
        c = pkl.load(open(element + "_Galah_impacts_lineless2_Ati_reduced.pkl", "rb"))
        if all:
            b = pkl.load(open(element + "_Galah_impacts2_Ati_reduced_all.pkl", "rb"))
            c = pkl.load(open(element + "_Galah_impacts_lineless2_Ati_reduced_all.pkl", "rb"))

    else:
        b = pkl.load(open(element + "_Galah_impacts2_Ati.pkl", "rb"))
        c = pkl.load(open(element + "_Galah_impacts_lineless2_Ati.pkl", "rb"))
        if all:
            b = pkl.load(open(element + "_Galah_impacts2_Ati_all.pkl", "rb"))
            c = pkl.load(open(element + "_Galah_impacts_lineless2_Ati_all.pkl", "rb"))
    print(len(b), len(c))
    feh = np.asarray([i[2] for i in c.values()])
    logg = np.asarray([i[1] for i in c.values()])
    teff = np.asarray([i[0] for i in c.values()])
    vmic = np.asarray([i[3] for i in c.values()])
    r_est_dr2 = np.asarray([i[5] for i in c.values()])
    snr_c2_iraf = np.asarray([i[4] for i in c.values()])
    # Indexes of stars close to our sun's parameters
    twindex = np.where(
        np.logical_and(
        np.logical_and(
        np.logical_and(
            abs(teff-5772) <= 60,
            abs(logg-4.438) <= 0.05),
            abs(vmic-1) <= 0.2),
            abs(feh) <= 0.05)

    )[0]

    index_list = np.asarray(list(b.keys()))[twindex]


    # we have indexes of solar twins,n ot the star number which is whbat we use in the dictioanry so now we convert
    keylist = []
    # lists to contain the lkte and nlte titanium abundances (lte from galah, nlte from us) for solar twins for us to
    # use line by line
    lte = {}
    nlte = {}
    mean = {}
    errors = {}
    ltesquare = {}
    for starn in index_list:
        teffstar = atom1[starn]['teff']
        fehstar = atom1['fe_h'][starn]
        gstar = atom1['logg'][starn]
        for line in b[starn] :
            if int(line) not in linekeep:
                continue
            #print("b", b[starn])
            if line in lte:
                lte[line].append(b[starn][line][-3])
                nlte[line].append(b[starn][line][-2])
                errors[line].append(b[starn][line][-4])
                if line == "5739_single":
                    ltesquare[line].append(atom1[starn]['ind_Ti5739_fe'])
                else:
                    key="ind_Ti"+str(line)+"_fe"
                    key2 = "ind_Ti"+str(int(line)+1)+"_fe"
                    try:
                        ltesquare[line].append(atom1[starn][key])
                    except KeyError:
                        ltesquare[line].append(atom1[starn][key2])

            else:
                lte[line]=[b[starn][line][-3]]
                nlte[line] = [b[starn][line][-2]]
                errors[line] = [b[starn][line][-4]]
                if line == "5739_single":
                    ltesquare[line]=([atom1[starn]['ind_Ti5739_fe']])
                else:
                    try :
                        ltesquare[line]=[atom1[starn]['ind_Ti' + str(line) + '_fe']]
                    except KeyError:
                        ltesquare[line]=[atom1[starn]['ind_Ti' + str(int(line) + 1) + '_fe']]



    total = []
    totalnlte = []
    totalerrors = []
    lenlist = []
    meanerrors = []
    for line in lte:
        total.extend(lte[line])

        print("length of ", line, ":", len(lte[line]))
        lenlist.append(len(lte[line]))
        totalnlte.extend(nlte[line])
        totalerrors.extend(errors[line])
        meanerrors.append(np.mean(errors[line]))
        print("std", line, np.std(lte[line]))

        if len(lte[line]) > 2:
            mean[line] = [np.average(lte[line], weights=errors[line]), np.average(nlte[line], weights=errors[line])]
            meangalah = np.average(ltesquare[line], weights=errors[line])
            print("square", line, meangalah)

        else:
            print("uhoh")
            exit()
        #print(line, mean[line])
    print(total)
    print("\nLTE and NLTE averaged from EACH line in EACH star, weighted as usual (But remember there are many more 4758/4759 lines and are thus biased)"
          , np.average(np.asarray(total), weights=np.asarray(totalerrors)), np.average(np.asarray(totalnlte), weights=np.asarray(totalerrors)))
    print("NLTE impact is thus: ", np.average(np.asarray(totalnlte), weights=np.asarray(totalerrors))- np.average(np.asarray(total), weights=np.asarray(totalerrors)))
    print("unweighted total avrg", np.mean(total), np.mean(totalnlte))
    print("error total:", np.std(total), np.std(totalnlte))
    lte = (np.average([i[0] for i in mean.values()], weights=meanerrors))
    nlte =(np.average([i[1] for i in mean.values()], weights=meanerrors))
    print("\n Following are lte and nlte values for averaged of each line, weighted by each line's mean error:")
    print("LTE:", lte, "NLTE:", nlte)
    print("LTE errors:", np.std([i[0] for i in mean.values()]))
    print("NLTE errors:", np.std([i[1] for i in mean.values()]))
    print("AVerage based on number of lines:", np.average(lte))
    corrections = (nlte-lte)
    print("Corrections", corrections)
    print("mean", mean)
    return mean, corrections, index_list

def run_all(element):
    star_mean_dict = {}
    redew_dict = {}
    counter = 1
    feh_list = []
    nlte_list = []
    """if reduced:
        star_count = pkl.load(open(element+"Accurate_param_indexes_galah.pkl", "rb"))
        abund_dict = pkl.load(open(element + "nlte_impact_dict_reduced.pkl", "rb"))
    else:
        abund_dict = pkl.load(open(element + "nlte_impact_dict.pkl", "rb"))
        star_count = pkl.load(open(element + "Accurate_param_indexes_galah.pkl", "rb"))
    """
    if all:
        pass
    abund_dict = pkl.load(open(element + "nlte_impact_dict_all.pkl", "rb"))
    ew_dict = pkl.load(open(element + "ew_dict_all.pkl", "rb"))
    star_count = pkl.load(open(element + "Accurate_param_indexes_galah.pkl", "rb"))
    """else:
        abund_dict = pkl.load(open(element + "nlte_impact_dict.pkl", "rb"))
        ew_dict = pkl.load(open(element + "ew_dict.pkl", "rb"))
        star_count = pkl.load(open(element + "Accurate_param_indexes_galah.pkl", "rb"))"""

    # This one we only want star info, including mean abundance nlte impact, but not line info as that clutters
    # it up . For that we use the full Galah dict.
    # Finds the mean A(Ti) of titanium for solar twins so we can more accurately find [Ti/H] which is better than Ti/fe
    # for reasons. Given as dict[line] = [lte, nlte]
    if "1" in element:
        if all:
            lines = ["4758", "4759", "4778", "4781", "4797", "4801", "4820", "5689", "5716", "5720", "5739", "5866", "6716", "7852"]
        else:
            lines = ["4758", "4759", "5689", "5716", "5739", "6599"]
            lines = ["4758", "4759", "5689", "5739"]
            #lines = ['4758','4759','4781','4801','4820','5739'] # galah ones
        ti_ab = {}

        for line_ab in lines:
            if line_ab == "4781":
                ti_ab[(line_ab)] = atom1["ind_Ti4782_fe"]
            elif line_ab == "4797":
                ti_ab[(line_ab)] = atom1["ind_Ti4798_fe"]
            elif line_ab == "5739_single":
                ti_ab[(line_ab)] = atom1["ind_Ti5739_fe"]
            elif line_ab == "4801":
                ti_ab[(line_ab)] = atom1["ind_Ti4802_fe"]
            elif line_ab == "6716":
                ti_ab[(line_ab)] = atom1["ind_Ti6717_fe"]
            elif line_ab == "7852":
                ti_ab[(line_ab)] = atom1["ind_Ti7853_fe"]
            else:
                ti_ab[(line_ab)] = atom1['ind_Ti' + str(line_ab) + "_fe"]

    elif "2" in element:
        if all:
            lines = ["4874", "4865", "4849", "4798", "4764", "4719"]

        else:
            lines = ["4798", "4719", "4874"]
        ti_ab = {}

        for line_ab in lines:
            if line_ab == "4865":
                ti_ab[(line_ab)] = atom1["ind_Ti4866_fe"]
            elif line_ab == "4798":
                ti_ab[(line_ab)] = atom1["ind_Ti4799_fe"]
            elif line_ab == "4764":
                ti_ab[(line_ab)] = atom1["ind_Ti4765_fe"]
            elif line_ab == "4719":
                ti_ab[(line_ab)] = atom1["ind_Ti4720_fe"]
            else:
                ti_ab[(line_ab)] = atom1['ind_Ti' + str(line_ab) + "_fe"]



    for starn in star_count:
        if starn not in redew_dict:
            redew_dict[starn] = {}
        counter += 1
        if counter % 500 == 0:
            print(counter, "/", len(star_count))
        # number of lines (after we restrict them) per star, remove is less than 3
        linecount = 0
        lines_mean = []
        errors = []
        for line_index in range(len(lines)):

            line = (lines[line_index])

            teff = (atom1[starn]['teff'])
            logg = (atom1[starn]['logg'])
            xi = (atom1[starn]['vmic'])
            feh = (atom1[starn]['fe_h'])
            # check line flags
            try:
                flag = atom1[starn]['ind_flag_Ti' + str(line)]
            except KeyError:
                flag = atom1[starn]['ind_flag_Ti' + str(int(line) + 1)]
            if flag == 1:
                continue

            # paramS 6324.3784 3.9686449 1.4542778 -0.38920593
            """if reduced:
                if all:
                    pklfile = open('grids/{}_reduced_all.pkl'.format(element), 'rb')
                else:
                    pklfile = open('grids/{}_reduced.pkl'.format(element), 'rb')
            else:"""
            if all:
                pklfile = open('grids/{}_all.pkl'.format(element), 'rb')
            else:
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
                if counter % 500 == 0:

                    print("                 Teff,      logg,     feh,       xi,      enh, starn, line")
                    print("Not in range:", teff, logg, feh, xi, ti_ab[line][starn], starn, line)

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
            ew_corrections = np.zeros([len(tg), len(gg), len(fg), len(xg), len(log_eps)]) + np.nan
            # MAkes a np array with all parameters in the correct slots to allow easy itnterpolation as the nlte corrections
            # are int he sane shape as the parameter list.
            for t in range(len(tg)):
                for g in range(len(gg)):
                    for f in range(len(fg)):
                        for x in range(len(xg)):
                            for a in range(len(log_eps)):
                                # Abund dict has the NLTE impact. Log eps is still just [Ti/Fe] at this point, change it to A(Ti) with +fg[f] + 4.9 (Star metalcitiy and solar abundance)
                                abundance_corrections[t, g, f, x, a] = ([abund_dict[tg[t]][gg[g]][fg[f]][xg[x]][np.round(log_eps[a]+fg[f]+4.9, 4)][int(line)]][0])
                                ew=(([ew_dict[tg[t]][gg[g]][fg[f]][xg[x]][np.round(log_eps[a]+fg[f]+4.9, 4)][int(line)]][0]))
                                ew_corrections[t, g, f, x, a] = ew
                                #print(line, t, g, f, x, a, [abund_dict[tg[t]][gg[g]][fg[f]][xg[x]][log_eps[a]][line]][0])
                                """if starn == 59329:
                                    print("\nLine      Teff       Logg      Feh       Vmic       Enhance      NLTE effect")
                                    print(line, tg[t], "    ", gg[g], "    ", fg[f], "      ", xg[x], "    ",
                                          log_eps[a], "            ", round(abund_dict[tg[t]][gg[g]][fg[f]][xg[x]][log_eps[a]][line], 4))
                                """
            # Log eps is still just [Ti/Fe] at this point, change it to A(Ti) with +fg[f] + 4.9
            nlte_impact = interpn([tg, gg, fg, xg, np.round(log_eps+fg+4.9, 4)], abundance_corrections, [teff, logg, feh, xi, ti_ab[line][starn]+feh+4.9], method="linear")
            interp_ew = interpn([tg, gg, fg, xg, np.round(log_eps+fg+4.9, 4)], ew_corrections, [teff, logg, feh, xi, ti_ab[line][starn]+feh+4.9], method="linear")
            if lines[line_index] == "5739_single":
                reduced_ew = (np.log10(((10 ** interp_ew) / 1000) / 5739))
            else:
                reduced_ew  = (np.log10(((10 ** interp_ew) / 1000) / int(line)))

            if line in redew_dict[starn]:
                redew_dict[starn][line].append(reduced_ew)
            else:
                redew_dict[starn][line] = []
                redew_dict[starn][line].append(reduced_ew)

            if reduced_ew > -4.9 or reduced_ew < -6.0:
            #if reduced_ew > -4.9:
                #print(reduced_ew, interp_ew, 10**interp_ew,(10 ** interp_ew) / 1000, line, teff, logg, feh, xi, ti_ab[line][starn])
                continue

            if np.isnan(nlte_impact):
                continue

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
            if line == "5739_single":
                error = 1/(atom1['ind_cov_e_Ti5739'][starn]**2)
            else:
                try:
                    error = 1/(atom1['ind_cov_e_Ti' + str(line)][starn]**2)
                except KeyError:
                    error = 1/(atom1['ind_cov_e_Ti' + str(int(line)+1)][starn]**2)
            full_galah_nlte[starn][line] = [teff, logg, feh, xi, atom1['snr_c2_iraf'], atom1['r_est_dr2'], error, lte_abundance ,nlte_abundance, round(nlte_abundance-lte_abundance, 4)]
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
            errors = [i[-4] for i in full_galah_nlte[starn].values()]
            mean = np.average([i[-1] for i in full_galah_nlte[starn].values()], weights=errors)
            # ltemean is the ACTUAL lte abundance for the star
            ltemean = np.average([i[-3] for i in full_galah_nlte[starn].values()], weights = errors)
            if np.isnan(mean):
                full_galah_nlte.pop(starn)
                continue
            star_mean_dict[starn] = [teff, logg, feh, xi, ltemean, mean]
            feh_list.append(feh)

            nlte_list.append(mean)
        except KeyError:
            pass

    if reduced:
        if all:
            pickle.dump(full_galah_nlte, open(element + "_Galah_impacts2_Ati_reduced_all.pkl", "wb"))
            pickle.dump(star_mean_dict, open(element + "_Galah_impacts_lineless2_Ati_reduced_all.pkl", "wb"))
        else:
            pickle.dump(full_galah_nlte, open(element + "_Galah_impacts2_Ati_reduced.pkl", "wb"))
            pickle.dump(star_mean_dict, open(element + "_Galah_impacts_lineless2_Ati_reduced.pkl", "wb"))


    else:
        if all:
            pickle.dump(full_galah_nlte, open(element + "_Galah_impacts2_Ati_all.pkl", "wb"))
            pickle.dump(star_mean_dict, open(element + "_Galah_impacts_lineless2_Ati_all.pkl", "wb"))
        else:
            pickle.dump(full_galah_nlte, open(element + "_Galah_impacts2_Ati.pkl", "wb"))
            pickle.dump(star_mean_dict, open(element + "_Galah_impacts_lineless2_Ati.pkl", "wb"))

    if all:
        pickle.dump(redew_dict, open(element+"Reducedews_all.pkl", "wb"))
    else:
        pickle.dump(redew_dict, open(element+"Reducedews.pkl", "wb"))

# Turn the A(ti) into Ti/H by using the new solar mean of A(Ti)_sun for each line
def turn_ti_tih(element):
    loadall = True
    all = False
    # indexes to remove from fullgalah dict after seeing they're empty or fail a requirement of some knd
    toremove = []
    star_mean_dict = {}
    # additional correction for NLTE impact figures to represent abundance differences in absolute abundance. Not in
    # other figures beyond HR and delta[Ti.H]

    feh_list = []
    nlte_list = []
    if reduced:
        if loadall:
            full_galah_nlte = pickle.load(open(element + "_Galah_impacts2_Ati_reduced_all.pkl", "rb"))
        else:
            full_galah_nlte = pickle.load(open(element + "_Galah_impacts2_Ati_reduced.pkl", "rb"))

    else:
        if loadall:
            full_galah_nlte = pickle.load(open(element + "_Galah_impacts2_Ati_all.pkl", "rb"))
        else:
            full_galah_nlte = pickle.load(open(element + "_Galah_impacts2_Ati.pkl", "rb"))
    
    
    if "1" in element:

        lines = [4758, 4759, 4778, 4781, 4797, 4801, 4820, 5689, 5716, 5720, 5739, 5866, 6716, 7852]

        lines = [4758, 4759, 5689, 5739, 5716]
        lines = [4758, 4759, 4778, 4781, 4797, 4801, 4820, 5689, 5716, 5720, 5739, 5866, 6716, 7852]

        linekeep = [4758, 4759, 4781, 4820, 5739]
        linekeep = [4758, 4759, 5689, 5716, 5739]

        #linekeep = [4758, 4759, 4782, 4820,  5739]
        #linekeep = [4758, 4759, 4781, 4820, 5739] # galah (without 4801) in their zero point
        solarlines = [4.7, 4.72, 5.04, 4.8, 4.82] #galah
        solarlines = [4.7, 4.72, 4.83, 4.82]
        nltecorrection = {4758: 0.0558, 4759: 0.0566, 5689: 0.0504, 5716: 0.04896, 5739: 0.0473}
    else:
        linekeep = [4719, 4764, 4798, 4865]
        linekeep = [4719, 4764, 4798, 4865, 4874]
        linekeep = [4719, 4798, 4874]
        #linekeep = [4719, 4764, 4798, 4865] # galah zero points
        solarlines = [5.12, 4.85, 4.85, 5.12]
        solarlines = [5.12, 4.85, 4.95]
        nltecorrection = {4798: -0.0013, 4719: -0.00155, 4874: -0.0415}

        #linekeep = [4798, 4719, 4874, 4865]
    solar_values, correction, twindex = solar_twin_lines(element, loadall, linekeep)
    for starn in full_galah_nlte:
        firstline = list(full_galah_nlte[starn].keys())[0]
        feh = full_galah_nlte[starn][firstline][2]
        if abs(feh) > 0.05:
            pass
        if starn not in twindex:
            pass
        for line in full_galah_nlte[starn]:
            if int(line) not in linekeep:

                continue
            # Same value in each line. Could just do it once but eh
            teff = full_galah_nlte[starn][line][0]
            logg = full_galah_nlte[starn][line][1]
            feh = full_galah_nlte[starn][line][2]
            xi = full_galah_nlte[starn][line][3]
            lineloc = np.where(int(line) == np.asarray(linekeep))[0][0]


            """if "1" in element:
                solar_line = [4.87, 4.9]
            else:
                solar_line = [4.92, 4.92]"""
            #lte_abundance = full_galah_nlte[starn][line][-3]-solar_values[line][0]
            #nlte_abundance = full_galah_nlte[starn][line][-2]-solar_values[line][1]
            lte_abundance = full_galah_nlte[starn][line][-3]-4.9
            nlte_abundance = full_galah_nlte[starn][line][-2]-4.9-nltecorrection[int(line)]
            print("here", 4.9+correction)
            nlte_difference = round(nlte_abundance-lte_abundance, 4)

            #print("Star and params:", starn, teff, logg, feh, xi)
            #print("NLTE effect and correction:", nlte_difference, correction)
            full_galah_nlte[starn][line] = [teff, logg, feh, xi, lte_abundance, nlte_abundance, nlte_difference]


        nlte_effect = []
        lte_effect = []
        for (line) in list(full_galah_nlte[starn].keys()):
            if int(line) not in linekeep:
                continue
            else:
                nlte_effect.append(full_galah_nlte[starn][line][-1])
                lte_effect.append(full_galah_nlte[starn][line][-3])
        # mean is the nlte impact NOT mean nlte abundance
        if len(nlte_effect) < 2:
            toremove.append(starn)
            continue
        mean = np.nanmean(nlte_effect)

        #mean = np.nanmean([i[-1] for i in full_galah_nlte[starn].values()])




        # ltemean is the ACTUAL lte abundance for the star
        #ltemean = np.nanmean([i[-3] for i in full_galah_nlte[starn].values()])
        ltemean = np.nanmean(lte_effect)
        if np.isnan(mean):
            toremove.append(starn)
            continue
        star_mean_dict[starn] = [teff, logg, feh, xi, ltemean, mean]
        feh_list.append(feh)
        # Removes the solar nlte correction added earlier, as we want to stick to Absolute abundance.
        nlte_list.append(mean+correction)
    for index_remove in toremove:
        full_galah_nlte.pop(index_remove)

    #full_galah_nlte[starn]['average'] = np.mean(full_galah_nlte[starn])
    save = True
    print("lengths (feh, star, full):", len(feh_list), len(star_mean_dict), len(full_galah_nlte))
    if save:

        if reduced:
            if all:
                pickle.dump(full_galah_nlte, open(element + "_Galah_impacts2_reduced_all.pkl", "wb"))
                pickle.dump(star_mean_dict, open(element + "_Galah_impacts_lineless2_reduced_all.pkl", "wb"))

                pickle.dump([feh_list, nlte_list], open(element + "_NLTE_list2_reduced_all.pkl", "wb"))
            else:
                pickle.dump(full_galah_nlte, open(element + "_Galah_impacts2_reduced.pkl", "wb"))
                pickle.dump(star_mean_dict, open(element + "_Galah_impacts_lineless2_reduced.pkl", "wb"))

                pickle.dump([feh_list, nlte_list], open(element + "_NLTE_list2_reduced.pkl", "wb"))

        else:
            if all:
                pickle.dump(full_galah_nlte, open(element + "_Galah_impacts2_all.pkl", "wb"))
                pickle.dump(star_mean_dict, open(element + "_Galah_impacts_lineless2_all.pkl", "wb"))

                pickle.dump([feh_list, nlte_list], open(element + "_NLTE_list2_all.pkl", "wb"))
            else:

                pickle.dump(full_galah_nlte, open(element + "_Galah_impacts2.pkl", "wb"))
                pickle.dump(star_mean_dict, open(element + "_Galah_impacts_lineless2.pkl", "wb"))

                pickle.dump([feh_list, nlte_list], open(element + "_NLTE_list2.pkl", "wb"))

    print("DONE")
    plt.scatter(feh_list, nlte_list, s=0.5, c="b")
    plt.xlabel("Feh")
    plt.ylabel("NLTE effect")
    plt.show()





element = "Ti2"

reliable     = pkl.load(open(element+"Accurate_param_indexes_galah.pkl", "rb"))
# Run once to make a dict of nlte and lte values for each star and line in A(Ti) form
#run_all(element)
# Then turn to Ti/H by taking the lte/nlte mean of solar twinms. We later turn it into Ti/Fe in other scripts
turn_ti_tih(element)
