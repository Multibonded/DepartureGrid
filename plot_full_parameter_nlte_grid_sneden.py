import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt

def average():
    element = "Ti1"
    reduced = False
    abund_dict = pkl.load(open(element + "nlte_impact_dict_sneden.pkl", "rb"))
    ew_dict = pkl.load(open(element + "ew_dict_sneden.pkl", "rb"))
    print(abund_dict.keys())
    fig, ax = plt.subplots(nrows=3, ncols=3)
    index = 0
    atmos = pkl.load(open("atmosphere_nodes.pkl", "rb"))
    fig.set_figwidth(12)
    fig.set_figheight(14)

    highest = 0
    lowest = 1000000
    for teff in abund_dict:
        index+=1
        if teff not in [4000., 4500, 5500,  5000., 6000., 6500, 7000.]:
            index += 1
            continue
        for logg in abund_dict[teff]:
            if logg in [0, 0.5, 5.5]:
                index += 1
                continue
            logglist = [i for i in abund_dict[teff].keys()]
            graph = (np.where(logg == logglist)[0][0])
            if graph<5:
                row = 0
            elif graph<8:
                row = 1
            elif graph<11:
                row = 2
            elif graph < 14:
                row = 3
            col = (graph-2)%3


            allfehnltelist = []
            allfehlte = []
            fehlist = []

            for feh in abund_dict[teff][logg]:
                if [teff, logg, feh] not in atmos:
                    print(teff, logg, feh)
                    continue
                ltelist = []
                nltelist = []
                for vmic in abund_dict[teff][logg][feh]:
                    if vmic != 2:
                        index += 1
                        continue
                    for eps in abund_dict[teff][logg][feh][vmic]:
                        # Turns abunds into comprable lists (eps? enhancement?)
                        epslist = [i for i in (abund_dict[teff][logg][feh][vmic].keys())]
                        abslist =  epslist-feh- 4.9
                        if eps-feh-4.9 not in [0]:
                            index += 1
                            continue
                        if len(abund_dict[teff][logg][feh][vmic][eps]) < 2:
                            index += 1
                            continue

                        for line in abund_dict[teff][logg][feh][vmic][eps]:
                            print("LINE", line)
                            impact = abund_dict[teff][logg][feh][vmic][eps][line]
                            ew = ew_dict[teff][logg][feh][vmic][eps][line]
                            red_ew = (np.log10(((10 ** ew) / 1000) / line))
                            if -6 > red_ew > -4.9:
                                continue
                            if abs(impact) < 1E-5:
                                continue
                            """if np.isnan(impact):
                                continue
                                impact = 0"""

                            ltelist.append(impact)
                            nltelist.append((impact))

                        index += 1
                        if len(nltelist) == 0:
                            continue
                        if impact > highest:
                            higheststar = [teff, logg, feh, vmic, eps, line, impact+0.05]
                            highest = impact
                        if impact < lowest:
                            loweststar = [teff, logg, feh, vmic, eps, line, impact+0.05]
                            lowest = impact
                allfehlte.append(np.nanmean(ltelist))

                allfehnltelist.append(np.nanmean(nltelist))
                fehlist.append(feh)
                print(np.nanmean(nltelist))
                print(fehlist)
            if teff == 4000:
                style = "-"
            elif teff == 4500:
                style = "-."
            elif teff == 5000:
                style = ":"
            elif teff == 5500:
                style = "--"
            elif teff == 6000:
                style = "-"
            elif teff == 6500:
                style = "-."
            elif teff == 7000:
                style = ":"
            else:
                style = "-"
            # 0.05 from the corrections to A(Ti)
            if "2" in element:
                ax[row, col].plot(fehlist, np.asarray(allfehnltelist), linestyle=style, label = str(int(teff))+"K")
            else:
                ax[row, col].plot(fehlist, np.asarray(allfehnltelist), linestyle=style, label=str(int(teff)) + "K")

            #ax[row, col].scatter(fehlist, allfehlte, label = str(int(teff))+"K", linestyle = "--")
            if row ==2 and col == 2:
                ax[row, col].legend()
            ax[row,col].title.set_text("log(g) = "+str(logg))
            if "2" in element:
                ax[row, col].set(ylim = [-0.05,0.7], xlim = [-5.5,0.5])

            else:
                ax[row, col].set(ylim=[0, 0.7], xlim=[-5.5, 0.5])
            if row == 2 and col ==1:
                ax[row,col].set_xlabel("[Fe/H]", fontsize=15)
            if col == 0 and row ==1:
                ax[row, col].set_ylabel("$\Delta$ A(Ti)$_{\mathrm{Ti\, I}}$", fontsize=15)
    if "2" in element:
        plt.savefig("graphs/NLTE_PARAMETER_SPACES_Ti2_sneden.png", bbox_inches="tight")
    else:
        plt.savefig("graphs/NLTE_PARAMETER_SPACES_sneden.png", bbox_inches="tight")
    print(higheststar, "\n", loweststar)
    plt.show()

def individual():
    element = "Ti2"
    reduced = False
    abund_dict = pkl.load(open(element + "nlte_impact_dict_sneden.pkl", "rb"))
    ew_dict = pkl.load(open(element + "ew_dict_sneden.pkl", "rb"))
    print(abund_dict.keys())
    if "1" in element:
        linelist = [3741, 4533, 5192, 5210]
    else:
        linelist = [4443, 4468, 4501, 4571]

    for liner in linelist:
        highest = 0
        lowest = 1000000
        fig, ax = plt.subplots(nrows=3, ncols=3)
        index = 0
        atmos = pkl.load(open("atmosphere_nodes.pkl", "rb"))
        fig.set_figwidth(12)
        fig.set_figheight(14)

        for teff in abund_dict:
            index += 1
            if teff not in [4000., 4500, 5500, 5000., 6000., 6500, 7000.]:
                index += 1
                continue
            for logg in abund_dict[teff]:
                if logg in [0, 0.5, 5.5]:
                    index += 1
                    continue
                logglist = [i for i in abund_dict[teff].keys()]
                graph = (np.where(logg == logglist)[0][0])
                if graph < 5:
                    row = 0
                elif graph < 8:
                    row = 1
                elif graph < 11:
                    row = 2
                elif graph < 14:
                    row = 3
                col = (graph - 2) % 3

                allfehnltelist = []
                allfehlte = []
                fehlist = []

                for feh in abund_dict[teff][logg]:
                    if [teff, logg, feh] not in atmos:
                        print(teff, logg, feh)
                        continue
                    ltelist = []
                    nltelist = []
                    for vmic in abund_dict[teff][logg][feh]:
                        if vmic != 2:
                            index += 1
                            continue
                        for eps in abund_dict[teff][logg][feh][vmic]:
                            # Turns abunds into comprable lists (eps? enhancement?)
                            epslist = [i for i in (abund_dict[teff][logg][feh][vmic].keys())]
                            abslist = epslist - feh - 4.9
                            if eps - feh - 4.9 not in [0]:
                                index += 1
                                continue
                            if len(abund_dict[teff][logg][feh][vmic][eps]) < 2:
                                index += 1
                                continue

                            for line in abund_dict[teff][logg][feh][vmic][eps]:
                                if line != liner:
                                    continue
                                impact = abund_dict[teff][logg][feh][vmic][eps][line]
                                ew = ew_dict[teff][logg][feh][vmic][eps][line]
                                red_ew = (np.log10(((10 ** ew) / 1000) / line))
                                print(red_ew)
                                if -6 > red_ew > -4.9:
                                    continue
                                if abs(impact) < 1E-5:
                                    continue
                                """if np.isnan(impact):
                                    continue
                                    impact = 0"""

                                ltelist.append(impact)
                                nltelist.append((impact))

                            index += 1
                            if len(nltelist) == 0:
                                continue
                            if impact > highest:
                                higheststar = [teff, logg, feh, vmic, eps, line, impact + 0.05]
                                highest = impact
                            if impact < lowest:
                                loweststar = [teff, logg, feh, vmic, eps, line, impact + 0.05]
                                lowest = impact
                    allfehlte.append(np.nanmean(ltelist))

                    allfehnltelist.append(np.nanmean(nltelist))
                    fehlist.append(feh)
                    print(np.nanmean(nltelist))
                    print(fehlist)
                if teff == 4000:
                    style = "-"
                elif teff == 4500:
                    style = "-."
                elif teff == 5000:
                    style = ":"
                elif teff == 5500:
                    style = "--"
                elif teff == 6000:
                    style = "-"
                elif teff == 6500:
                    style = "-."
                elif teff == 7000:
                    style = ":"
                else:
                    style = "-"
                # 0.05 from the corrections to A(Ti)
                if "2" in element:
                    ax[row, col].plot(fehlist, np.asarray(allfehnltelist), linestyle=style, label=str(int(teff)) + "K")
                else:
                    ax[row, col].plot(fehlist, np.asarray(allfehnltelist), linestyle=style, label=str(int(teff)) + "K")

                # ax[row, col].scatter(fehlist, allfehlte, label = str(int(teff))+"K", linestyle = "--")
                if row == 2 and col == 2:
                    ax[row, col].legend()
                ax[row, col].title.set_text("log(g) = " + str(logg))
                if "2" in element:
                    ax[row, col].set(ylim=[-0.05, 0.7], xlim=[-5.5, 0.5])

                else:
                    ax[row, col].set(ylim=[0, 0.7], xlim=[-5.5, 0.5])
                if row == 2 and col == 1:
                    ax[row, col].set_xlabel("[Fe/H]", fontsize=15)
                if col == 0 and row == 1:
                    ax[row, col].set_ylabel("$\Delta$ A(Ti)$_{\mathrm{Ti\, I}}$", fontsize=15)
        plt.suptitle(liner)
        if "2" in element:
            plt.savefig("graphs/NLTE_PARAMETER_SPACES_Ti2_sneden" + str(liner) + ".png", bbox_inches="tight")
        else:
            plt.savefig("graphs/NLTE_PARAMETER_SPACES_sneden" + str(liner) + ".png", bbox_inches="tight")
        print(higheststar, "\n", loweststar)
        plt.show()