import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
all = True
from matplotlib.colors import LinearSegmentedColormap

element = "Ti2"

all = False
def hr_diag_absolute():

    for element in ["Ti1", "Ti2"]:
        if "1" in element:
            lines = ["4758", "4759", "5689", "5739"]

            lines = [int(x) for x in lines]
            lines = np.asarray(lines)
        else:
            lines = ["4798", "4719", "4874"]
            lines = [int(x) for x in lines]
            lines = np.asarray(lines)

        fig = plt.figure()

        maxab = -100
        minab = 1000

        for line in lines:
            whr = np.where(line == lines)[0][0]
            print(whr)
            if all:
                if "1" in element:
                    number = 5
                    number2 = 3
                else:
                    number = 3
                    number2 = 1
            else:
                if "1" in element:
                    number = 4
                    number2 = 1
                else:
                    number = 3
                    number2= 1

            if whr == 0:
                ax = fig.add_subplot(number, number2, 1)
            elif whr == 1:
                ax = fig.add_subplot(number, number2, 2)
            elif whr == 2:
                ax = fig.add_subplot(number, number2, 3)
            elif whr == 3:
                ax = fig.add_subplot(number, number2, 4)
            elif whr == 4:
                ax = fig.add_subplot(number, number2, 5)
            elif whr == 5:
                ax = fig.add_subplot(number, number2, 6)
            elif whr == 6:
                ax = fig.add_subplot(number, number2, 7)
            elif whr == 7:
                ax = fig.add_subplot(number, number2, 8)
            elif whr == 8:
                ax = fig.add_subplot(number, number2, 9)
            elif whr == 9:
                ax = fig.add_subplot(number, number2, 10)
            elif whr == 10:
                ax = fig.add_subplot(number, number2, 11)
            elif whr == 11:
                ax = fig.add_subplot(number, number2, 12)
            elif whr == 12:
                ax = fig.add_subplot(number, number2, 13)
            elif whr == 13:
                ax = fig.add_subplot(number, number2, 14)
            elif whr == 14:
                ax = fig.add_subplot(number, number2, 15)
            elif whr == 15:
                ax = fig.add_subplot(number, number2, 16)
            elif whr == 16:
                ax = fig.add_subplot(number, number2, 17)

            # Teff logg feh vmic lteabund nlteIMPACT, _Galah_impacts2_reduced_galah
            impact = pkl.load(open(element + "_Galah_impacts_lineless2_Ati_reduced_vmic1.pkl", "rb"))
            impact2 = pkl.load(open(element + "_Galah_impacts_lineless2_Ati_reduced_vmic2.pkl", "rb"))
            impact = pkl.load(open(element + "_Galah_impacts2_Ati_reduced_vmic1.pkl", "rb"))
            impact2 = pkl.load(open(element + "_Galah_impacts2_Ati_reduced_vmic2.pkl", "rb"))
            nlteimpact = []
            lte_imbalance_list = []
            logg_list = []
            feh_list = []
            temp_list = []
            """load = False
            if load:
                lists = pkl.load(open("ion_imbalance_feh_nlte_lte2.pkl", "rb"))
                nlte_imbalance_list = lists[1]
                feh_list = lists[0]
                lte_imbalance_list=lists[2]
            else:"""
            # Make a list of imbalance in non lte and lte, including the feh and logg as we limit to stars with both
            # ti2 and ti1 lines (which is not done previously)
            for starn in impact:
                if starn in impact and starn in impact2:
                    pass
                else:
                    continue
                if line in impact[starn] and line in impact2[starn]:
                    pass
                else:
                    continue

                ti1lte = impact[starn][line][-2]
                ti1lte2 = impact2[starn][line][-2]
                diff = ti1lte2 - ti1lte
                if diff > 0.1:
                    print("YAAA", diff)
                if diff > maxab:
                    maxab = diff
                if diff < minab:
                    minab = diff
                ti1nlteimpact = impact[starn][line][-1]
                feh_list.append(impact[starn][line][2])
                nlteimpact.append(diff)
                logg_list.append(impact[starn][line][1])
                temp_list.append(impact[starn][line][0])
            """plt.scatter(temp_list, logg_list, c=nlteimpact, s= 2, cmap=LinearSegmentedColormap.from_list("white_viridis", [
                (0, '#ffffff'),
                (1e-20, '#ffffff'),
                (0.25, '#1f32db'),
                (0.5, '#740699'),
                (1, '#e21616'),
            ], N=1000))"""
            vmin = -0.01
            vmax = 0.01
            im = ax.scatter(temp_list, logg_list, c=nlteimpact, s=2,
                        cmap=LinearSegmentedColormap.from_list("white_viridis", [
                            (0, '#dea600'),
                            (1e-20, '#dea600'),
                            (0.1, '#dea600'),

                            (0.27, '#c5cae3'),
                            (0.5, '#1f32db'),
                            (0.75, '#740699'),
                            (1, '#e21616'),
                        ], N=1000), vmin=vmin, vmax=vmax)


            plt.gca().invert_yaxis()
            plt.gca().invert_xaxis()
            if "1" in element:
                plt.title("Ti I " + str(line)+"$\AA$")
            else:
                plt.title("Ti II " + str(line)+"$\AA$")
            # den = ax.scatter_density(feh, nlte, cmap=LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N).reversed())

            # den = ax.scatter_density(feh, nlte, cmap="Reds")


            #pkl.dump([feh_list, nlte_imbalance_list, lte_imbalance_list, logg_list], open("ion_imbalance_feh_nlte_lte2.pkl", "wb"))


            #plt.plot(trendfehgiants, trendgiants, label = "Trend line of Giant STD", color="blue")
            #plt.plot(trendfehdwarfs, trenddwarfs, label = "Trend line of Dwarf STD", color="black")
            if whr == number-1:
                plt.xlabel("T$_{eff}$ / K")
            else:
                ax.set_xticklabels([])
            if whr == 1:
                plt.ylabel("log(g / cm s$^{-2}$)")
        cbar_ax = fig.add_axes([0.82, 0.15, 0.025, 0.7])
        # fig.colorbar((den), cax=cbar_ax)
        fig.subplots_adjust(right=0.8)

        plt.savefig("graphs//" + str(element) + "HRdiagNLTE2_absolute_vmic.png", bbox_inches="tight")

        plt.colorbar(im, cax=cbar_ax, label="A(Ti)$_{\\xi = 2km s^{-1}}$ - A(Ti)$_{\\xi = 1km s^{-1}}$")
        print("max, min: ", maxab, minab)
        plt.show()
hr_diag_absolute()