"""see impact of vmic on stars with red EW above -4.9"""
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
all = False
element = "Ti1"
if all:
    abund_dict = pkl.load(open(element + "vmic_abund_dict_grid_all.pkl", "rb"))
    vmic_dict = pkl.load(open(element + "vmic_impact_dict_all.pkl", "rb"))
    ew_dict = pkl.load(open(element + "vmic_ew_dict_all.pkl", "rb"))

else:
    abund_dict = pkl.load(open(element + "vmic_abund_dict_grid.pkl", "rb"))
    vmic_dict = pkl.load(open(element + "vmic_impact_dict.pkl", "rb"))
    ew_dict = pkl.load(open(element + "vmic_ew_dict.pkl", "rb"))
sub49 = []
plus49 = []
subew = []
plusew = []
subabs = []
plusabs = []

fig = plt.figure()
if "1" in element:
    lines = ["5739", "4758", "4759", "4778", "4781", "4797", "4801", "4820", "5689", "5716", "5720",
             "5866", "6716", "7852"]
    lines = ["4758", "4759", "5689", "5739"]


    lines = [int(x) for x in lines]
    lines = np.asarray(lines)
else:
    lines = ["4874", "4865", "4849", "4798", "4764", "4719"]
    lines = ["4719", "4798", "4874"]

    lines = [int(x) for x in lines]
    lines = np.asarray(lines)

for line_X in lines:
    sub49 = []
    plus49 = []
    subew = []
    plusew = []
    subabs = []
    plusabs = []

    whr = np.where(line_X == lines)[0][0]
    print(whr, line_X)
    if "1" in element:
        elementstring = "Ti I"
        number = 5
        number2= 3
        if not all:
            number = 4
            number2 = 1

    else:
        elementstring = "Ti II"

        number = 3
        number2 = 1
        if not all:
            number = 3
            number2 = 1

    if whr == 0:
        ax = fig.add_subplot(number, number2,  1)
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

    line = int(line_X)

    for teff in abund_dict:
        for logg in abund_dict[teff]:
            for feh in abund_dict[teff][logg]:

                for vmic in abund_dict[teff][logg][feh]:

                    if vmic != 1.:
                        continue
                    zerokey = round(0.75+feh+4.9, 2)
                    abundances = list(abund_dict[teff][logg][feh][vmic])

                    if (line) in abund_dict[teff][logg][feh][vmic][zerokey]:
                        red1 = (np.log10(((10 ** ew_dict[teff][logg][feh][1.][zerokey][line]) / 1000) / line))
                        red2 = (np.log10(((10 ** ew_dict[teff][logg][feh][2.][zerokey][line]) / 1000) / line))
                        abund1 = abund_dict[teff][logg][feh][1.][zerokey][line]
                        abund2 = abund_dict[teff][logg][feh][2.][zerokey][line]
                        
                        vmic1 = vmic_dict[teff][logg][feh][1.][zerokey][line]
                        vmic2 = vmic_dict[teff][logg][feh][2.][zerokey][line]
                        diff = abund2 - abund1

                        if np.isnan(diff):
                            continue
                        if teff == 5500 and logg == 1.0 and feh == -.75 and line == 4759:
                            testab1 = abund_dict[teff][logg][feh][1.][zerokey][line]
                            testab2 = abund_dict[teff][logg][feh][2.][zerokey][line]
                            nltetest1 = vmic_dict[teff][logg][feh][1.][zerokey][line]
                            nltetest2= vmic_dict[teff][logg][feh][2.][zerokey][line]
                            print("HERE SIR", testab1, testab2, testab1 - nltetest1, testab2 - nltetest2,
                                  diff, 10 **(ew_dict[teff][logg][feh][1.][zerokey][line])/1000, (10 **ew_dict[teff][logg][feh][2.][zerokey][line])/1000, red1, red2)
                            print(list(abund_dict[teff][logg][feh][vmic]), 3.65-feh-4.9)
                        if red1 > -4.9:
                            plus49.append((vmic1))
                            plusew.append(red1)
                            plusabs.append((abund1))
                        else:
                            sub49.append((vmic1))
                            subew.append(red1)
                            subabs.append((abund1))

    print(sub49)
    print("max", max(sub49))
    print("max above 4.9", max(plus49))
    type = "delta"
    if type == "delta":
        if all:
            if str(line) in ["4758", "4759", "5689", "5739"]:
                ax.scatter(subew, sub49, s=1, c="green", label = "Line: " + str(line) + "$\AA$")

            else:
                ax.scatter(subew, sub49, s=1, label = "W$_{red,\, \lambda}$ < -4.9")
        else:
            ax.scatter(subew, sub49, s=1, label = "W$_{red,\, \lambda}$ < -4.9")

        ax.scatter(plusew, plus49, s =1,  label = "W$_{red,\, \lambda}$ > -4.9")
        plt.legend()

        if whr == 3:
            plt.xlabel("W$_{red,\, \lambda}$")
        if "2" in element:
            if whr == 2:
                plt.xlabel("W$_{red,\, \lambda}$")

        #plt.title("Impact of Vturb for line" + str(line))
        if whr == 2 and "1" in element:
            plt.ylabel("A(Ti)$_{\\xi = 1km s^{-1}}$ - A(Ti)$_{\\xi = 2km s^{-1}}$")
        elif "2" in element and whr == 1:
            plt.ylabel("A(Ti)$_{\\xi = 1km s^{-1}}$ - A(Ti)$_{\\xi = 2km s^{-1}}$")

        plt.ylim(-2,2)
        plt.xlim(-7,-4)
        plt.savefig("graphs/vmictest" + str(element) + "png")

    if type == "abundance":
        ax.scatter(subew, subabs, s =1, label = "Line: " + str(line) + "$\AA$")
        ax.scatter(plusew, plusabs, s =1)
        if whr >= 11:

            plt.xlabel("Reduced EW")
        if whr == 0 or whr == 2:

            plt.ylabel("A(Ti)")
        plt.legend()
        plt.savefig("graphs/vabundtest" + str(element) + str(line)+".png")
    if whr == 0:
        plt.title(elementstring+": "+ str(line) + "$\AA$")
    else:
        plt.title(str(line) + "$\AA$")

plt.show()
