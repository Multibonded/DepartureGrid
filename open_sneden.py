import numpy as np
import os
def find_lines_in_all_stars():
    ti1lines = []
    ti2lines = []
    ti1all = []
    ti1dict = {}
    ti2all = []
    ti2dict = {}
    for filename in os.listdir("/home/jama2357/Documents/TiFeAtoms/DepartureGrid/CoG/sneden_stars/"):
        ti1this = []
        ti2this = []
        print(filename)
        with open("/home/jama2357/Downloads/"+filename) as file:
            for line in file:
                splitline  = (line.split())
                if splitline[0] =="Ti" and splitline[1] == "I":
                    ti1lines.append(splitline[2])
                    ti1this.append(splitline[2])
                if splitline[0] =="Ti" and splitline[1] == "II":
                    ti2lines.append(splitline[2])
                    ti2this.append(splitline[2])

        ti1dict[filename] = ti1this
        ti2dict[filename] = ti2this
    print(ti1dict.values())
    print(set(ti1lines))
    print(set(ti2lines))
    print(len(set(ti1lines)))
    print(len(set(ti2lines)))
    print(len(ti1dict.keys()))
    for line in ti1dict["g18128.lin"]:
        fail = False
        failure = 0
        for key in ti1dict.keys():
            if line not in ti1dict[key]:
                failure +=1
                if failure > 20:
                    fail = True
        if fail == False and line not in ti1all:
            ti1all.append(line)

    print(ti1all)

    for line in ti2dict["g18128.lin"]:
        fail = False
        for key in ti2dict.keys():
            if line not in ti2dict[key]:
                fail = True
        if fail == False and line not in ti2all:
            ti2all.append(line)

    print("2", ti2all)






def save_abunds():
    for element in ["Ti1", "Ti2"]:
        if "1" in element:
            lines = ['3354.63', '3653.49', '3741.06', '4533.24', "5192.97"]
        else:
            lines = ['3409.81', '4443.80', '4468.49', '4501.27', '4571.97']

        tidict = {}
        for filename in os.listdir("/home/jama2357/Documents/TiFeAtoms/DepartureGrid/CoG/sneden_stars/"):
            star = (filename.replace(".lin", ""))
            with open("/home/jama2357/Downloads/"+filename) as file:
                for line in file:
                    splitline  = (line.split())
                    if "/" in splitline[0]:
                        print(splitline)
                        teff = float(splitline[0].replace("/", ""))
                        logg = float(splitline[1].replace("/", ""))
                        feh = float(splitline[2])
                        vmic = float(splitline[-1])
                    if splitline[0] =="Ti" and splitline[2] in lines:

                        print(int(np.floor(float(splitline[2]))))
                        print(splitline)
                        epsilon = float(splitline[-1])
                        if star not in tidict:
                            tidict[star] = {}
                        tidict[star][int(np.floor(float(splitline[2])))] = [teff, logg, feh, vmic, epsilon]
        print(tidict)
        import pickle
        pickle.dump(tidict, open(element+"sneden_abundances.pkl", "wb"))
save_abunds()
