
import numpy as np
more = open("/home/jama2357/Documents/TiFeAtoms/DepartureGrid/CoG/Ti/Ti1/Ti1_5739/LTE_EW.txt", "r").read()
one = open("/home/jama2357/Documents/TiFeAtoms/DepartureGrid/CoG/Ti/Ti1/Ti1_5739_single/LTE_EW.txt", "r").read()
mored={}
oned={}
for x in more.split("\n"):
    sp = x.split()

    if len(sp) == 0:
        continue
    star = sp[-1]

    mored[star] = float(sp[1])
for x in one.split("\n"):
    sp = x.split()

    if len(sp) == 0:
        continue
    star = sp[-1]
    oned[star] = float(sp[1])

combd = {}
diffl = []
for k in oned.keys():
    onea = oned[k]/1000
    morea = mored[k]/1000

    diff = np.log10(morea)-np.log10(onea)
    if diff > 0.05:
        print(onea, diff, np.log10(onea/5739))
    if np.isnan(diff):
        continue
    if onea < 1E-4:
        continue
    diffl.append(diff)
print(sorted(diffl))
