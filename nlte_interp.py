import os
import sys
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from scipy.interpolate import interpn
from scipy.interpolate import interp1d
import warnings

#file,element,teff,logg,feh,xi = sys.argv[1],sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),float(sys.argv[6])

star = "Sun"
file = star+"_Ti1_ew"
element="Ti1"


"""
Star    Teff    Logg    Fe/H    Mic
Sun     5772    4.44    0       0.9
Arcturus 4286   1,64    -0.53   1.3
84937   6356    4.06    -2.06   1.2
140283  5792    3.65    -2.36   1.3
122563  4646    1.4     -2,5    1.8

"""
teff = 6324.3784
logg = 3.9686449
feh = -0.38920593
xi = 1.4542778
# paramS 6324.3784 3.9686449 1.4542778 -0.38920593

data = pd.read_csv(file,delim_whitespace=True)
index = data.index
ews = data["ew"]
ers = data["err"]
line = data["line"]
print("line", line)
pklfile = open('grids/{}.pkl'.format(element),'rb')
grid=pkl.load(pklfile)

if element == 'Na':
    solar=6.17
if element == 'Mg':
    solar=7.53
if element == 'Al':
    solar=6.37
if element == "Ti" or element == "Ti1" or element == "Ti2":
    solar = 4.90
tg = grid["tg"]
gg = grid["gg"]
xg = np.array(grid["xg"])
fg = np.array(grid["fg"])
ag = grid["ag"]+solar+feh
wg = grid["wg"]

wl = grid["wl"]
wn = grid["wn"]
ti = np.where(tg > teff)
ti = min(ti[0])-1
gi = np.where(gg > logg)
gi = min(gi[0])-1
fi = np.where(fg > feh)
fi = min(fi[0])-1
print(tg, gg, fg)

tg=tg[ti:ti+2]
gg=gg[gi:gi+2]
fg=fg[fi:fi+2]
print("suf", tg, gg, fg)

print("A", [ti, ti+2])
print("B", wn.shape)
print("C", wn[np.isfinite(wn)], len(wn[np.isfinite(wn)]))
print("C2", len(wl[np.isfinite(wl)]), len(wn[np.isnan(wn)]))
isin = wn[np.isfinite(wn)]

wl=wl[ti:ti+2,gi:gi+2,fi:fi+2,:,:]
print("D",ti)
wn=wn[ti:ti+2,gi:gi+2,fi:fi+2,:,:]

print("E", wn[np.isfinite(wn)], len(wn[np.isfinite(wn)]))
print("F", wl[np.isfinite(wl)], len(wl[np.isfinite(wl)]))
print("F2", len(wl[np.isnan(wl)]))
print("G", wl[0][0][0][0][0])
print("H", np.where(wl == -99), len(np.where(wl == -99)))

for windex in index:
    alte = -9
    anlte = -9
    elte = 9
    enlte = 9
    
    ew=ews[windex]
    print("EW?", ew, len(ews), np.log10(ew))
    er=ers[windex]
    print("fin", wl[np.isfinite(wl)])

    wli=wl[:,:,:,:,:,windex]
    wni=wn[:,:,:,:,:,windex]
    print("fin", wli[np.isfinite(wli)])
    wla=[]
    wna=[]
    print("ag", ag)
    for a in range(0,len(ag)):
        wla.extend(interpn([tg,gg,fg,xg,ag],wli,[teff,logg,feh,xi,ag[a]],method="linear",fill_value=None,bounds_error=False))
        print(wli.shape)
        print("paramas", [tg, gg, fg, xg, ag], [teff, logg, feh, xi, ag[a]])
        print("inter", ag[a], interpn([tg,gg,fg,xg,ag],wli,[teff,logg,feh,xi,ag[a]],method="linear",fill_value=None,bounds_error=False))
        print("inter2", ag[a], interpn([tg,gg,fg,xg,ag],wni,[teff,logg,feh,xi,ag[a]]))
        wna.extend(interpn([tg,gg,fg,xg,ag],wni,[teff,logg,feh,xi,ag[a]],method="linear",fill_value=None,bounds_error=False))
    wla=np.array(wla)
    wna=np.array(wna)
    print(1, wna[np.isfinite(wna)])

    #    if (ew !=0):
#        print(10**wla)
#        print(10**wna)
#        print(ag)
    if (ew != 0) & (len(wla[np.isfinite(wla)]) > 1):
        agl=ag[(np.isfinite(wla))]
        wll=wla[(np.isfinite(wla))]
        delta=10**wla[1:]-10**wla[0:-1]
        print("min", min(wll), "max", max(wll))
        print("AGL", agl)
        if min(delta)< -0.1:
            cut=min(np.where(delta < -0.1)[0])
            agl=agl[0:cut+1]
            wll=wll[0:cut+1]
            print("Cutting LTE",wg[windex],cut,10**wll[-1])
        if (np.log10(ew) < min(wll)):
            print("Ew not in LTE range for {}, adopting lower limit".format(wg[windex]))
            alte=min(agl)
            elte=9
        if (np.log10(ew) > max(wll)):
            print("Ew not in LTE range for {}, adopting upper limit".format(wg[windex]))
            alte=max(agl)
            elte=9
        if (np.log10(ew) > min(wll)) & (np.log10(ew) < max(wll)):
            print("ew within range1")
            flte=interp1d(wll,agl,kind='linear')
            alte=flte(np.log10(ew))
            elte=alte-flte(np.log10(ew-er))

    print("ew", ew, wna[np.isfinite(wna)])
    if (ew != 0 ) & (len(wna[np.isfinite(wna)]) > 1):

        agn=ag[(np.isfinite(wna))]
        wnn=wna[(np.isfinite(wna))]
        delta=10**wna[1:]-10**wna[0:-1]
        if min(delta)< -0.1:
            cut=min(np.where(delta < -0.1)[0])
            agn=agn[0:cut+1]
            wnn=wnn[0:cut+1]
#            print("Cutting NLTE",wg[windex],cut,10**wnn[-1])
        if (np.log10(ew) < min(wnn)):
            print("Ew not in NLTE range for {}, adopting lower limit".format(wg[windex]))
            anlte=min(agn)
            enlte=9
        if (np.log10(ew) > max(wnn)):
            print("Ew not in NLTE range for {}, adopting upper limit".format(wg[windex]))
            anlte=max(agn)
            enlte=9
        if (np.log10(ew) > min(wnn)) & (np.log10(ew) < max(wnn)):
            print("Ew within range")
            fnlte=interp1d(wnn,agn,kind='linear')
            anlte=fnlte(np.log10(ew))
            enlte=anlte-fnlte(np.log10(ew-er))

    print("Line,      EW,       EW error,  A(LTE),    Error,     A(NLTE),     Error,   Ion Imb")
    if (ew !=0) : #& (len(wla[np.isfinite(wla)]) > 1) & (len(wna[np.isfinite(wna)]) > 1):
        # line, equiv width, error in width,
        print("{} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f}".format(wg[windex],ew,er,alte,elte,anlte,enlte,anlte-alte))
    """else:
        if (ew !=0):
            print(10**wla,10**wna,ew)

            plt.plot(ag,wla)
            plt.show()"""
#plt.plot(flte(np.log10(ew)),np.log10(ew),color='red',marker='.')
#plt.show()
