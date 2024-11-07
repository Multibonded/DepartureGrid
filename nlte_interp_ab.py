import os
import sys
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from scipy.interpolate import interpn
from scipy.interpolate import interp1d
import warnings

file,element,teff,logg,feh,xi,alind = sys.argv[1],sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),float(sys.argv[6]),int(sys.argv[7])

data = pd.read_csv(file,delim_whitespace=True)
index = data.index
ews = data["ew"]
ers = data["err"]
line = data["line"]
        
pklfile = open('grids/{}.pkl'.format(element),'rb')
grid=pkl.load(pklfile)

if element == 'Na':
    solar=6.17
if element == 'Mg':
    solar=7.53
if element == 'Al':
    solar=6.37

tg = grid["tg"]
gg = grid["gg"]
xg = np.array(grid["xg"])
fg = np.array(grid["fg"])
ag = grid["ag"][0:12]+solar+feh
wg = grid["wg"]
wl = grid["wl"][:,:,:,:,0:12,:]
wn = grid["wn"][:,:,:,:,0:12,:]

ti = np.where(tg > teff)
ti = min(ti[0])-1
gi = np.where(gg > logg)
gi = min(gi[0])-1
fi = np.where(fg > feh)
fi = min(fi[0])-1

tg=tg[ti:ti+2]
gg=gg[gi:gi+2]
fg=fg[fi:fi+2]

wl=wl[ti:ti+2,gi:gi+2,fi:fi+2,:,:]
wn=wn[ti:ti+2,gi:gi+2,fi:fi+2,:,:]

for windex in index: 

    alte = -9
    anlte = -9
    elte = 9
    enlte = 9
    
    ew=ews[windex]
    er=ers[windex]
    
    wli=wl[:,:,:,:,:,windex]
    wni=wn[:,:,:,:,:,windex]

    wla=[]
    wna=[]

    for a in range(0,len(ag)): 
        wla.extend(interpn([tg,gg,fg,xg,ag],wli,[teff,logg,feh,xi,ag[a]],method="linear",fill_value=None,bounds_error=False))
        wna.extend(interpn([tg,gg,fg,xg,ag],wni,[teff,logg,feh,xi,ag[a]],method="linear",fill_value=None,bounds_error=False))

    wla=np.array(wla)
    wna=np.array(wna)

    if (ew !=0):
        ew=10**wla[alind]
#        print(10**wla)
#        print(10**wna)
#        print(ag)

    if (ew != 0) & (len(wla[np.isfinite(wla)]) > 1):
        agl=ag[(np.isfinite(wla))]
        wll=wla[(np.isfinite(wla))]
        delta=10**wla[1:]-10**wla[0:-1]
        if min(delta)< 0:
            cut=min(np.where(delta < 0)[0])
            agl=agl[0:cut+1]
            wll=wll[0:cut+1]
#            print("Cutting ",wg[windex],cut,10**wll[-1])
        if (np.log10(ew) < min(wll)):
#            print("Ew not in LTE range for {}, adopting lower limit".format(wg[windex]))
            alte=min(agl)
            elte=9
        if (np.log10(ew) > max(wll)):
#            print("Ew not in LTE range for {}, adopting upper limit".format(wg[windex]))
            alte=max(agl)
            elte=9
        if (np.log10(ew) > min(wll)) & (np.log10(ew) < max(wll)):
            flte=interp1d(wll,agl,kind='linear')
            alte=flte(np.log10(ew))
            elte=alte-flte(np.log10(ew-er))

    if (ew != 0) & (len(wna[np.isfinite(wna)]) > 1):
        agn=ag[(np.isfinite(wna))]
        wnn=wna[(np.isfinite(wna))]
        delta=10**wna[1:]-10**wna[0:-1]
        if min(delta)< 0:
            cut=min(np.where(delta < 0)[0])
            agn=agn[0:cut+1]
            wnn=wnn[0:cut+1]
#        print(delta,10**wnn)
        if (np.log10(ew) < min(wnn)):
#            print("Ew not in NLTE range for {}, adopting lower limit".format(wg[windex]))
            anlte=min(agn)
            enlte=9
        if (np.log10(ew) > max(wnn)):
#            print("Ew not in NLTE range for {}, adopting upper limit".format(wg[windex]))
            anlte=max(agn)
            enlte=9
        if (np.log10(ew) > min(wnn)) & (np.log10(ew) < max(wnn)):
            fnlte=interp1d(wnn,agn,kind='linear')
            anlte=fnlte(np.log10(ew))
            enlte=anlte-fnlte(np.log10(ew-er))
            
    if (ew !=0) : #& (len(wla[np.isfinite(wla)]) > 1) & (len(wna[np.isfinite(wna)]) > 1):
        print("{} {:10.3f} {:10.3f} {: .3f} {: .3f} {: .3f} {: .3f} {: .3f}".format(wg[windex],ew,er,alte,elte,anlte,enlte,anlte-alte))
#    else:
#        if (ew !=0):
#            print(10**wla,10**wna,ew)

            #plt.plot(ag,wla)
#plt.plot(flte(np.log10(ew)),np.log10(ew),color='red',marker='.')
#plt.show()
