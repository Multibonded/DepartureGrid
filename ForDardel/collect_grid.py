import os
import sys
import numpy as np
import pandas as pd
from scipy.integrate import simps
import pickle as pkl


#windex,element,ion = int(sys.argv[1]),sys.argv[2],sys.argv[3]

windex=0
element="Ti"
ion="1"




#home_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/NaMgAl/'
home_dir = '/home/jama2357/Documents/TiFeAtoms/DepartureGrid/CoG/'

tg = np.concatenate((np.linspace(3800,4000,3),np.linspace(4250,8000,16)))
gg = np.linspace(0.0,5.5,12)
xg = np.asarray([1.0,2.0])
fg = np.asarray([-5.,-4.,-3.,-2.5,-2.,-1.5,-1.,-0.75,-0.5,-0.25,0.0,0.25,0.5])
ag = np.linspace(-2,2,17)

data=pd.read_csv('{}{}_lines'.format(element,ion),dtype='str')
wg=np.asarray(data['ws'])


wl=np.zeros([len(tg),len(gg),len(fg),len(xg),len(ag)])+np.nan
wn=np.zeros([len(tg),len(gg),len(fg),len(xg),len(ag)])+np.nan
# wg[windex] is just the index appleid to the lines in the line file to represent thew avelength of the line we're looking at.
work_dir = home_dir + element + "/" + element + ion + '/{}{}_{}/'.format(element, ion, wg[windex])
print(work_dir)
lte=pd.read_csv(work_dir+'LTE_EW.txt',dtype='str',delim_whitespace=True,names=['x1','EW','IND','x2','x3','x4'])
nlte=pd.read_csv(work_dir+'NLTE_EW.txt',dtype='str',delim_whitespace=True,names=['x1','EW','IND','x2','x3','x4'])
lte_ew=np.asarray(lte['EW'],dtype='float64')
lte_ind=np.asarray(lte['IND'])
nlte_ew=np.asarray(nlte['EW'],dtype='float64')
nlte_ind=np.asarray(nlte['IND'])

ind = 0
for g in gg:
    for t in tg:
        if (g==0.0) & (t>=4500): 
            break
        if (g==0.5) & (t>=5250): 
            break
        if (g==1.0) & (t>=6000): 
            break
        if (g==1.5) & (t>=6750): 
            break
        if (g==2.0) & (t>=7500): 
            break
        if (g==5.5) & (t>=4000): 
            break
        for f in fg:
            for x in xg:
                for a in ag:
                    # When the index of the star we are looking at matches the ind counter (can be optimised easily)
                    # np where because there will usually be more thanb 1 index of a star I guess
                    #tupel, therefore we index to [0] for just array
                    #i = np.where(lte_ind == str(ind).zfill(5))
                    i = np.where(lte_ind == str(ind))[0]
                    if i.size:
                        if lte_ew[i] > 0.0:
                            print("here",(lte_ew[i]), t, g, f, x, a)
                            wl[np.where(tg == t),np.where(gg == g),np.where(fg == f),np.where(xg == x),np.where(ag == a)]=np.log10(lte_ew[i])
                        else:
                            wl[np.where(tg == t),np.where(gg == g),np.where(fg == f),np.where(xg == x),np.where(ag == a)]=-99

                    #i = np.where(nlte_ind == str(ind).zfill(5))
                    i = np.where(nlte_ind == str(ind))[0]

                    if i.size:

                        if nlte_ew[i] > 0.0:
                            print("there",(nlte_ew[i]), t, g, f, x, a)
                            print("share", np.where(tg == t))
                            wn[np.where(tg == t),np.where(gg == g),np.where(fg == f),np.where(xg == x),np.where(ag == a)]=np.log10(nlte_ew[i])
                        else:
                            wn[np.where(tg == t),np.where(gg == g),np.where(fg == f),np.where(xg == x),np.where(ag == a)]=-99
                    ind += 1
                    if (ind % 1000 == 0):
                        print(ind)
grid = {}
grid["tg"]= tg
grid["gg"]= gg
grid["xg"]= xg
grid["fg"]= fg
grid["ag"]= ag
grid["wg"]= wg
grid["wl"]= wl
grid["wn"]= wn
print(np.where(np.isfinite(wl)))
print(wl[np.where(np.isfinite(wl))])
print("{}{}_{}.pkl".format(element,ion,wg[windex]))
print(tg.shape)
outfile=open("grids/{}{}/{}{}_{}.pkl".format(element, ion, element,ion,wg[windex]),"wb")
pkl.dump(grid,outfile)
outfile.close()
