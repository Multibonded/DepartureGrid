import os
import sys
import numpy as np
import pandas as pd
from scipy.integrate import simps
import pickle as pkl

windex,element,ion = int(sys.argv[1]),sys.argv[2],sys.argv[3]

home_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/NaMgAl/'

tg = np.concatenate((np.linspace(3800,4000,3),np.linspace(4250,8000,16)))
gg = np.linspace(0.0,5.5,12)
xg = np.asarray([1.0,2.0])
fg = np.asarray([-5.,-4.,-3.,-2.5,-2.,-1.5,-1.,-0.75,-0.5,-0.25,0.0,0.25,0.5])

data=pd.read_csv('{}{}_lines'.format(element,ion),dtype='str')
wg=np.asarray(data['ws'])

wl=np.zeros([len(tg),len(gg),len(fg),len(xg)])+np.nan
wn=np.zeros([len(tg),len(gg),len(fg),len(xg)])+np.nan

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
                w=wg[windex]
                filename='{}{}_{}/{}_{:.0f}_{}_{}_{}_LTE_compressed_5.txt'.format(element,ion,w,str(ind).zfill(5),t,float(g),float(f),float(x))
                file=os.path.join(home_dir,filename)
                if os.path.isfile(file):
                    data = pd.read_csv(file)
                    if len(data['wave']) > 4:
                        wave=np.asarray(data['wave'])
                        wave=wave[1:-1]
                        synth=np.asarray(data['synth'])
                        synth=synth[1:-1]
                        ew=simps(1-synth,wave)*1e3
#                        ew=simps(1-data['synth'],data['wave'])*1e3
                        if ew > 0:
                            wl[np.where(tg == t),np.where(gg == g),np.where(fg == f),np.where(xg == x)]=np.log10(ew)
                            
                filename='{}{}_{}/{}_{:.0f}_{}_{}_{}_NLTE_compressed_5.txt'.format(element,ion,w,str(ind).zfill(5),t,float(g),float(f),float(x))
                file=os.path.join(home_dir,filename)
                if os.path.isfile(file):
                    data = pd.read_csv(file)
                    if len(data['wave']) > 2:
                        ew=simps(1-data['synth'],data['wave'])*1e3
                        wn[np.where(tg == t),np.where(gg == g),np.where(fg == f),np.where(xg == x)]=np.log10(ew)
                ind += 1
                if (ind % 100 == 0):
                    print(ind)
                        

grid = {}
grid["tg"]= tg
grid["gg"]= gg
grid["xg"]= xg
grid["fg"]= fg
grid["wg"]= wg
grid["wl"]= wl
grid["wn"]= wn

outfile=open("{}{}_{}.pkl".format(element,ion,wg[windex]),"wb")
pkl.dump(grid,outfile)
outfile.close()
