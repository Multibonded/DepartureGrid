import os
import sys
import glob
import numpy as np
import pandas as pd
import pickle as pkl

#windex = int(sys.argv[1])
home_dir='/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/NaMgAl'

element='Fe'
tg = np.concatenate((np.linspace(3800,4000,3),np.linspace(4250,8000,16)))
gg = np.linspace(0.0,5.5,12)
xg = [1,2]
fg = [-5,-4,-3,-2.5,-2,-1.5,-1,-0.75,-0.5,-0.25,0.0,0.25,0.5]
#io = np.ones(32)
#ev = np.zeros(32)
#gf = np.zeros(32)
data=pd.read_csv('{}_lines'.format(element),delim_whitespace=True)
wg=np.asarray(data['ws'],dtype='str')
#wg=wg.strip()

wl=np.zeros([len(tg),len(gg),len(fg),len(xg),len(wg)])+np.nan
wn=np.zeros([len(tg),len(gg),len(fg),len(xg),len(wg)])+np.nan

ind = 0
for w in wg:
    filename=glob.glob('grids/{}[12]_{}.pkl'.format(element,w))
    filename=filename[0]
    file=os.path.join(home_dir,filename)
    if os.path.isfile(file):
        print(file)
        linefile = open(file,'rb')
        data=pkl.load(linefile)
        wlline=data['wl']
        wl[:,:,:,:,ind]=wlline[:,:,:,:]
        wnline=data['wn']
        wn[:,:,:,:,ind]=wnline[:,:,:,:]
    ind += 1
    print(ind)
                        

grid = {}
grid["tg"]= tg
grid["gg"]= gg
grid["xg"]= xg
grid["fg"]= fg
#grid["io"]= io
#grid["ev"]= ev
#grid["gf"]= gf
grid["wg"]= wg
grid["wl"]= wl
grid["wn"]= wn

outfile=open("grids/{}.pkl".format(element),"wb")
pkl.dump(grid,outfile)
outfile.close()
