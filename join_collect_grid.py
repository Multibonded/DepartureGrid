import os
import sys
import glob
import numpy as np
import pandas as pd
import pickle as pkl

#home_dir='/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/NaMgAl'
home_dir='/home/jama2357/Documents/TiFeAtoms/DepartureGrid/CoG/'

#element=sys.argv[1]
all = True
element="Ti1"
reduced = False
tg = np.asarray(np.concatenate((np.linspace(3800,4000,3),np.linspace(4250,8000,16))))
gg = np.asarray(np.linspace(0.0,5.5,12))
xg = np.asarray([1,2])
fg = np.asarray([-5,-4,-3,-2.5,-2,-1.5,-1,-0.75,-0.5,-0.25,0.0,0.25,0.5])
ag = np.asarray(np.linspace(-2,2,17))
if all:
    data=pd.read_csv('{}_lines_all'.format(element),delim_whitespace=True)

else:
    data=pd.read_csv('{}_lines'.format(element),delim_whitespace=True)
wg = np.asarray(data['line'],dtype='str')
print("wg", wg)
wl = np.zeros([len(tg),len(gg),len(fg),len(xg),len(ag),len(wg)])+np.nan
wn = np.zeros([len(tg),len(gg),len(fg),len(xg),len(ag),len(wg)])+np.nan
ind = 0
# w is the line in AA
for w in wg:
    # NOT A TI1 LINE
    print("w", w)
    #filename=glob.glob('grids/{}[12]_{}.pkl'.format(element,w))
    if all:
        filename = glob.glob('grids/{}/{}_{}_all.pkl'.format(element, element, w))
    else:
        filename = glob.glob('grids/{}/{}_{}.pkl'.format(element, element, w))
    print('grids/{}[12]_{}.pkl'.format(element,w))
    print("Hhh", filename)
    filename=filename[0]
    file=os.path.join(home_dir,filename)
    print(os.path.isfile(file))
    if os.path.isfile(file):
        print(file)
        linefile = open(file,'rb')
        data=pkl.load(linefile)
        wlline=data['wl']
        wl[:,:,:,:,:,ind]=wlline[:,:,:,:,:]
        wnline=data['wn']
        wn[:,:,:,:,:,ind]=wnline[:,:,:,:,:]

    else:
        print(file, "is not a file, are you sure this is correct?\n--------------------")
    ind += 1
    print(ind)
# what is ag? Why do we do this? limits to -1 to 0.75. Enhancement I think ,so abundance grid = ag?
ag = ag[4:12]

print("ag", ag)
wl = wl[:,:,:,:,4:12,:]
wn = wn[:,:,:,:,4:12,:]

grid = {}
grid["tg"]= tg
grid["gg"]= gg
grid["xg"]= xg
grid["fg"]= fg
grid["ag"]= ag
grid["wg"]= wg
grid["wl"]= wl
grid["wn"]= wn

if all:
    outfile = open("grids/{}_all.pkl".format(element), "wb")
else:
    outfile=open("grids/{}.pkl".format(element),"wb")

pkl.dump(grid,outfile)
outfile.close()

original_stdout = sys.stdout # Save a reference to the original standard output
max = []
nlte = 0
print("Writing out now")
with open('grids/{}.txt'.format(element), 'w') as file:
    sys.stdout = file
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
                        for w in wg:
                            finwl=np.isfinite(wl[np.where(tg == t),np.where(gg == g),np.where(fg == f),np.where(xg == x),np.where(ag == a),np.where(wg == w)])
                            finwn=np.isfinite(wn[np.where(tg == t),np.where(gg == g),np.where(fg == f),np.where(xg == x),np.where(ag == a),np.where(wg == w)])
                            poswl=wl[np.where(tg == t),np.where(gg == g),np.where(fg == f),np.where(xg == x),np.where(ag == a),np.where(wg == w)] > 0
                            poswn=wn[np.where(tg == t),np.where(gg == g),np.where(fg == f),np.where(xg == x),np.where(ag == a),np.where(wg == w)] > 0

#                            print(finwl,finwn,poswl,poswn)
                            if finwl and finwn and poswl and poswn:
                                print("{} {:7s} {:4.0f} {:4.1f} {:5.2f} {:4.1f} {:4.1f} {:5.3f} {:5.3f}".format(element,w,t,g,f,x,a,*wl[np.where(tg == t),np.where(gg == g),np.where(fg == f),np.where(xg == x),np.where(ag == a),np.where(wg == w)][0],*wn[np.where(tg == t),np.where(gg == g),np.where(fg == f),np.where(xg == x),np.where(ag == a),np.where(wg == w)][0]))

sys.stdout = original_stdout

