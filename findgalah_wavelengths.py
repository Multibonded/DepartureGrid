
"""Outdated, use Open_Galah_linelist (in pytools Jack) instead as that comes frm their own source"""
import os
import sys
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from scipy.interpolate import interpn
from scipy.interpolate import interp1d
import warnings

from astropy.io import fits
atomname = r"/home/jama2357/Documents/TiFeAtoms/DepartureGrid/galah_master_v5.2.fits"
atom = fits.open(atomname)
atom1 = (atom[1].data)
atom0 = (atom[0].data)
print(repr(atom[1].header))
lam = (atom1['lambda'])
name = atom1['name']
ion = atom1['ion']
y = 0
lambda_list = []
energy_list = []

tilam = []
for x in range(len(atom1['name'])):
    if name[x][0] == "Ti" and not name[x][1] and not name[x][2]:
        tilam.append(lam[x])
        y+=1

print(y)
starname = r"/home/jama2357/Documents/TiFeAtoms/DepartureGrid/GALAH_DR3_main_allstar_v2.fits"
"""stars = fits.open(starname)
star1= (stars[1].data)
star0 = (stars[0].data)"""

with open("/home/jama2357/Documents/TiFeAtoms/DepartureGrid/ges_clean_v6.txt") as f:
    for line in f:
        if line.split()[0] == "Ti" and "Y" in line.split()[5]:
            #print(line.split())
            lambda_list.append(float(line.split()[2]))



print(len(lambda_list))
lambda_list=np.asarray(lambda_list)

# check which ges clean lines are in the galah main file. Fewer than thought
double_lambda = []
for l in lambda_list:

    #print(l, min(tilam, key=lambda x:abs(x-l)), l in tilam)
    if l in tilam and l != 5702.66 and l!= 5716.45:
        double_lambda.append(l)

print(double_lambda)
print(len(double_lambda))
for x in double_lambda:
    a = np.where(lam==x)
    print(atom1[a])
    spec = atom1[a]
    species = spec['name'][0]
    atom_number=22
    ionization = spec['ion']
    wlcent = spec['lambda']
    excit = spec['E_LOW']
