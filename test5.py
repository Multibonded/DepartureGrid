"""Applies NLTE effects to GALAH stars via interpolation. Produces a few files that has the info for all
galah stars and their nlte impact """
import os
import pickle
import sys
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from scipy.interpolate import interpn
from scipy.interpolate import interp1d
import warnings

from astropy.io import fits

atomname = r"/home/jama2357/Downloads/galahdr3_abundance_zeropoints.fits"
atom = fits.open(atomname)
atom1 = atom[1].data
print(repr(atom[1].header))
print(atom1['A_Ti4758'])

print(atom1['A_Ti4759'])
print(atom1['A_Ti5689'])
print(atom1['A_Ti5716'])
print(atom1['A_Ti5739'])


print(atom1['A_Ti4720'])

print(atom1['A_Ti4799'])
print(atom1['A_Ti4874'])

full_galah_nlte = {}
