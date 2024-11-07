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
file = "Sun_Ti1_ew"
element="Ti1"
teff = 5772
logg = 4.44
feh = 0.01
xi = 0.9

# Reading in equivalent width file of format line, ew, error
data = pd.read_csv(file,delim_whitespace=True)
index = data.index
ews = data["ew"]
ers = data["err"]
line = data["line"]
        
pklfile = open('grids/{}.pkl'.format(element),'rb')
grid=pkl.load(pklfile)
if element == 'Fe':
    solar=7.45

tg = grid["tg"]
gg = grid["gg"]
xg = np.array(grid["xg"])
fg = np.array(grid["fg"])
wg = grid["wg"]
# I assume wl is wavelength. There's a tonne of "nan". Grid of wavelengths for all parameters? We later index to only
# our chosen parameters for Ti and Gi, but not fi?? The non iron file uses Fi, so that's weird.
# What's the 2ndtofinal column? Perhaps microturb. We always take all of them.
# Spoiler: wl is not wavelength. It contains a huge amount of info. 19, 12, 13, 2, 18, so the parameters possible
# and then the index of the lines. I don't know if windex actually DOES anything but maybe for another star.
# So what do the values themselves represent
wl = grid["wl"]
wn = grid["wn"]

# ti will give the index of the minimum teff for interpolation (so start point and later, ti+2 is the ceiling for
# maximum teff (Due to how arrays work, it then takes the values from index ti and ti+1)
# Same for gravity (g) and metallicity (f)
ti = np.where(tg > teff)
ti = min(ti[0])-1
gi = np.where(gg > logg)
gi = min(gi[0])-1
fi = np.where(fg > feh)
fi = min(fi[0])-1
print(tg)
tg=tg[ti:ti+2]
print(tg)
gg=gg[gi:gi+2]
#fg=fg[fi:fi+2]
print(tg,gg,fg)

# Takes metallicity limitations for non iron scripts. Final column or 2ndlast may be microturb.
# Wn has the majority of the "nan" values./ Can't find any real numbers on a quick check.
# Has 5 nested lists. Temp, Grav, Metal, micro?, what's the last one? Wijndex, so the indx of the line?
wl=wl[ti:ti+2,gi:gi+2,:,:,:]
wn=wn[ti:ti+2,gi:gi+2,:,:,:]
#wl=wl[ti:ti+2,gi:gi+2,fi:fi+2,:,:]
#wn=wn[ti:ti+2,gi:gi+2,fi:fi+2,:,:]

for windex in index: 

    flte = -9
    elte = 9
    
    ew=ews[windex]
    er=ers[windex]
    # How odd that we never indexed this column previously to match our teff etc, but now we index it based on... what?
    # Line number? Windex is just a run through of index, which I think is just the number of the line in question (not
    # wavelength)
    wli=wl[:,:,:,:,windex]
    wla=[]
    # fg is the list of metallicities. Long in the iron list, limited in non-iron
    for f in range(0,len(fg)): # What are we interpolating here?
        # We use the parameters as the grid points
        # Wli are the values for these grid points. So param=x and wli=y?
        # the coordinates of the array is [teff, logg, fg[f], xi], so where to sample the data from..? What does that
        # even MEAN. Wla/wll is used to find minimum and maximum LTE limits for ew...
        # AH OK it's getting wli, wavelengths that we've indexed to find the lines we have in our input file (e.g for
        # the sun) and finding the values at the parameters needed. wla contains logged equiv widths!!
        wla.extend(interpn([tg,gg,fg,xg],wli,[teff,logg,fg[f],xi],method="linear",fill_value=None,bounds_error=False))

    wla=np.array(wla)
#    if (ew !=0):
#        print(10**wla)
#        print(fg)

    # Where we have an EW calculated (Why so many nans?)
    if (ew != 0) & (len(wla[np.isfinite(wla)]) > 1):
        # fininite lists for metalicity and wavelengths
        fgl=fg[(np.isfinite(wla))]
        wll=wla[(np.isfinite(wla))]
        # So we have wll/wla as the range of equivalent widths for the range of parameters between fg and fg+2
        # If no metallicity did this, we take the min value (poor thing to do)
        if (np.log10(ew) < min(wll)):
            print("Ew not in LTE range for {}, adopting lower limit".format(wg[windex]))
            flte=min(fgl)
            elte=9
        if (np.log10(ew) > max(wll)):
            print("Ew not in LTE range for {}, adopting upper limit".format(wg[windex]))
            flte=max(fgl)
            elte=9
        # If within the range, we then have the equivalent width for the star at these parameters?
        if (np.log10(ew) > min(wll)) & (np.log10(ew) < max(wll)):
            # and now we interpolate to find these values at the metallicity we want
            funclte=interp1d(wll,fgl,kind='linear')
            # log it to get logged ew and find the metallicity that would give the ew we see? in lte though?
            flte=funclte(np.log10(ew))
            # error in logged ew lte
            elte=flte-funclte(np.log10(ew-er))
    # Now finally print the line, its observed ew and er, and the lte interpolated value of it. But I think EW isn't
    # logged so that's dumb
    if (ew !=0) : #& (len(wla[np.isfinite(wla)]) > 1) & (len(wna[np.isfinite(wna)]) > 1):
        print("{} {:10.3f} {:10.3f} {: .3f} {: .3f}".format(wg[windex],ew,er,flte,elte))
