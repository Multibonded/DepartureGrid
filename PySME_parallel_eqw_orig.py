import sys 
import os.path  
# Use sys.path.insert to /src/ if .bash_profile is not modified with: export PYTHONPATH=/proj/snic2020-16-23/private/USERNAME/SME-master/src/:$PYTHONPATH
from pysme import sme as SME #Standard sme imports
from pysme.abund import Abund 
from pysme.linelist.vald import ValdFile

from pysme.synthesize import synthesize_spectrum #Use to just synthesize spectra without fitting, otherwise sme.solve
from pysme import util #Used for logging file

import numpy as np #Using numpy to set up sme.wran and sme.wave, not needed if working from input files
import glob #Using glob to identify the correct linelist, needed because the names are not formatted identically; not all linelists have hfs and I need to remember
import pandas as pd

#from eqwidth import eqwidth

from sme_compressor import compress
from scipy.integrate import simps

"______________________________________________________"
# Arguments from bash script
windex = int(sys.argv[1])

element, stars, NLTE, ion, line= sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]
home_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/NaMgAl/'
work_dir = home_dir + '{}{}_{}/'.format(element,ion,line)

# Trying out this pandas method to use seq for parallel, if it's too slow (approx 0.06 sec) I'll just use cat instead
# _params_parallel_filt are the parameters that don't need interpolation
params = pd.read_csv(home_dir + stars +  '_params_parallel.txt',delim_whitespace = True,names=['i','teff','logg','monh','vmic','enh'])

our_params = params.loc[params.i == windex].values[0]
index, teff, logg, monh, vmic, enh = int(our_params[0]), int(our_params[1]), round(our_params[2],1), round(our_params[3],2), round(our_params[4],1), round(our_params[5],2)
EleIon = element + ion

# target specifies filenames, index used to keep track of progress
target = "{}_{}_{}_{}_{}_{}".format(str(index).zfill(5),teff,logg,monh,vmic,enh)

if NLTE == 'y':
    target = target + "_NLTE"
    print("\n Synthesizing {} in NLTE for parameters {}".format(element,str(index).zfill(5)))
else:
    target = target + "_LTE"
    print("\n Synthesizing {} in LTE for parameters {}".format(element,str(index).zfill(5)))

line_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/Linelists/'
NLTE_dir = "/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/NLTE_grids/"

out_file = os.path.join(work_dir, f"{target}.sme")
npy_file = os.path.join(work_dir, f"{target}.npy")
log_file = os.path.join(work_dir, f"{target}.log")
ew_file = os.path.join(work_dir, f"{target}.csv")
#util.start_logging(log_file)

"______________________________________________________"
sme = SME.SME_Structure()

sme.teff, sme.logg, sme.monh, sme.vmic = teff, logg, monh, vmic

# The 0 sets abundances to solar, Abund(0,"grevesse12007") gets the same solar abundances used to compute the 2014 MARCS grid
solar_abund = Abund(0,"grevesse2007").__getitem__('{}'.format(element))
solar_fe = Abund(0,"grevesse2007").__getitem__('Fe')

# Sets abundances of all elements according to the metallicity adn Grevesse 2007. Calculates the element abundance in the A12 format from the [A/Fe] formatted enh parameter. 
sme.abund = Abund(monh,"grevesse2007")
stellar_fe = sme.abund.__getitem__('Fe')
stellar_abund = round(enh + solar_abund,2) # Using round to avoid NLTE grid extrapolation
#setitem subtracts monh 
sme.abund.__setitem__('{}'.format(element),stellar_abund)
#sme.abund.__setitem__('Ca',sme.abund.__getitem__('Ca')-enh*2)
#print(sme.abund.__getitem__('{}'.format(element)))s
print(sme.abund)

#sme.h2broad=False

# Stellar atmosphere model is temporary until we can get MARCS 2014 working in PySME. Microturbulence and atmosphere geometery is different for the dwarf/giant model atmosphere grids
if logg >= 4.0:
    sme.atmo.geom = "PP"
else:
    sme.atmo.geom = "SPH"
 
sme.atmo.method = "grid"
sme.atmo.source = "marcs2014.sav"

# Finds the linelist file, the naming is slightly messy becaus of the 
#LLname = glob.glob(line_dir + '{}{}*.lin'.format(element,ion))
LLname = glob.glob(line_dir + '{}{}_{}.lin'.format(element,ion,line))
try:
    sme.linelist = ValdFile(LLname[0])
except IndexError:
    raise SystemExit('Linelist for {}{}_{} not found'.format(element,ion,line))

# Broadening, radial velocity and continuum. All set to zero as we aren't synthesizing a real star
sme.vsini = 0
sme.vrad = 0
sme.vmac = 0
sme.vrad_flag = "none"
sme.cscale_flag = "none"
sme.cscale_type = "mask"

# Accuracy parameters, needs to be at least 10**-5 for most elements, some require 10**-6.
#sme.accwi = 10**-prec
#sme.accrt = 10**-prec
sme.accwi = 10**-7
sme.accrt = 10**-3

# NLTE, easier to find than linelists due to more consistent naming scheme
if NLTE == 'y':
    sme.nlte.set_nlte(element, grid = NLTE_dir + "nlte_{}_scatt_pysme.grd".format(element))
    sme.nlte.solar = 'grevesse2007'

# Start and endpoints of spectra. sme.wave sets the resolution, currently using 0.01 AA
#sme.wran = np.array([[3000.000,30000.000]])
#sme.wran = np.array([[float(line)-1000.000,float(line)+1000.000]])
#sme.wave = np.linspace(3000,30000,2700001)
#sme.wave = np.concatenate([np.flip(float(line)-np.logspace(-2,3,500)),float(line)+np.logspace(-2,3,500)])
#sme.wave = 5688.2+np.concatenate([-np.flip(np.logspace(-3,3,100)),[0.],np.logspace(-3,3,100)])
sme.wran=np.array([float(line)-100,float(line)+100])

# The action part
sme = synthesize_spectrum(sme)
#sme.save(out_file)
#np.save(npy_file,sme)
#print(sme.nlte.flags == True)

# Saves a compressed spectra, rounds of decimals and removes the continuum. Filenames from target
#compress(sme,5,target,work_dir)
smeEW = simps(1-sme.synth,sme.wave)*1e3
#print("\n EW {:10.3f} {} {} {} {}".format(smeEW,str(index).zfill(5),NLTE,element,line))
print("\n EW {} {} {} {} {}".format(smeEW,str(index).zfill(5),NLTE,element,line))
#print(smeEW)
#EWs = pd.DataFrame([smeEW],columns=[line])
#EWs.to_csv(f'{ew_file}',index=False)

