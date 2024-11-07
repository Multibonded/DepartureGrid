import sys 
import os.path  
# Use sys.path.insert to /src/ if .bash_profile is not modified with: export PYTHONPATH=/proj/snic2020-16-23/private/USERNAME/SME-master/src/:$PYTHONPATH
from pysme import sme as SME #Standard sme imports
from pysme.abund import Abund 
from pysme.linelist.vald import ValdFile
import pysme
print(pysme.__file__)
from pysme.synthesize import synthesize_spectrum #Use to just synthesize spectra without fitting, otherwise sme.solve
from pysme import util #Used for logging file

import numpy as np #Using numpy to set up sme.wran and sme.wave, not needed if working from input files
import glob #Using glob to identify the correct linelist, needed because the names are not formatted identically; not all linelists have hfs and I need to remember
import pandas as pd

#from eqwidth import eqwidth
# JACK   REMOVED THIS BECAUSE DOESN'T EXIST
#from sme_compressor import compress
from scipy.integrate import simps

def ew(windex, NLTE, ion, line, element):
    "______________________________________________________"
    # Arguments from bash script
    stars = "grid"
    #home_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/NaMgAl/'
    home_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/Jack/Departures/DepartureGrid/CoG/'
    #work_dir = home_dir + '{}{}_{}/'.format(element,ion,line)
    work_dir = home_dir + element + "/" + element + ion + '/{}{}_{}/'.format(element, ion, line)

    # Trying out this pandas method to use seq for parallel, if it's too slow (approx 0.06 sec) I'll just use cat instead
    # _params_parallel_filt are the parameters that don't need interpolation
    # Params are like the MARCS paramaters, available in large steps that we must interpolate between.
    params = pd.read_csv(home_dir + stars +  '_params_parallel.txt',delim_whitespace = True,names=['i','teff','logg','monh','vmic','enh'])
    # Index of the star we choose, temp, logg, monh, vmic, enh

    """xl = []
    for x in range(len(params.values)): 
        if 6000 >= int(params.values[x][1]) >= 5750 and 4 >= round(params.values[x][2], 1) >= 3.5 and -2 >= round(params.values[x][3], 2) >= -2.5: #and 2>= round(params.values[x][4], 2) >= 1.00 and (0.75== round(params.values[x][5], 2) or round(params.values[x][5], 2)  == -0.75):
            print(x, params.values[x])
            xl.append(x)
    print(xl)
    exit()"""
    """xl = []
    for x in range(len(params.values)):
        if int(params.values[x][1]) == 5750 and 3.50 == round(params.values[x][2], 1)  and -2.50 == round(params.values[x][3], 2) and 1.00== round(params.values[x][4], 2) and  round(params.values[x][5], 2)  == 0.5:
            print(x, params.values[x])
            xl.append(x)
    print(xl)
    exit()"""
    our_params = params.loc[params.i == windex].values[0]
    print("params", our_params)

    index, teff, logg, monh, vmic, enh = int(our_params[0]), int(our_params[1]), round(our_params[2],1), round(our_params[3],2), round(our_params[4],1), round(our_params[5],2)
    EleIon = element + ion

    # target specifies filenames, index used to keep track of progress
    target = "{}_{}_{}_{}_{}_{}".format(str(index).zfill(5),teff,logg,monh,vmic,enh)
    if NLTE == 'y':
        target = target + "_NLTE"
        print("\n Synthesizing {} in NLTE for parameters {}".format(element,target))
    else:
        target = target + "_LTE"
        print("\n Synthesizing {} in LTE for parameters {}".format(element,target))

    #line_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/Linelists/'
    #NLTE_dir = "/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/NLTE_grids/"
    line_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/Jack/Departures/DepartureGrid/CoG/linelists/'
    NLTE_dir = "/cfs/klemming/projects/snic/pdc-bus-2022-4/Jack/Departures/DepartureGrid/CoG/NLTE_grids/"

    out_file = os.path.join(work_dir, f"{target}.sme")
    npy_file = os.path.join(work_dir, f"{target}.npy")
    log_file = os.path.join(work_dir, f"{target}.log")
    ew_file = os.path.join(work_dir, f"{target}.csv")
    ew_dict = os.path.join(work_dir, f"{target}.pkl")
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
    #sme.h2broad=False

    # Stellar atmosphere model is temporary until we can get MARCS 2014 working in PySME. Microturbulence and atmosphere geometery is different for the dwarf/giant model atmosphere grids
    if logg >= 4.0:
        sme.atmo.geom = "PP"
    else:
        sme.atmo.geom = "SPH"
    if logg < 1:
        return 0
    sme.atmo.method = "grid"
    sme.atmo.source = "marcs2014.sav"

    # Finds the linelist file, the naming is slightly messy because of the
    #LLname = glob.glob(line_dir + '{}{}*.lin'.format(element,ion))
    LLname = glob.glob(line_dir + '{}{}_{}.lin'.format(element,ion,line))
    try:
        sme.linelist = ValdFile(LLname[0])
    except IndexError:
        raise SystemExit('Linelist for {}{}_{} not found'.format(element,ion,line))

    """sme.linelist['term_upper'] = "y1D*:      3d2.(1D).4s.4p.(1P*)"
    sme.linelist['term_lower'] = "1P:      3d3.(2P).5s"""


    """sme.linelist['species'] = "Ti 1"
    sme.linelist['wlcent'] = 1"""

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
        sme.nlte.set_nlte(element, grid = r"/cfs/klemming/projects/snic/pdc-bus-2022-4/Jack/Departures/DepartureGrid/pysme/nlte_Ti_ama51_Nov2022_pysme.grd")
        sme.nlte.solar = 'grevesse2007'



    # Start and endpoints of spectra. sme.wave sets the resolution, currently using 0.01 AA
    #sme.wran = np.array([[3000.000,30000.000]])
    #sme.wran = np.array([[float(line)-1000.000,float(line)+1000.000]])
    #sme.wave = np.linspace(3000,30000,2700001)
    #sme.wave = np.concatenate([np.flip(float(line)-np.logspace(-2,3,500)),float(line)+np.logspace(-2,3,500)])
    #sme.wave = 5688.2+np.concatenate([-np.flip(np.logspace(-3,3,100)),[0.],np.logspace(-3,3,100)])
    sme.wran=np.array([float(line)-100,float(line)+100])

    # The action part
    try:
        sme = synthesize_spectrum(sme)
    except pysme.atmosphere.atmosphere.AtmosphereError:
        with open(home_dir + element + "/" + element + ion + "/AtmosphereErrors.txt", "a") as f:
            print(our_params, file=f)

        return -99
    #sme.save(out_file)
    #np.save(npy_file,sme)

    # Saves a compressed spectra, rounds of decimals and removes the continuum. Filenames from target
    #compress(sme,5,target,work_dir)
    smeEW = simps(1-sme.synth,sme.wave)*1e3
    print("smeEW", smeEW)
    #print("\n EW {:10.3f} {} {} {} {}".format(smeEW,str(index).zfill(5),NLTE,element,line))
    #print("\n EW {} {} {} {} {}".format(smeEW,str(index).zfill(5),NLTE,element,line))
    EWs = pd.DataFrame([smeEW],columns=[line])
    #EWs.to_csv(f'{ew_file}',index=False)
    #import pickle
    #pickle.dump(smeEW, open(ew_dict, "wb"))
    return smeEW

#ew(10000, "y", "1", "5896")

