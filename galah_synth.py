import sys 
import os.path  
from pysme import sme as SME 
from pysme.abund import Abund 
from pysme.linelist.vald import ValdFile
from pysme.synthesize import synthesize_spectrum 
from pysme import util 
import numpy as np 
import glob 
import pandas as pd

home_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/NaMgAl/'
work_dir = home_dir + 'Sun'
line_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/Linelists/'
NLTE_dir = "/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/NLTE_grids/"

npy_file = os.path.join(work_dir, f"Sun.npy")
log_file = os.path.join(work_dir, f"Sun.log")
util.start_logging(log_file)

teff, logg, monh, vmic = 5777, 4.44, 0.00, 1.0

sme = SME.SME_Structure()

sme.teff, sme.logg, sme.monh, sme.vmic = teff, logg, monh, vmic

sme.abund = Abund(monh,"grevesse2007")
sme.atmo.geom = "PP"
sme.atmo.method = "grid"
sme.atmo.source = "marcs2014.sav"

LLname = glob.glob(line_dir + 'galah_master_v5.2.lin')
sme.linelist = ValdFile(LLname[0])

sme.vsini = 0
sme.vrad = 0
sme.vmac = 0
sme.vrad_flag = "none"
sme.cscale_flag = "none"
sme.cscale_type = "mask"

sme.accwi = 10**-7
sme.accrt = 10**-3

sme.nlte.set_nlte('Mg', grid = NLTE_dir + "nlte_Mg_scatt_pysme.grd")
sme.nlte.solar = 'grevesse2007'
sme.wran = np.array([[6550.000,6570.000]])

sme = synthesize_spectrum(sme)
np.save(npy_file,sme)
