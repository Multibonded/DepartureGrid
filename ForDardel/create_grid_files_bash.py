import pandas as pd
import os
element = "Ti"
stars = "grid"
NLTE =  ["n", "y"]
ion= "2"
#lines= ["4758", "4759", "4778", "4781", "4797", "4801", "4820", "5689", "5716", "5720", "5739", "5866", "6716", "7852"]
#home_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/NaMgAl/'

import PySME_parallel_eqw_bash
# Index of star with the given list of stellar parameters. Comes from     home_dir + stars +  '_params_parallel.txt'


# Here we will collect EWs for all stars and lines involved.
ew_dict = {}
# The reason we index by line before star is that seems to be the EW file that we need to make. E.g one for each line
# that includes all stars examined
windex = int(sys.argv[1])
line=  sys.argv[2]

index = int(sys.argv[1])

ew_dict[line] = {}
ew_dict[line][index] = {}
for TE in NLTE:


    home_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/Jack/Departures/DepartureGrid/CoG/'
    work_dir = home_dir + element + "/" + element + ion + '/{}{}_{}/'.format(element,ion,line)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    # produce the equivalent width of given line for given star (index)
    # Originally, karin seems to make a csv file for each equivalent width but we will collect them into one
    # dict for each star and line
    ew = PySME_parallel_eqw_bash.ew(index, TE, ion, line, element)
    if ew == -99:
        continue
    ew_dict[line][index][TE] = ew
    print("line")
    print("ewdict", ew_dict)
    print("EW ", ew, index, TE, element, line)
    print("deets", "{:<5} {:<20} {:<10} {:<5} {:<5} {:<5}".format("EW", ew, index, TE, element, line))
    if TE == "y":
        with open(work_dir+"/NLTE_EW.txt", "a") as f:
            print("{:<5} {:<20} {:<10} {:<5} {:<5} {:<5}".format("EW", ew, index, TE, element, line), file=f)

    elif TE == "n":
        with open(work_dir+"/LTE_EW.txt", "a") as f:
            print("{:<5} {:<20} {:<10} {:<5} {:<5} {:<5}".format("EW", ew, index, TE, element, line), file=f)
