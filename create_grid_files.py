import pandas as pd
import os
element = "Ti"
stars = "grid"
NLTE =  ["n", "y"]
NLTE= ["n"]
ion= "1"

lines= ["4758", "4759", "4778", "4781", "4797", "4801", "4820", "5022", "5689", "5716", "5720", "5739", "5866", "6716", "7852"]
lines = ["4758", "4759", "5689", "5739"] # For Galah
lines = ["4759"]
if ion == "2":
    lines = ["4719", "4798", "4874"]
#lines = ["3635"]
#lines= [4719, 4764, 4798, 4849, 4865, 4874]
#home_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/klind/NaMgAl/'
import PySME_parallel_eqw
# Index of star with the given list of stellar parameters. Comes from     home_dir + stars +  '_params_parallel.txt'

indices = [42976, 42977, 42978, 42979, 42980, 42981, 42982, 42983, 42984, 42985, 42986, 42987, 42988, 42989, 42990, 42991, 42992, 42993, 42994, 42995, 42996, 42997, 42998, 42999, 43000, 43001, 43002, 43003, 43004, 43005, 43006, 43007, 43008, 43009, 43010, 43011, 43012, 43013, 43014, 43015, 43016, 43017, 43018, 43019, 43020, 43021, 43022, 43023, 43024, 43025, 43026, 43027, 43028, 43029, 43030, 43031, 43032, 43033, 43034, 43035, 43036, 43037, 43038, 43039, 43040, 43041, 43042, 43043, 43418, 43419, 43420, 43421, 43422, 43423, 43424, 43425, 43426, 43427, 43428, 43429, 43430, 43431, 43432, 43433, 43434, 43435, 43436, 43437, 43438, 43439, 43440, 43441, 43442, 43443, 43444, 43445, 43446, 43447, 43448, 43449, 43450, 43451, 43452, 43453, 43454, 43455, 43456, 43457, 43458, 43459, 43460, 43461, 43462, 43463, 43464, 43465, 43466, 43467, 43468, 43469, 43470, 43471, 43472, 43473, 43474, 43475, 43476, 43477, 43478, 43479, 43480, 43481, 43482, 43483, 43484, 43485, 51374, 51375, 51376, 51377, 51378, 51379, 51380, 51381, 51382, 51383, 51384, 51385, 51386, 51387, 51388, 51389, 51390, 51391, 51392, 51393, 51394, 51395, 51396, 51397, 51398, 51399, 51400, 51401, 51402, 51403, 51404, 51405, 51406, 51407, 51408, 51409, 51410, 51411, 51412, 51413, 51414, 51415, 51416, 51417, 51418, 51419, 51420, 51421, 51422, 51423, 51424, 51425, 51426, 51427, 51428, 51429, 51430, 51431, 51432, 51433, 51434, 51435, 51436, 51437, 51438, 51439, 51440, 51441, 51816, 51817, 51818, 51819, 51820, 51821, 51822, 51823, 51824, 51825, 51826, 51827, 51828, 51829, 51830, 51831, 51832, 51833, 51834, 51835, 51836, 51837, 51838, 51839, 51840, 51841, 51842, 51843, 51844, 51845, 51846, 51847, 51848, 51849, 51850, 51851, 51852, 51853, 51854, 51855, 51856, 51857, 51858, 51859, 51860, 51861, 51862, 51863, 51864, 51865, 51866, 51867, 51868, 51869, 51870, 51871, 51872, 51873, 51874, 51875, 51876, 51877, 51878, 51879, 51880, 51881, 51882, 51883]
indices = [8658, 8641]
#indices = list(range(0, 73372))
# odd departure grid coeffs of 1
#indices = [36913]
# working departure grids
#indices = [36947]
#indices = [23910]
#06191_4250_1.0_-5.0_1.0_-1.25_LTE
# Here we will collect EWs for all stars and lines involved.
ew_dict = {}
# The reason we index by line before star is that seems to be the EW file that we need to make. E.g one for each line
# that includes all stars examined

newrun = True
for line in lines:
    print("linehere", line)
    reset = True

    ew_dict[line] = {}
    for index in indices:
        print("index", index)
        ew_dict[line][index] = {}
        for TE in NLTE:


            home_dir = '/home/jama2357/Documents/TiFeAtoms/DepartureGrid/CoG/'
            work_dir = home_dir + element + "/" + element + ion + '/{}{}_{}/'.format(element,ion,line)
            if not os.path.exists(work_dir):
                os.makedirs(work_dir)
            if reset:
                #open(work_dir+"/NLTE_EW.txt", "w")
                #open(work_dir+"/LTE_EW.txt", "w")
                reset = False
            # produce the equivalent width of given line for given star (index)
            # Originally, karin seems to make a csv file for each equivalent width but we will collect them into one
            # dict for each star and line
            ew = PySME_parallel_eqw.ew(index, TE, ion, line, element, newrun)
            newrun = False
            if ew == -99:
                print("Incorrect EW from Atmosphere")
                continue
            if ew == 0:
                print("Check here create_grid_files")
                exit()
            ew_dict[line][index][TE] = ew
            print("line")
            print("ewdict", ew_dict)
            print("EW ", ew, index, TE, element, line)
            print("deets", "{:<5} {:<20} {:<10} {:<5} {:<5} {:<5}".format("EW", ew, index, TE, element, line))
            if TE == "y":
                print("{:<5} {:<20} {:<10} {:<5} {:<5} {:<5}".format("EW", ew, index, TE, element, line))

            elif TE == "n":
                print("{:<5} {:<20} {:<10} {:<5} {:<5} {:<5}".format("EW", ew, index, TE, element, line))
