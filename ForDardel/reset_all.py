
element = "Ti"
ion = "1"

line = int(sys.argv[1])



home_dir = '/cfs/klemming/projects/snic/pdc-bus-2022-4/Jack/Departures/DepartureGrid/CoG/'
work_dir = home_dir + element + "/" + element + ion + '/{}{}_{}/'.format(element, ion, line)

with open(work_dir + "/AtmosphereErrors.txt", "w") as f:
    pass
with open(work_dir + "/NLTE_EW.txt", "w") as g:
    pass
with open(work_dir + "/LTE_EW.txt", "w") as h:
    pass

