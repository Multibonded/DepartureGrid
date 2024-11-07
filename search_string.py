import os
atomname = r"/home/jama2357/Documents/TiFeAtoms/Jack_Home/PYTHONATOM.fits"

#atomname = r"D:\PhD\TiFeAtoms\Jack\atom.Fe463_nq1541144_wH.fits"


user_input = "/home/jama2357/Documents/TiFeAtoms/DepartureGrid/CoG/"
#user_input = "/home/jama2357/Documents/TiFeAtoms/Jack_Home/"
directory = os.listdir(user_input)

searchstring = r"NLTE_PARAMETER_SPACES"
fileswith = []
for fname in directory:
    if os.path.isfile(user_input + os.sep + fname):

        # Full path
        if ".py" not in fname or ".pyc" in fname or "search_string" in fname:
            continue

        f = open(user_input + os.sep + fname, 'r')

        if searchstring in f.read():
            print('found string in file %s' % fname)
            fileswith.append(fname)
        else:
            print('string not found')
        f.close()

for folname in directory:
    if os.path.isdir(user_input + os.sep + folname):
        directory = os.listdir(user_input + os.sep + folname)
        for fname in directory:

            if os.path.isfile(user_input + os.sep + folname + os.sep + fname):

                # Full path
                if ".py" not in fname or ".pyc" in fname or "search_string" in fname:
                    continue
                f = open(user_input + os.sep + folname + os.sep + fname, 'r')

                if searchstring in f.read():
                    print('found string in file %s' % fname)
                    fileswith.append(folname+os.sep+fname)
                else:
                    print('string not found')
                f.close()
print("files with it in:", fileswith)
