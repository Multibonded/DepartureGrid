"""find similar stars and check diff. between abundance/nlte effect and micro turb"""

import pickle as pkl
element = "Ti1"
galah     = pkl.load(open(element+"_Galah_impacts2_reduced.pkl", "rb"))
reliable2     = pkl.load(open(element+"_Galah_impacts_lineless2_reduced.pkl", "rb"))
reliable3     = pkl.load(open(element+"_NLTE_list2_reduced.pkl", "rb"))


print(galah[676216][4758])