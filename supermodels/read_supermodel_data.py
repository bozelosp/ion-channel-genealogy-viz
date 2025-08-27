import pickle
from pprint import pprint

# (base) ➜  supermodels lsd
# total 62152
#  2976 -rw-r--r--@ 1 pbozelos   1.5M 22 Sep 14:35 icg-channels-IH.pkl
# 23400 -rw-r--r--@ 1 pbozelos    11M 22 Sep 14:35 icg-channels-K.pkl
#  2768 -rw-r--r--@ 1 pbozelos   1.3M 22 Sep 14:35 icg-channels-KCa.pkl
# 21664 -rw-r--r--@ 1 pbozelos    11M 22 Sep 14:35 icg-channels-Na.pkl
# 11344 -rw-r--r--@ 1 pbozelos   5.5M 22 Sep 14:42 icg-channels-Ca.pkl

# Load all the icg-channels-*.pkl files
supermodels = {}

for ion_class in ["IH", "K", "KCa", "Na", "Ca"]:
    
    with open(f"icg-channels-{ion_class}.pkl", "rb") as f:
        supermodels[ion_class] = pickle.load(f)

# Save the supermodels dictionary
with open("supermodel_data.pkl", "wb") as f:
    pickle.dump(supermodels, f)