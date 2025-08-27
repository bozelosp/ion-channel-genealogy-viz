import pickle
from collections import defaultdict

# Load the comprehensive dictionary
with open("icg_comprehensive_dict.pkl", "rb") as f:
    comprehensive_dict = pickle.load(f)

# Initialize a defaultdict to store ICG entries for each modelDB_mod_ID
modelDB_mod_ID_to_ICG_entries = defaultdict(list)

# Iterate through the dictionary to populate modelDB_mod_ID_to_ICG_entries
for key, value in comprehensive_dict.items():
    modelDB_mod_ID = value.get("modelDB_mod_ID", None)
    if modelDB_mod_ID is not None:
        modelDB_mod_ID_to_ICG_entries[modelDB_mod_ID].append(value)

# Save the new dictionary as a Pickle file
with open("modelDB_mod_ID_to_ICG_entries.pkl", "wb") as f:
    pickle.dump(modelDB_mod_ID_to_ICG_entries, f)

print("The data has been re-oriented by modelDB_mod_ID.")
