import pickle
from collections import defaultdict, Counter

# Load the comprehensive dictionary
with open("icg_comprehensive_dict.pkl", "rb") as f:
    comprehensive_dict = pickle.load(f)

# Initialize a Counter to count the occurrences of each modelDB_mod_ID
modelDB_mod_ID_counter = Counter()

# Initialize a defaultdict to store keys for each modelDB_mod_ID
modelDB_mod_ID_to_keys = defaultdict(list)

# Iterate through the dictionary to populate the Counter and defaultdict
for key, value in comprehensive_dict.items():
    modelDB_mod_ID = value.get("modelDB_mod_ID", None)
    if modelDB_mod_ID is not None:
        modelDB_mod_ID_counter[modelDB_mod_ID] += 1
        modelDB_mod_ID_to_keys[modelDB_mod_ID].append(key)

# Find and print the modelDB_mod_IDs that have a count of >= 2
for modelDB_mod_ID, count in modelDB_mod_ID_counter.items():
    if count >= 2:
        print(f"{modelDB_mod_ID}: {count}")
        print(f"Keys: {modelDB_mod_ID_to_keys[modelDB_mod_ID]}")
        print('-' * 40)