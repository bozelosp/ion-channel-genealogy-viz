import pickle

# Load the ICG super comprehensive dictionary
with open("ICG_super_comprehensive_dict.pkl", "rb") as f:
    icg_super_comprehensive_dict = pickle.load(f)

# Initialize a new dictionary to hold the swapped entries
swapped_dict = {}

# Loop through the original dictionary and swap the keys with the modelDB_model values
for original_key, original_value in icg_super_comprehensive_dict.items():
    new_key = original_value.get('modelDB_model', None)
    original_value['ICG_model_name'] = original_key

    # delete the ModelDB_mod_ID key
    del original_value['ModelDB_mod_ID']
    
    # Make sure the new key is not None before proceeding
    if new_key is not None:
        swapped_dict[new_key] = original_value

# Save the swapped dictionary back to disk
with open("swapped_ICG_super_comprehensive_dict.pkl", "wb") as f:
    pickle.dump(swapped_dict, f)