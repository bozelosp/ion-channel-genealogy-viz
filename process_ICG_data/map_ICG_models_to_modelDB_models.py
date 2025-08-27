import re
import json
import pickle
from rapidfuzz import process, fuzz

# Load the processed ICG model data
with open("ICG_processed_mod_files.pkl", "rb") as f:
	icg_data = pickle.load(f)

# Load the processed ModelDB model data
with open("modelDB_processed_mod_files.pkl", "rb") as f:
	modelDB_data = pickle.load(f)

# Initialize a dictionary to hold the mappings
icg_to_modelDB_mapping_with_ratio = {}

# Loop through each ICG model
for icg_mod_filename, icg_mod_data in icg_data.items():

	mod_filename = icg_mod_data.get("mod_filename", "")
	mod_filename = mod_filename[:-4]
	icg_key = icg_mod_data.get("modelDB_dir", "")
	icg_stripped_comments_whitespaces_source_code = icg_mod_data.get("stripped_comments_whitespaces_source_code", "")

	# try:
	# 	icg_key = re.search(r"^(\d+)_", mod_filename).group(1)
	# except Exception as e:
	# 	print(f"Could not process {mod_filename} due to {e}")
	# 	continue

	found_exact_match = False

	# Loop through each ModelDB model
	modelDB_models_with_ICG_key = [modelDB_mod_key for modelDB_mod_key in modelDB_data if modelDB_data[modelDB_mod_key]["modelDB_dir"] == icg_key]
	
	for modelDB_mod_key in modelDB_models_with_ICG_key:
		
		modelDB_mod_data = modelDB_data[modelDB_mod_key]
		
		modelDB_stripped_comments_whitespaces_source_code = modelDB_mod_data.get("stripped_comments_whitespaces_source_code", "")

		# Check if the "stripped_comments_whitespaces" properties match
		if icg_stripped_comments_whitespaces_source_code == modelDB_stripped_comments_whitespaces_source_code:
			icg_to_modelDB_mapping_with_ratio[mod_filename] = {"modelDB_model": modelDB_mod_key, "ratio": 100}
			found_exact_match = True
			break  # No need to continue this loop if we've found a match
			
	if not found_exact_match:

		# If no exact match was found, find the closest match
		modelDB_stripped_texts = [modelDB_data[modelDB_mod_key].get("stripped_comments_whitespaces_source_code", "") for modelDB_mod_key in modelDB_models_with_ICG_key]

		closest_match = process.extractOne(icg_stripped_comments_whitespaces_source_code, modelDB_stripped_texts, scorer=fuzz.ratio)
		
		# print('Closest match:', closest_match)

		if closest_match:
			index_of_closest_match = modelDB_stripped_texts.index(closest_match[0])
			closest_modelDB_mod_key = modelDB_models_with_ICG_key[index_of_closest_match]
			icg_to_modelDB_mapping_with_ratio[mod_filename] = {

					"modelDB_model": closest_modelDB_mod_key, 
					"ratio": closest_match[1]

				}

# Save the mapping to a pickle file
with open("icg_to_modelDB_mapping_with_ratio.pkl", "wb") as f:
	pickle.dump(icg_to_modelDB_mapping_with_ratio, f)

# Save the mapping to a JSON file
with open("icg_to_modelDB_mapping_with_ratio.json", "w") as f:
	json.dump(icg_to_modelDB_mapping_with_ratio, f, indent=4)

print("Finished creating and saving the mapping.")