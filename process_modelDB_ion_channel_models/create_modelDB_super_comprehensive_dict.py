import sys
import pickle
import csv
import re
import random
from pprint import pprint

# Load the icg_processed_mod_files.pkl dictionary
with open("ICG_processed_mod_files.pkl", "rb") as f:
	icg_processed_mod_files = pickle.load(f)

# Load the modelDB_processed_mod_files.pkl dictionary
with open("modelDB_processed_mod_files.pkl", "rb") as f:
	modelDB_processed_mod_files = pickle.load(f)

# Read the pickle file back in
with open("modelDB_to_icg_mapping_with_ratio_sorted.pkl", "rb") as f:
	modelDB_to_icg_mapping_with_ratio_sorted = pickle.load(f)

# Load the pools of identicals
with open("pools_of_identicals.pkl", "rb") as f:
	pools_of_identicals = pickle.load(f)

# Load model_data.csv and create a dictionary icg_to_modelDB_mapping_with_ratio model_id to year
modelDB_dir_to_year = {}
with open("model_data.csv", "r") as f:
	csv_reader = csv.reader(f)
	next(csv_reader)  # Skip the header row
	for row in csv_reader:
		modelDB_dir_to_year[row[1]] = row[3]  # Map model_id to year

# Create a reverse mapping from ModelDB_mod_ID to pool_id
model_to_pool = {}
for pool_id, models in pools_of_identicals.items():
	for model in models:
		model_to_pool[model] = pool_id

modeldb_science_url = "https://modeldb.science/"
github_base_url = "https://github.com/ModelDBRepository/"

# Iterate over each model in modelDB_processed_mod_files to enrich it with additional properties
for model_name, model_data in modelDB_processed_mod_files.items():

	#print(f"Processing {model_name}")

	model_data["pool_id"] = model_to_pool.get(model_name, None)

	modelDB_dir = model_data.get("modelDB_dir")
	if modelDB_dir:
		model_data["modelDB.science_url"] = f"{modeldb_science_url}{modelDB_dir}"
		model_data["modelDB_github_url"] =f"{github_base_url}{modelDB_dir}"

	if model_name in modelDB_to_icg_mapping_with_ratio_sorted:
	
		model_data["ICG"] = True
		model_data["ICG_entries"] = modelDB_to_icg_mapping_with_ratio_sorted[model_name]

		corresponding_icg_models = [ entry["icg_model"] for entry in model_data["ICG_entries"] ]
		corresponding_ion_classes = [ icg_processed_mod_files[entry]["ion_class"] for entry in corresponding_icg_models ]
		ion_class = corresponding_ion_classes[0]
		model_data["ion_class"] = ion_class

	else:
		# print(f"Could not find {model_name} in icg_to_modelDB_mapping_with_ratio")
		model_data["ICG"] = False
		model_data["ICG_entries"] = None
		model_data["ion_class"] = None

	# Randomly assign Supermodel 1 and Supermodel 2 properties
	model_data["Supermodel 1"] = random.random() < 0.7
	model_data["Supermodel 2"] = random.random() < 0.7

	# Use regex to extract model_id and look up the year
	modelDB_dir = model_data["modelDB_dir"]
	try:
		model_data["Year"] = int(modelDB_dir_to_year.get(modelDB_dir, None))  # Convert to int
	except ValueError:
		model_data["Year"] = None

	# Add ModelDB icg_to_modelDB_mapping_with_ratio with ratio if it exists in icg_to_modelDB_mapping_with_ratio

# Save the enriched dictionary back to a Pickle file
with open("modelDB_comprehensive_dict.pkl", "wb") as f:
	pickle.dump(modelDB_processed_mod_files, f)

print("The modelDB_processed_mod_files.pkl dictionary has been enriched with additional properties.")