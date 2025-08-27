import sys
import pickle
import csv
import re
import random
from pprint import pprint

# Load the icg_processed_mod_files.pkl dictionary
with open("ICG_processed_mod_files.pkl", "rb") as f:
	icg_processed_mod_files = pickle.load(f)

# Load the class_to_model_name_to_ICG_detailed_page.pkl dictionary
with open("class_to_model_name_to_ICG_detailed_page.pkl", "rb") as f:
	class_to_model_name_to_ICG_detailed_page = pickle.load(f)

# Load model_data.csv and create a dictionary icg_to_modelDB_mapping_with_ratio model_id to year
modelDB_dir_to_year = {}
with open("model_data.csv", "r") as f:
	csv_reader = csv.reader(f)
	next(csv_reader)  # Skip the header row
	for row in csv_reader:
		modelDB_dir_to_year[row[1]] = row[3]  # Map model_id to year

# Load the ICG to ModelDB icg_to_modelDB_mapping_with_ratio with ratio
with open("icg_to_modelDB_mapping_with_ratio.pkl", "rb") as f:
	icg_to_modelDB_mapping_with_ratio = pickle.load(f)

# Load the pools of identicals
with open("pools_of_identicals.pkl", "rb") as f:
	pools_of_identicals = pickle.load(f)

# Create a reverse mapping from ModelDB_mod_ID to pool_id
model_to_pool = {}
for pool_id, models in pools_of_identicals.items():
	for model in models:
		model_to_pool[model] = pool_id

idx = 0
modelDB_dir_list = []
# Iterate over each model in icg_processed_mod_files to enrich it with additional properties
for model_name, model_data in icg_processed_mod_files.items():

	#print(f"Processing {model_name}")

	ion_class = model_data["ion_class"]
	modelDB_dir = model_data["modelDB_dir"]

	# Get the ICG detailed page URL for this model
	icg_friendly_model_name = re.sub(r"^(\d+)_", r"\1-", model_name)
	icg_detailed_page_url = class_to_model_name_to_ICG_detailed_page[ion_class].get(icg_friendly_model_name, None)
	model_data["ICG_detailed_page_url"] = icg_detailed_page_url

	# the github pages look like that https://github.com/icgenealogy/icg-channels-K/blob/master/3344_kadg.mod/3344_kadg.mod
	# so we need to transfrom this property 'mod_filepath': 'ICG-github/icg-channels-K/3344_kadg.mod/3344_kadg.mod'
	# and save it as 'ICG_github_url': ' ...

	mod_filepath = model_data["mod_filepath"]
	# replace the first ICG-github with https://github.com/icgenealogy
	github_url = mod_filepath.replace("ICG-github", "https://github.com/icgenealogy")

	# now /blob/master after https://github.com/icgenealogy/icg-channels-([a-zA-Z]+)/
	github_url = re.sub(r"https://github.com/icgenealogy/icg-channels-([a-zA-Z]+)/", r"https://github.com/icgenealogy/icg-channels-\1/blob/master/", github_url)

	model_data["ICG_github_url"] = github_url

	# Add ICG property, they are all ICG models, so set to True
	model_data["ICG"] = True

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
	if model_name in icg_to_modelDB_mapping_with_ratio:
		model_data["modelDB_mod_ID"] = icg_to_modelDB_mapping_with_ratio[model_name]["modelDB_model"]
		model_data["levenshtein_ratio"] = icg_to_modelDB_mapping_with_ratio[model_name]["ratio"]
		model_data["pool_id"] = model_to_pool.get(model_data["modelDB_mod_ID"], None)
	else:
		idx += 1
		modelDB_dir_list.append(modelDB_dir)
		print(f"Could not find {model_name} in icg_to_modelDB_mapping_with_ratio")
		model_data["modelDB_mod_ID"] = None
		model_data["levenshtein_ratio"] = None
		model_data["pool_id"] = None

print(f"Could not find {idx} models in icg_to_modelDB_mapping_with_ratio")

# Save the enriched dictionary back to a Pickle file
with open("icg_comprehensive_dict.pkl", "wb") as f:
	pickle.dump(icg_processed_mod_files, f)

print("The icg_processed_mod_files.pkl dictionary has been enriched with additional properties.")
