import pickle
import csv
import re
import random

# Step 1: Load the comprehensive ion_channel_model_to_class dictionary
with open("ion_channel_model_to_class_comprehensive.pkl", "rb") as f:
    ion_channel_model_to_class = pickle.load(f)

# Step 2: Load the model_data.csv and create a dictionary mapping model_id to year
model_to_year = {}
with open("model_data.csv", "r") as f:
    csv_reader = csv.reader(f)
    next(csv_reader)  # Skip the header row
    for row in csv_reader:
        model_to_year[row[1]] = row[3]  # Map model_id to year

# Step 3: Load the ICG to ModelDB mapping
with open("icg_to_modelDB_mapping.pkl", "rb") as f:
    icg_to_modelDB_mapping = pickle.load(f)

# remove .mod from the keys
icg_to_modelDB_mapping = {re.sub(r"\.mod$", "", k): v for k, v in icg_to_modelDB_mapping.items()}

# Initialize the new dictionary
ion_channel_model_to_class_ICG_supermodels = {}

# Regex pattern for extracting model_id
pattern = re.compile(r"^(\d+)")

# Step 4: Iterate over each model in the comprehensive dictionary
for model_name, model_class in ion_channel_model_to_class.items():

    # Initialize properties
    properties = {"Class": model_class, "ICG": True}

    # Step 5: Randomly assign Supermodel 1 and Supermodel 2 properties
    properties["Supermodel 1"] = random.random() < 0.7
    properties["Supermodel 2"] = random.random() < 0.7

    # Step 6: Use regex to extract model_id and look up the year
    match = pattern.match(model_name)
    if match:
        model_id = match.group(1)
        properties["Year"] = model_to_year.get(model_id, None)

    # Step 7: Add ModelDB mapping if it exists in icg_to_modelDB_mapping
    if model_name in icg_to_modelDB_mapping:
        properties["ModelDB_mod_ID"] = icg_to_modelDB_mapping[model_name]
    else:
        properties["ModelDB_mod_ID"] = None

    # Add the model with its properties to the new dictionary
    ion_channel_model_to_class_ICG_supermodels[model_name] = properties

# Step 8: Save the new dictionary as a Pickle file
with open("ion_channel_model_to_class_year_ICG_supermodels.pkl", "wb") as f:
    pickle.dump(ion_channel_model_to_class_ICG_supermodels, f)