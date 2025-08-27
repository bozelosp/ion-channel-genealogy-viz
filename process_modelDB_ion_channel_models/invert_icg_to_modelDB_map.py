import pickle
from pprint import pprint
import sys

# Load the ICG to ModelDB icg_to_modelDB_mapping_with_ratio with ratio
with open("icg_to_modelDB_mapping_with_ratio.pkl", "rb") as f:
	icg_to_modelDB_mapping_with_ratio = pickle.load(f)

# pprint(icg_to_modelDB_mapping_with_ratio)
#  '9889_gen': {'modelDB_model': '9889_gen-ID-1', 'ratio': 100},
#  '9889_nmda': {'modelDB_model': '9889_nmda-ID-1', 'ratio': 100},
#  '9889_passiv': {'modelDB_model': '9889_passiv-ID-1', 'ratio': 100},
#  '9889_presyn': {'modelDB_model': '9889_presyn-ID-1', 'ratio': 100},
#  '9889_pulse': {'modelDB_model': '9889_pulse-ID-1', 'ratio': 100},
#  '9889_rand': {'modelDB_model': '9889_rand-ID-1', 'ratio': 99.38282489860694}}

# Initialize an empty dictionary to store the inverted mapping
modelDB_to_icg_mapping_with_ratio = {}

# Iterate through the original icg_to_modelDB_mapping_with_ratio dictionary
for icg_model, modelDB_data in icg_to_modelDB_mapping_with_ratio.items():
    modelDB_model = modelDB_data["modelDB_model"]

    # Check if the ModelDB model already exists in the new dictionary
    if modelDB_model not in modelDB_to_icg_mapping_with_ratio:
        # If not, initialize an empty list
        modelDB_to_icg_mapping_with_ratio[modelDB_model] = []

    # Append the ICG model to the list corresponding to the ModelDB model
    modelDB_to_icg_mapping_with_ratio[modelDB_model].append({
        "icg_model": icg_model,
        "ratio": modelDB_data["ratio"]
    })

    # Sort the list by ratio in descending order
    modelDB_to_icg_mapping_with_ratio[modelDB_model].sort(key=lambda x: x["ratio"], reverse=True)

# Save the inverted dictionary as a pickle file
with open("modelDB_to_icg_mapping_with_ratio_sorted.pkl", "wb") as f:
    pickle.dump(modelDB_to_icg_mapping_with_ratio, f)