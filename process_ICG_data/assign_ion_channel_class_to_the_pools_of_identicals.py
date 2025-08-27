# Import required modules
import pickle
import re
import json

from pprint import pprint

# Step 1: Load the existing ion channel model-to-class and class-to-ion-channel-model mappings from ICG
with open("ion_channel_model_to_class.pkl", "rb") as f:
    model_to_class = pickle.load(f)

with open("class_to_ion_channel_model.pkl", "rb") as f:
    class_to_ion_channel_model = pickle.load(f)

# Step 2: Load the identified pools
with open("identified_pools.pkl", "rb") as f:
    identified_pools = pickle.load(f)

# Initialize the result dictionary and class counter
result = {}
class_counter = {}

# Step 3: Iterate over each pool to identify their classes
for pool_id, models in identified_pools.items():
    pool_dict = {"pool": models, "class": None}

    for model in models:
        clean_model_name = re.sub(r'\.mod-ID-\d+$', '', model)

        if clean_model_name in model_to_class:
            pool_dict["class"] = model_to_class[clean_model_name]
            break

    if pool_dict["class"] in class_counter:
        class_counter[pool_dict["class"]] += 1
    else:
        class_counter[pool_dict["class"]] = 1

    if pool_dict["class"] is not None:
        result[pool_id] = pool_dict

        # Update the model_to_class and class_to_ion_channel_model dictionaries
        for model in models:
            clean_model_name = re.sub(r'\.mod-ID-\d+$', '', model)
            model_to_class[clean_model_name] = pool_dict["class"]

            if pool_dict["class"] not in class_to_ion_channel_model:
                class_to_ion_channel_model[pool_dict["class"]] = []
            class_to_ion_channel_model[pool_dict["class"]].append(clean_model_name)

# Print the statistics
print("Statistics of ion channel classes:")
pprint(class_counter)

# Save the comprehensive mappings
with open("ion_channel_model_to_class_comprehensive.pkl", "wb") as f:
    pickle.dump(model_to_class, f)

with open("class_to_ion_channel_model_comprehensive.pkl", "wb") as f:
    pickle.dump(class_to_ion_channel_model, f)