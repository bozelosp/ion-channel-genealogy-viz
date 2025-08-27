import os
import json
import pickle

# Initialize dictionaries to store ion channel models
class_to_ion_channel_model = {}
ion_channel_model_to_class = {}

# Specify the path to the root directory (replace with your actual path)
root_dir = "ICG-github"

# Loop through the directories and subdirectories
for dirpath, dirnames, filenames in os.walk(root_dir):
    # Skip the root directory
    if dirpath == root_dir:
        continue
    
    # Extract ion channel class from directory path
    ion_channel_class = os.path.basename(os.path.normpath(dirpath)).split('-')[-1]
    
    # Skip directories that are not ion channel classes
    if ion_channel_class not in ['Ca', 'Na', 'KCa', 'K', 'Ih', 'Other']:
        continue
    
    # Initialize list to store model names for this ion channel class
    if ion_channel_class not in class_to_ion_channel_model:
        class_to_ion_channel_model[ion_channel_class] = []
    
    # Loop through subdirectories and append model names to the list
    for sub_dir in dirnames:
        model_name = sub_dir.split('/')[-1]
        model_name = model_name.replace('.mod', '')  # Remove the .mod extension
        class_to_ion_channel_model[ion_channel_class].append(model_name)
        ion_channel_model_to_class[model_name] = ion_channel_class

# Save the dictionaries to Pickle files
with open("class_to_ion_channel_model.pkl", 'wb') as f:
    pickle.dump(class_to_ion_channel_model, f)

with open("ion_channel_model_to_class.pkl", 'wb') as f:
    pickle.dump(ion_channel_model_to_class, f)

print("The ion channel models and mappings have been saved in both JSON and Pickle formats.")