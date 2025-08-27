import pickle
import os

# Portable absolute path
absolute_path = os.path.expanduser("~/Dropbox/icg-Chai-Panos/visualizer/")

# Create main directory and subdirectories
main_dir = os.path.join(absolute_path, "static/modelDB/")
sub_dirs = ["source_code", "stripped_comments", "stripped_comments_whitespaces"]

# Create directories if they don't exist
if not os.path.exists(main_dir):
    os.makedirs(main_dir)

for sub_dir in sub_dirs:
    full_path = os.path.join(main_dir, sub_dir)
    if not os.path.exists(full_path):
        os.makedirs(full_path)

# Load the modelDB_comprehensive_dict.pkl dictionary
with open("modelDB_comprehensive_dict.pkl", "rb") as f:
    modelDB_comprehensive_dict = pickle.load(f)

# Iterate through the dictionary to save the code into appropriate directories
for model_name, model_data in modelDB_comprehensive_dict.items():
    
    unique_modelDB_mod_id = model_data.get("unique_modelDB_mod_id", None)

    if unique_modelDB_mod_id is None:
        print(f"Could not find unique_modelDB_mod_id for {model_name}")
        continue

    # Original source code
    source_code = model_data.get("source_code", None)
    if source_code:
        with open(os.path.join(main_dir, "source_code", f"{unique_modelDB_mod_id}.mod"), "w") as f:
            f.write(source_code)

    # Stripped comments source code
    stripped_comments_source_code = model_data.get("stripped_comments_source_code", None)
    if stripped_comments_source_code:
        with open(os.path.join(main_dir, "stripped_comments", f"{unique_modelDB_mod_id}.mod"), "w") as f:
            f.write(stripped_comments_source_code)

    # Stripped comments and whitespaces source code
    stripped_comments_whitespaces_source_code = model_data.get("stripped_comments_whitespaces_source_code", None)
    if stripped_comments_whitespaces_source_code:
        with open(os.path.join(main_dir, "stripped_comments_whitespaces", f"{unique_modelDB_mod_id}.mod"), "w") as f:
            f.write(stripped_comments_whitespaces_source_code)

print("Source code files have been successfully saved.")




# import pickle
# import os
# from pprint import pprint
# import sys

# # # Directory where you want to save the source code files
# # output_dir = "static/modelDB/"

# # # Create the output directory if it doesn't exist
# # if not os.path.exists(output_dir):
# #     os.makedirs(output_dir)

# # Load the modelDB_comprehensive_dict.pkl dictionary
# with open("modelDB_comprehensive_dict.pkl", "rb") as f:
#     modelDB_comprehensive_dict = pickle.load(f)

# # Iterate through the dictionary to assess the properties and retrieve the code
# for model_name, model_data in modelDB_comprehensive_dict.items():

#     # pprint(model_data)
#     # print(type(model_data))
#     # <class 'dict'>

#     # print a list of the keys
#     # print(model_data.keys())
#     # dict_keys(['mod_filepath', 'mod_filename', 'unique_modelDB_mod_id', 'modelDB_dir', 'source_code', 'stripped_comments_source_code', 'stripped_comments_whitespaces_source_code', 'pool_id', 'ICG', 'ICG_entries', 'ion_class', 'Supermodel 1', 'Supermodel 2', 'Year'])

#     # pprint every key value except the code
#     for key, value in model_data.items():

#         if key not in ["source_code", "stripped_comments_source_code", "stripped_comments_whitespaces_source_code"]:
#             print(f"{key}: {value}")

#     sys.exit()

#     print(f"Processing {model_name}")

#     # Assess the properties
#     pool_id = model_data.get("pool_id", None)
#     is_icg = model_data.get("ICG", False)
#     icg_entries = model_data.get("ICG_entries", None)
#     ion_class = model_data.get("ion_class", None)
#     is_supermodel_1 = model_data.get("Supermodel 1", False)
#     is_supermodel_2 = model_data.get("Supermodel 2", False)
#     year = model_data.get("Year", None)

#     # Output the properties for this model
#     print(f"  Pool ID: {pool_id}")
#     print(f"  Is ICG: {is_icg}")
#     print(f"  ICG Entries: {icg_entries}")
#     print(f"  Ion Class: {ion_class}")
#     print(f"  Is Supermodel 1: {is_supermodel_1}")
#     print(f"  Is Supermodel 2: {is_supermodel_2}")
#     print(f"  Year: {year}")

#     # Retrieve the source code (assuming it's stored in the 'code' key)
#     source_code = model_data.get("code", None)

#     if source_code:
#         # Save the source code into a file in the output directory
#         with open(os.path.join(output_dir, f"{model_name}.txt"), "w") as code_file:
#             code_file.write(source_code)

# print("Processing complete.")
