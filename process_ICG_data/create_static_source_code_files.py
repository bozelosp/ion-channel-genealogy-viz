import pickle
import os

# Portable absolute path
absolute_path = os.path.expanduser("~/Dropbox/icg-Chai-Panos/visualizer/")

# Create main directory and subdirectories for ICG
main_dir = os.path.join(absolute_path, "static/ICG/")
sub_dirs = ["source_code", "stripped_comments", "stripped_comments_whitespaces"]

# Create directories if they don't exist
if not os.path.exists(main_dir):
    os.makedirs(main_dir)

for sub_dir in sub_dirs:
    full_path = os.path.join(main_dir, sub_dir)
    if not os.path.exists(full_path):
        os.makedirs(full_path)

# Load the icg_comprehensive_dict.pkl dictionary
with open("icg_comprehensive_dict.pkl", "rb") as f:
    icg_comprehensive_dict = pickle.load(f)

# Iterate through the dictionary to save the code into appropriate directories
for model_name, model_data in icg_comprehensive_dict.items():
    
    mod_filename = model_data.get("mod_filename", None)

    if mod_filename is None:
        print(f"Could not find mod_filename for {model_name}")
        continue

    # Original source code
    source_code = model_data.get("source_code", None)
    if source_code:
        with open(os.path.join(main_dir, "source_code", f"{mod_filename}.mod"), "w") as f:
            f.write(source_code)

    # Stripped comments source code
    stripped_comments_source_code = model_data.get("stripped_comments_source_code", None)
    if stripped_comments_source_code:
        with open(os.path.join(main_dir, "stripped_comments", f"{mod_filename}.mod"), "w") as f:
            f.write(stripped_comments_source_code)

    # Stripped comments and whitespaces source code
    stripped_comments_whitespaces_source_code = model_data.get("stripped_comments_whitespaces_source_code", None)
    if stripped_comments_whitespaces_source_code:
        with open(os.path.join(main_dir, "stripped_comments_whitespaces", f"{mod_filename}.mod"), "w") as f:
            f.write(stripped_comments_whitespaces_source_code)

print("Source code files for ICG have been successfully saved.")