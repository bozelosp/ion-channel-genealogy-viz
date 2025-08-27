import sys
import os
import pickle
import re
import glob

# Function to strip comments from NMODL source code
def strip_nmodl_comments(my_source_code):

	my_source_code = re.sub(r'^TITLE.*$', '', my_source_code, flags=re.MULTILINE)

	# Step 1: Identify VERBATIM blocks
	verbatim_blocks = re.findall(r'VERBATIM(.*?)ENDVERBATIM', my_source_code, flags=re.DOTALL)
	
	# Step 2: Replace VERBATIM blocks with placeholders
	placeholders = []
	for i, block in enumerate(verbatim_blocks):
		placeholder = f"__PLACEHOLDER__{i}__"
		placeholders.append(placeholder)
		my_source_code = re.sub(f"VERBATIM{re.escape(block)}ENDVERBATIM", placeholder, my_source_code, flags=re.DOTALL)
		
	# Step 3: Strip NMODL-specific comments
	# Remove single-line comments starting with :
	my_source_code = re.sub(r':.*$', '', my_source_code, flags=re.MULTILINE)
	# Remove multi-line comments (COMMENT ... ENDCOMMENT)
	my_source_code = re.sub(r'COMMENT.*?ENDCOMMENT', '', my_source_code, flags=re.DOTALL)
	
	# Step 4: Restore VERBATIM blocks and remove C-style comments within them
	for i, (block, placeholder) in enumerate(zip(verbatim_blocks, placeholders)):
		# Remove C-style single-line comments (// ...)
		block = re.sub(r'//.*$', '', block, flags=re.MULTILINE)
		# Remove C-style multi-line comments (/* ... */)
		block = re.sub(r'/\*.*?\*/', '', block, flags=re.DOTALL)
		# Restore the VERBATIM block
		my_source_code = my_source_code.replace(placeholder, f"VERBATIM{block}ENDVERBATIM")
	
	my_source_code = my_source_code.strip()

	return my_source_code

def strip_nmodl_whitespaces(source_code):
	return ''.join(source_code.split())

# Initialize
root_dir = os.path.expanduser("ICG-github/")  # Root directory containing all ion channel classes

# Initialize dictionary to store model data
mod_source_code_dict = {}

# Loop through each class directory
for ion_class in os.listdir(root_dir):

	# (base) ➜  process_ICG_data lsd ICG-github/icg-channels-Ca/
	# icg-channels-Ca/     icg-channels-IH/     icg-channels-K/      icg-channels-KCa/    icg-channels-Na/     icg-channels-Other/

	ion_class_name = ion_class.split("-")[-1]

	ion_class_dir = os.path.join(root_dir, ion_class)
	
	# Skip non-directory files
	if not os.path.isdir(ion_class_dir):
		continue

	# Loop through each model within the class directory
	for model_name in os.listdir(ion_class_dir):

		modelDB_dir_ID=model_name.split("_")[0]

		model_dir = os.path.join(ion_class_dir, model_name)

		# Only interested in directories, skip files like README.md
		if os.path.isdir(model_dir):

			# Search for the first .mod file within the model directory
			mod_files = glob.glob(os.path.join(model_dir, "*.mod"))
			mod_files = [mod_file for mod_file in mod_files if "sm1_" not in mod_file and "sm2_" not in mod_file]
			if mod_files:

				mod_file_path = mod_files[0]
				
				# Extract source code
				try:
					with open(mod_file_path, 'r', encoding='ISO-8859-1') as f:
						source_code = f.read()
					
					# Strip comments and whitespaces
					stripped_comments_source_code = strip_nmodl_comments(source_code)
					stripped_comments_whitespaces_source_code = strip_nmodl_whitespaces(stripped_comments_source_code)
					
					# Extract subdir_name
					subdir_name = os.path.basename(os.path.dirname(mod_file_path))
					icg_file_name = subdir_name[:-4]

					# Extract mod_filename
					mod_filename = os.path.basename(mod_file_path)

					full_file_path = os.path.join(model_dir, mod_filename)
					
					# Add to dictionary
					mod_source_code_dict[icg_file_name] = {
						"mod_filename": mod_filename,
						"mod_filepath": full_file_path,
						"modelDB_dir": modelDB_dir_ID,
						"source_code": source_code,
						"stripped_comments_source_code": stripped_comments_source_code,
						"stripped_comments_whitespaces_source_code": stripped_comments_whitespaces_source_code,
						"ion_class": ion_class_name,
					}
				
				except Exception as e:
					print(f"Could not process {model_name} due to {e}")

# Save processed model data
with open("ICG_processed_mod_files.pkl", "wb") as f:
	pickle.dump(mod_source_code_dict, f)

print("Finished all tasks and saved the updated data in 'ICG_processed_mod_files.pkl'")
