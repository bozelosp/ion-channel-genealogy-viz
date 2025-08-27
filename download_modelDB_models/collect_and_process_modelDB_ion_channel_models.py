import os
import pickle
import json
import re

# Function to strip comments from NMODL source code
def strip_nmodl_comments(source_code):

	source_code = re.sub(r'^TITLE.*$', '', source_code, flags=re.MULTILINE)

	# Step 1: Identify VERBATIM blocks
	verbatim_blocks = re.findall(r'VERBATIM(.*?)ENDVERBATIM', source_code, flags=re.DOTALL)
	
	# Step 2: Replace VERBATIM blocks with placeholders
	placeholders = []
	for i, block in enumerate(verbatim_blocks):
		placeholder = f"__PLACEHOLDER__{i}__"
		placeholders.append(placeholder)
		source_code = re.sub(f"VERBATIM{re.escape(block)}ENDVERBATIM", placeholder, source_code, flags=re.DOTALL)
		
	# Step 3: Strip NMODL-specific comments
	# Remove single-line comments starting with :
	source_code = re.sub(r':.*$', '', source_code, flags=re.MULTILINE)
	# Remove multi-line comments (COMMENT ... ENDCOMMENT)
	source_code = re.sub(r'COMMENT.*?ENDCOMMENT', '', source_code, flags=re.DOTALL)
	
	# Step 4: Restore VERBATIM blocks and remove C-style comments within them
	for i, (block, placeholder) in enumerate(zip(verbatim_blocks, placeholders)):
		# Remove C-style single-line comments (// ...)
		block = re.sub(r'//.*$', '', block, flags=re.MULTILINE)
		# Remove C-style multi-line comments (/* ... */)
		block = re.sub(r'/\*.*?\*/', '', block, flags=re.DOTALL)
		# Restore the VERBATIM block
		source_code = source_code.replace(placeholder, f"VERBATIM{block}ENDVERBATIM")
	
	source_code = source_code.strip()

	return source_code

# [3rd script] Function to strip white spaces from NMODL source code
def strip_nmodl_whitespaces(source_code):
	return ''.join(source_code.split())


root_dir = "modelDB_unzipped_models"  # Replace with your actual root directory

mod_source_code_dict = {}

for dirpath, dirnames, filenames in os.walk(root_dir):

	# Skip hidden directories
	dirnames[:] = [d for d in dirnames if not d.startswith('.')]
	
	# Skip hidden files
	filenames = [f for f in filenames if not f.startswith('.')]
	
	for filename in filenames:
		if filename.endswith(".mod"):
			model_folder_name = dirpath.split("/")[1]

			filename_without_mod_extension = filename[:-4]
			
			# Create a unique identifier for each .mod file within a ModelDB directory
			unique_mod_id = 1
			while f"{model_folder_name}_{filename_without_mod_extension}-ID-{unique_mod_id}" in mod_source_code_dict:
				unique_mod_id += 1
			
			unambiguous_mod_name = f"{model_folder_name}_{filename_without_mod_extension}-ID-{unique_mod_id}"
			full_file_path = os.path.join(dirpath, filename)
			
			try:
				with open(full_file_path, 'r', encoding='ISO-8859-1') as f:
					source_code = f.read()

				# Strip comments
				stripped_comments_source_code = strip_nmodl_comments(source_code)

				# Strip whitespaces
				stripped_comments_whitespaces_source_code = strip_nmodl_whitespaces(stripped_comments_source_code)
				
				mod_source_code_dict[unambiguous_mod_name] = {
					"mod_filepath": full_file_path,
					"mod_filename": filename,
					"unique_modelDB_mod_id": unambiguous_mod_name,
					"modelDB_dir": model_folder_name,
					"source_code": source_code,
					"stripped_comments_source_code": stripped_comments_source_code,
					"stripped_comments_whitespaces_source_code": stripped_comments_whitespaces_source_code
				}
			except Exception as e:
				print(f"Could not read {unambiguous_mod_name} due to {e}")

# Save the dictionary
with open("modelDB_processed_mod_files.pkl", "wb") as f:
	pickle.dump(mod_source_code_dict, f)

# 11379 models were processed.