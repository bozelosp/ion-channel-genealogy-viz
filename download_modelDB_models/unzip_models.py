import zipfile
import os

# Directory containing the zipped model files
zipped_models_dir = "modelDB_zipped_models"

# Directory to save the unzipped models
unzipped_models_dir = "modelDB_unzipped_models"

# Create the directory if it doesn't exist
if not os.path.exists(unzipped_models_dir):
	os.mkdir(unzipped_models_dir)

# Loop through each file in the zipped models directory
for filename in os.listdir(zipped_models_dir):
	if filename.endswith(".zip"):
		# Full path to the zipped file
		zipped_file_path = os.path.join(zipped_models_dir, filename)
		
		# Model ID (obtained from the filename by removing the ".zip" part)
		model_id = filename[:-4]
		
		# Directory to unzip the model into (named after the model ID)
		unzip_dir = os.path.join(unzipped_models_dir, model_id)
		
		# Create the directory if it doesn't exist
		if not os.path.exists(unzip_dir):
			os.mkdir(unzip_dir)
		
		try:
			# Unzip the model
			with zipfile.ZipFile(zipped_file_path, 'r') as zip_ref:
				zip_ref.extractall(unzip_dir)
			#print(f"Unzipped {filename} into {unzip_dir}/")
		except zipfile.BadZipFile:
			print(f"Corrupted zip file detected: {filename}")