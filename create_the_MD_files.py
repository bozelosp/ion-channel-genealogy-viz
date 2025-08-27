import os
import pickle

icg_dir = "ICG/"  # directory where pickle files are located
md_dir = "md_files/"  # directory to save markdown files

# get the list of pickle files in the directory
pickle_files = [f for f in os.listdir(icg_dir) if f.startswith('gvals')]

for file in pickle_files:
	# create a new subdirectory for each pickle file
	sub_dir = os.path.join(md_dir, file.replace(".pkl", ""))
	if not os.path.exists(sub_dir):
		os.makedirs(sub_dir)
	
	# load the pickle file
	with open(os.path.join(icg_dir, file), 'rb') as f:
		data = pickle.load(f, encoding='latin1')

	# create a markdown file for each entry in the new subdirectory
	for key, value in data.items():
		md_file = key.replace(".mod", ".md")
		with open(os.path.join(sub_dir, md_file), 'w') as f:
			f.write(f"# {key}\n\n")
			for sub_key, sub_value in value.items():
				if isinstance(sub_value, dict):
					f.write(f"- **{sub_key}**:\n")
					for k, v in sub_value.items():
						f.write(f"  - {k}: {v}\n")
				else:
					f.write(f"- **{sub_key}**: {sub_value}\n")