import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
import sys
import re

# Define the directory where you want to save the markdown files
output_directory = 'output_markdown_files'

# Create the output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Load the supermodel_data.pkl file
with open("supermodel_data.pkl", "rb") as f:
	supermodels = pickle.load(f)

# Function to create markdown content for each entry
def create_markdown_content(entry_name, entry_data, images_folder):
	print(entry_name)
	markdown_content = f"# {entry_name}\n\n"  # Start with a title
	plot_status = False
	for key, value in entry_data.items():

		# SM1_PARAMS_SS, SM2_PARAMS_SS, SM3_PARAMS_SS, ...
		is_it_sm_params_ss=re.match(r"SM\d_PARAMS_SS", key)

		# SM1_ERROR\d+_SS, SM2_ERROR\d+_SS, SM3_ERROR\d+_SS, ...
		is_it_sm_error_ss=re.match(r"SM\d_ERROR\d+_SS", key)

		# SM4_PARAMS_TAU, SM3_PARAMS_TAU, SM1_PARAMS_TAU, SM5_PARAMS_TAU
		is_it_sm_params_tau=re.match(r"SM\d_PARAMS_TAU", key)

		# SM4_ERROR\d+_TAU, SM3_ERROR\d+_TAU, SM1_ERROR\d+_TAU, SM5_ERROR\d+_TAU
		is_it_sm_error_tau=re.match(r"SM\d_ERROR\d+_TAU", key)

		if is_it_sm_params_ss or is_it_sm_error_ss:
			continue
		
		if is_it_sm_params_tau or is_it_sm_error_tau:
			continue

		if key in ["RATE_VALS_SS", "RATE_VALS_TAU"]:

			if plot_status == False:
				markdown_content += f"## RATE_VALS\n\n"
			else:
				continue
		else:
			markdown_content += f"## {key}\n\n"
		
		if key == "RATE_VALS_V":
			continue
		elif key in ["RATE_VALS_SS", "RATE_VALS_TAU"]:
			
			# for sub_key, sub_value in value.items():
			# 	markdown_content += f"- **{sub_key}**: {sub_value}\n"
			
			if entry_data.get('RATES', False):
				states = entry_data.get('STATES', [])
				for state in states:
					try:
						image_path = plot_and_save_graphs(state, entry_data, images_folder)
					except Exception as e:
						print("====================================")
						print(entry_name)
						print(f"Error plotting and saving graphs: {e}")
						print("====================================")
						continue
					markdown_content += f"![{state} Plot]({image_path})\n"
			plot_status = True
		elif isinstance(value, dict):
			for sub_key, sub_value in value.items():
				markdown_content += f"- **{sub_key}**: {sub_value}\n"
		elif isinstance(value, list):
			for item in value:
				markdown_content += f"- {item}\n"
		elif isinstance(value, bool) or value is None:
			markdown_content += f"- {value}\n"
		else:
			markdown_content += f"{value}\n"
		
		markdown_content += "\n"

	states = entry_data.get('STATES', [])

	sm_fit_bools = [entry_data.get(f"SM{i}_FIT", False) for i in range(1, 6)]

	if len(states) > 1 and all(sm_fit_bools):
		# print(entry_name)
		# print(f"Number of states: {len(states)}")
		pass
	else:
		return markdown_content

	sigmoid_mapping = {'SM1_PARAMS_SS': "sigmoid", 'SM2_PARAMS_SS': "non-sigmoid"}

	map_value_numbers_to_abcd = {1: 'a', 2: 'b', 3: 'c', 4: 'd', 5: 'e', 6:'f', 7:'g', 8:'h', 9:'i', 10:'j', 11:'k', 12:'l'}

	for state in states:

		# Create a header for the table
		markdown_content += f"## {state} State Parameters\n\n"
		
		state_sm1_params_ss = entry_data.get('SM1_PARAMS_SS', {}).get(state, [])
		state_sm2_params_ss = entry_data.get('SM2_PARAMS_SS', {}).get(state, [])

		# Find the maximum length to make sure we have a well-formed table
		max_len = max(len(state_sm1_params_ss), len(state_sm2_params_ss))

		# The header row
		markdown_content += "| Type | " + " | ".join(f"{map_value_numbers_to_abcd[i]}" for i in range(1, max_len + 1)) + " |\n"

		# The alignment row
		markdown_content += "| --- | " + " | ".join("---" for _ in range(max_len)) + " |\n"

		# The 'Sigmoid' row
		markdown_content += "| Sigmoid | " + " | ".join(str(value) for value in state_sm1_params_ss) + " |\n"

		# The 'Non-Sigmoid' row
		markdown_content += "| Non-Sigmoid | " + " | ".join(str(value) for value in state_sm2_params_ss) + " |\n"

		markdown_content += "\n"  # Add an extra newline for readability after the table

	map_value_numbers_to_abcd = {1: 'a', 2: 'b', 3: 'c', 4: 'd', 5: 'e', 6: 'f', 7: 'g', 8: 'h', 9: 'i', 10: 'j'}  # Extend this dictionary as needed

	for state in states:

		# Create a header for the table
		markdown_content += f"## {state} State Errors\n\n"
		
		# find all SM1_ERROR\d+_SS and SM2_ERROR\d+_SS keys
		sm1_errors = {key: entry_data[key][state] for key in entry_data.keys() if re.match(fr"SM1_ERROR\d+_SS", key)}
		sm2_errors = {key: entry_data[key][state] for key in entry_data.keys() if re.match(fr"SM2_ERROR\d+_SS", key)}

		# Find the maximum length to make sure we have a well-formed table
		max_len = max(len(sm1_errors), len(sm2_errors))

		# The header row
		markdown_content += "| Type | " + " | ".join(f"Error {i}" for i in range(1, max_len + 1)) + " |\n"

		# The alignment row
		markdown_content += "| --- | " + " | ".join("---" for _ in range(max_len)) + " |\n"

		# The 'Sigmoid' row
		markdown_content += "| Sigmoid | " + " | ".join(str(sm1_errors.get(f"SM1_ERROR{i}_SS", '')) for i in range(1, max_len + 1)) + " |\n"

		# The 'Non-Sigmoid' row
		markdown_content += "| Non-Sigmoid | " + " | ".join(str(sm2_errors.get(f"SM2_ERROR{i}_SS", '')) for i in range(1, max_len + 1)) + " |\n"

		markdown_content += "\n"  # Add an extra newline for readability after the table


	sm_params_to_tau_mapping = {
		"SM4_PARAMS_TAU": "tau 1",
		"SM3_PARAMS_TAU": "tau 2",
		"SM1_PARAMS_TAU": "tau 3",
		"SM5_PARAMS_TAU": "tau 4"
	}

	params_to_name_mapping = {
		1: "Vh",
		2: "A",
		3: "b1",
		4: "b2",
		5: "c1",
		6: "c2",
		7: "d1",
		8: "d2",
		9: "e1",
		10: "e2"
	}

	for state in states:
		# Create a header for the table
		markdown_content += f"## {state} Tau Data\n\n"
		
		# find all tau lists
		taus = {
			sm_params_to_tau_mapping[key]: entry_data.get(key, {}).get(state, []) 
			for key in sm_params_to_tau_mapping
		}

		# Find the maximum length to make sure we have a well-formed table
		max_len = max(len(tau) for tau in taus.values())

		# The header row
		markdown_content += "| Type | " + " | ".join(params_to_name_mapping[i] for i in range(1, max_len + 1)) + " |\n"

		# The alignment row
		markdown_content += "| --- | " + " | ".join("---" for _ in range(max_len)) + " |\n"

		# Rows for each tau
		for tau_name, tau_values in taus.items():
			markdown_content += f"| {tau_name} | " + " | ".join(str(value) for value in tau_values) + " |\n"

		markdown_content += "\n"  # Add an extra newline for readability after the table

	sm_params_to_tau_mapping = {
		"SM4": "tau 1",
		"SM3": "tau 2",
		"SM1": "tau 3",
		"SM5": "tau 4"
	}

	error_names_mapping = {
		1: "Error 1",
		2: "Error 2",
		3: "Error 3"
	}

	for state in states:
		# Create a header for the table
		markdown_content += f"## {state} Tau Error Data\n\n"
		
		# Create a dictionary to hold all tau errors
		tau_errors = {}
		for key, tau_name in sm_params_to_tau_mapping.items():
			tau_errors[tau_name] = [
				entry_data.get(f"{key}_ERROR{error_index}_TAU", {}).get(state, [])
				for error_index in range(1, 4)
			]

		#pprint(tau_errors)

		# Find the maximum length to make sure we have a well-formed table
		max_len = max(len(error_lists) for error_lists in tau_errors.values())

		# The header row
		markdown_content += "| Type | " + " | ".join(error_names_mapping[i] for i in range(1, max_len + 1)) + " |\n"

		# The alignment row
		markdown_content += "| --- | " + " | ".join("---" for _ in range(max_len)) + " |\n"

		# Rows for each tau
		for tau_name, error_list in tau_errors.items():    
			markdown_content += f"| {tau_name} | " + " | ".join(str(value) for value in error_list) + " |\n"

		markdown_content += "\n"  # Add an extra newline for readability after the table

	
	return markdown_content


def plot_and_save_graphs(state, entry_data, images_folder):

	rate_vals_v = entry_data.get('RATE_VALS_V', {}).get(state, {})
	# print(type(rate_vals_v), len(rate_vals_v)), 
	# it's <class 'numpy.ndarray'> 56
	# we want to check if it is empty
	if len(rate_vals_v) == 0:
		rate_vals_v = default_rate_vals_v

	fig, axs = plt.subplots(1, 2, figsize=(10, 5))
	fig.suptitle(f"{state} for different temperatures")

	for i, rate_key in enumerate(["RATE_VALS_SS", "RATE_VALS_TAU"]):
		ax = axs[i]
		ax.set_title(rate_key)
		ax.set_xlabel('Membrane Potential (mV)')
		ax.set_ylabel('Value')

		rate_values = entry_data.get(rate_key, {})
		for name, values in rate_values.items():
			if f"{state}_" in name:
				temp_suffix = name.split('_')[-1]
				ax.plot(rate_vals_v, values, label=f"{temp_suffix}°C")
				ax.legend()

	image_path = f"{images_folder}/{state}.png"
	fig.savefig(image_path)
	plt.close(fig)
	
	# Return the absolute path of the image
	return os.path.abspath(image_path)

default_rate_vals_v = [-100.,  -90.,  -80.,  -78.,  -76.,  -74.,  -72.,  -70.,  -68.,  -66.,  -64.,  -62.,
  -60.,  -58.,  -56.,  -54.,  -52.,  -50.,  -48.,  -46.,  -44.,  -42.,  -40.,  -38.,
  -36.,  -34.,  -32.,  -30.,  -28.,  -26.,  -24.,  -22.,  -20.,  -18.,  -16.,  -14.,
  -12.,  -10.,   -8.,   -6.,   -4.,   -2.,    0.,    2.,    4.,    6.,    8.,   10.,
   20.,   30.,   40.,   50.,   60.,   70.,   80.,  100.]

user_path = os.path.expanduser('~')

for ion_class, models in supermodels.items():
	ion_class_directory = os.path.join(output_directory, ion_class)
	os.makedirs(ion_class_directory, exist_ok=True)
	
	for model_name, model_data in models.items():
		model_folder = f"{ion_class_directory}/{model_name}"
		os.makedirs(model_folder, exist_ok=True)
		
		images_folder = f"{model_folder}/images"
		os.makedirs(images_folder, exist_ok=True)
		
		markdown_content = create_markdown_content(model_name, model_data, images_folder)
		file_name = f"{model_folder}/{model_name}.md"
		
		with open(file_name, 'w') as f:
			f.write(markdown_content)

print(f"Markdown files and graphs have been written to {output_directory}/")