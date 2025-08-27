import os
import sys
import json
import pickle

from pprint import pprint

# Load the icg_comprehensive_dict.pkl
with open('icg_comprehensive_dict.pkl', 'rb') as f:
	icg_comprehensive_dict = pickle.load(f)

# # get the first item and pprint it
# for key, value in icg_comprehensive_dict.items():
#     pprint(value)
#     break

# {'ICG': True,
#  'Supermodel 1': True,
#  'Supermodel 2': True,
#  'Year': 2005,
#  'ion_class': 'KCa',
#  'levenshtein_ratio': None,
#  'mod_filename': '45539_kahp_slower.mod',
#  'mod_filepath': 'ICG-github/icg-channels-KCa/45539_kahp_slower.mod/45539_kahp_slower.mod',
#  'modelDB_dir': '45539',
#  'modelDB_mod_ID': None,
#  'pool_id': None,
#  'source_code': 'TITLE Potasium AHP (slower) type current for RD Traub, et al '
#                 '2005\n'
#                 '\n'
#                 '\t}\n'
#                 '\tbeta = 0.001\n'
#                 '}\n'
#                 '\n'
#                 'UNITSON\n',
#  'stripped_comments_source_code': 'UNITS { \n'
#                                   '\t(mV) = (millivolt) \n'
#                                   '\t(mA) = (milliamp) \n'
#                                   '\t(mM) = (milli/liter)\n'
#                                   '} \n'
#                                   '\tbeta = 0.001\n'
#                                   '}\n'
#                                   '\n'
#                                   'UNITSON',
#  'stripped_comments_whitespaces_source_code': "UNITS{(mV)=(mil0.01}beta=0.001}UNITSON"}

# Load the modelDB_comprehensive_dict.pkl
with open('modelDB_comprehensive_dict.pkl', 'rb') as f:
	modelDB_comprehensive_dict = pickle.load(f)

# # get the first item and pprint it
# for key, value in modelDB_comprehensive_dict.items():
#     pprint(value)
#     break

# 'ICG': True,
# 'ICG_entries': [{'icg_model': '9889_rand',
#                  'ratio': 99.38282489860694}],
#  'Supermodel 1': False,
#  'Supermodel 2': True,
#  'Year': 2021,
#  'ion_class': None,
#  'mod_filename': 'distr.mod',
#  'mod_filepath': 'modelDB_unzipped_models/266864/purkinje_pf_source_code/mod/distr.mod',
#  'modelDB_dir': '266864',
#  'pool_id': 1,
#  'source_code': 'TITLE ...just to store peak membrane voltage\n'
#                 ': M.Migliore June 2001\n'
#                 ': T Morse February 2010 added times of occurrence\n'
#                 '\ttmax=0\n'
#                 '\t}\n'
#                 '\t\n'
#                 '}\n',
#  'stripped_comments_source_code': 'UNITS {\n'
#                                   '\t(mA) = (milliamp)\n'
#                                   '\t(mV) = (millivolt)\n'
#                                   '\ttmax=0\n'
#                                   '\t}\n'
#                                   '\t\n'
#                                   '}',
#  'stripped_comments_whitespaces_source_code': 'UNITS{(mA)=(mil70tmax=0}}',
#  'unique_modelDB_mod_id': '266864_distr-ID-1'}

# Load the pools_of_identicals.pkl
with open('pools_of_identicals.pkl', 'rb') as f:
	pools_of_identicals = pickle.load(f)

# # get the first item and pprint it
# for key, value in pools_of_identicals.items():
#     print(key,value)
#     break

# 4471 ['183300_ampanmda-ID-3', '185864_ampanmda-ID-1', '183300_ampanmda-ID-1', '183300_ampanmda-ID-2']

# Initialize unified_comprehensive_dict
unified_comprehensive_dict = {}

# Iterate over pools_of_identicals
for pool_id, models in pools_of_identicals.items():
	
	model_details = []
	highest_ratio_icg_entry = None
	
	# Iterate over each model in the pool
	for model in models:
		
		# Check if the model exists in modelDB_comprehensive_dict
		if model in modelDB_comprehensive_dict:

			# Get the model_info
			modelDB_info = modelDB_comprehensive_dict[model]
			
			# Remove unwanted keys from ModelDB info
			unwanted_keys_modelDB = ['source_code', 'stripped_comments_source_code', 
									 'stripped_comments_whitespaces_source_code', 'pool_id']
			for key in unwanted_keys_modelDB:
				modelDB_info.pop(key, None)

			icg_entries = modelDB_info.get('ICG_entries', None)

			enriched_icg_info = []

			if icg_entries:

				for icg_entry in icg_entries:
					icg_entry_name = icg_entry['icg_model']
					icg_model_info = icg_comprehensive_dict[icg_entry_name]
					
					# Remove unwanted keys from ICG info
					unwanted_keys_ICG = ['ion_class', 'Supermodel 1', 'Supermodel 2', 
										 'modelDB_dir', 'source_code', 'stripped_comments_source_code', 
										 'stripped_comments_whitespaces_source_code', 'ICG', 
										 'Year', 'pool_id']
					for key in unwanted_keys_ICG:
						icg_model_info.pop(key, None)

					icg_entry['info'] = icg_model_info

					enriched_icg_info.append(icg_entry)

					if (highest_ratio_icg_entry is None or 
						icg_entry['ratio'] > highest_ratio_icg_entry['ratio']):
						highest_ratio_icg_entry = icg_entry

				modelDB_info['ICG_entries'] = enriched_icg_info

		model_details.append(modelDB_info)

	# Sort models by year, accommodating None values
	sorted_by_year_model_details = sorted(model_details, key=lambda k: (k['Year'] is None, k['Year']))

	# Associate with pool_id
	unified_comprehensive_dict[pool_id] = {
		'modelDB_models': sorted_by_year_model_details,
		'ICG_model': highest_ratio_icg_entry['info'] if highest_ratio_icg_entry else None
	}

# Save unified_comprehensive_dict
with open('unified_comprehensive_dict.pkl', 'wb') as f:
	pickle.dump(unified_comprehensive_dict, f)

# Save unified_comprehensive_dict as json
with open('unified_comprehensive_dict.json', 'w') as f:
	json.dump(unified_comprehensive_dict, f, indent=4)