import pickle
import json
import numpy as np
import networkx as nx
from networkx.readwrite import json_graph

# Function to calculate temporal distance based on the publication year
def temporal_distance_based_on_year(year_1, year_2):
	diff = abs(year_1 - year_2)
	return np.log(diff) if diff != 0 else 9

# Load unified comprehensive dict containing both ModelDB and ICG data
with open("unified_comprehensive_dict.pkl", "rb") as f:
	unified_comprehensive_dict = pickle.load(f)

# Load Levenshtein ratios
with open("levenshtein_ratios.pkl", "rb") as f:
	levenshtein_ratios = pickle.load(f)

# Initialize the graph
G = nx.Graph()

# Add nodes to the graph based on pools
for pool_id, pool_data in unified_comprehensive_dict.items():
	original_model = pool_data['modelDB_models'][0]  # The oldest model is the first one

	# node.original_model.ion_class if ion class is IH make it Ih
	if original_model['ion_class'] == 'IH':
		original_model['ion_class'] = 'Ih'

	identical_models = pool_data['modelDB_models'][1:]  # The rest are identical models
	ICG_model = pool_data['ICG_model']  # The representative ICG_model
	
	G.add_node(
		pool_id,
		original_model=original_model,
		identical_models=identical_models,
		num_of_identicals=len(identical_models),
		ICG_model=ICG_model,
		label=original_model['unique_modelDB_mod_id']
	)

# Add edges based on Levenshtein ratios (between pool IDs)
for (pool_1, pool_2), ratio in levenshtein_ratios.items():
	
	# If both pool IDs are in the graph and the ratio is greater than 75%
	if pool_1 in G.nodes and pool_2 in G.nodes and ratio >= 75:
		
		node_1 = G.nodes[pool_1]['original_model']
		node_2 = G.nodes[pool_2]['original_model']
		
		year_1 = node_1['Year']
		year_2 = node_2['Year']

		try:
			temporal_diff = temporal_distance_based_on_year(year_1, year_2)
		except:
			temporal_diff = 9
		
		# G.add_edge(
		# 	pool_1, pool_2,
		# 	weight=ratio,
		# 	temporal_diff=temporal_diff
		# )

		# we need to add an id too which should be str(pool_1) + '-' + str(pool_2)

		G.add_edge(
			pool_1, pool_2,
			weight=ratio,
			temporal_diff=temporal_diff,
			id=str(pool_1) + '-' + str(pool_2)
		)

# Save the graph to a JSON file
data = json_graph.node_link_data(G)
json.dump(data, open('network_data.json', 'w'))