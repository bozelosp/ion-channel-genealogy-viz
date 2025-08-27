import os
import sys
import json
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import numpy as np

def calculate_exponent(G):
    degrees = [d for n, d in G.degree()]
    degree_counts = nx.degree_histogram(G)
    x = np.array(range(1, len(degree_counts)))
    y = np.array(degree_counts[1:])
    mask = y > 0
    x = x[mask]
    y = y[mask]
    coefficients = np.polyfit(np.log(x), np.log(y), 1)
    return coefficients[0]

# Load the JSON data
with open('network_data.json', 'r') as f:
    data = json.load(f)

# Create a NetworkX graph
G = nx.Graph()

# Create a dictionary to map node IDs to node items
node_dict = {node['id']: node for node in data['nodes']}

# Add nodes to the graph
for node_id in node_dict:
    G.add_node(node_id, **node_dict[node_id])

# Define similarity thresholds
thresholds = [99, 96, 93, 90, 87, 84]

# Get unique ion classes
ion_classes = set(node_dict[link['source']]['original_model']['ion_class'] for link in data['links'])
ion_classes.discard(None)
ion_classes.discard('Other')

# Count the number of nodes for each ion class
ion_class_counts = {}
for node_id in node_dict:
    ion_class = node_dict[node_id]['original_model']['ion_class']
    if ion_class in ion_classes:
        ion_class_counts[ion_class] = ion_class_counts.get(ion_class, 0) + 1

# Sort the ion classes based on the number of nodes in descending order
sorted_ion_classes = sorted(ion_classes, key=lambda x: ion_class_counts[x], reverse=True)

# Set up the plot using Seaborn
sns.set(style="whitegrid", font_scale=1.2)
fig, ax = plt.subplots(figsize=(10, 8))

# Define color palette
palette = sns.color_palette("pastel", len(sorted_ion_classes))

# Plot lines for each ion class
for i, ion_class in enumerate(sorted_ion_classes):
    exponents = []
    for threshold in thresholds:
        G_threshold = nx.Graph()
        for node_id in node_dict:
            if node_dict[node_id]['original_model']['ion_class'] == ion_class:
                G_threshold.add_node(node_id, **node_dict[node_id])
        for link in data['links']:
            if link['weight'] >= threshold and node_dict[link['source']]['original_model']['ion_class'] == ion_class and node_dict[link['target']]['original_model']['ion_class'] == ion_class:
                G_threshold.add_edge(link['source'], link['target'], **link)
        num_components = nx.number_connected_components(G_threshold)
        num_edges = G_threshold.number_of_edges()
        num_nodes = G_threshold.number_of_nodes()
        print(f"Ion Class: {ion_class}, Threshold: {threshold}, Components: {num_components}, Edges: {num_edges}, Nodes: {num_nodes}")
        exponent = calculate_exponent(G_threshold)
        exponents.append(exponent)
    ax.plot(thresholds, exponents, marker='o', label=f'{ion_class} ({ion_class_counts[ion_class]} nodes)', color=palette[i], linewidth=2)

ax.set_xlabel('Similarity Threshold', fontsize=14, labelpad=20)
ax.set_ylabel('Power-law Exponent', fontsize=14, labelpad=20)
ax.set_xticks(thresholds[::2])
ax.set_xticklabels(thresholds[::2])
ax.tick_params(axis='both', labelsize=12)
ax.grid(True, linestyle='--', alpha=0.7)
ax.legend(title='Ion Channel Class', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0., fontsize=12)
ax.set_title("Power-law Exponent vs Similarity Thresholds", fontsize=18, fontweight='bold', pad=20)

# Adjust the layout to make space for the legend and increase margins
plt.tight_layout()
plt.subplots_adjust(right=0.75, top=0.9, bottom=0.15)

fname = 'ion_class_exponent_threshold.png'
fig.savefig(fname, dpi=600)
if sys.platform == 'darwin':
    os.system(f'open {fname}')
elif sys.platform == 'linux':
    os.system(f'xdg-open {fname}')