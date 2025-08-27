import json
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import numpy as np
import os
import sys

def plot_shortest_path_length(G, thresholds, ax, ion_class, ion_class_counts, color, node_dict):
    shortest_path_lengths = []
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
        print(f"Threshold: {threshold}, Components: {num_components}, Edges: {num_edges}, Nodes: {num_nodes}")

        if not nx.is_connected(G_threshold):
            # If the graph is not connected, calculate the average shortest path length for each connected component
            shortest_path_length_sum = 0
            num_components = 0
            for component in nx.connected_components(G_threshold):
                subgraph = G_threshold.subgraph(component)
                shortest_path_length_sum += nx.average_shortest_path_length(subgraph)
                num_components += 1
            shortest_path_length = shortest_path_length_sum / num_components
        else:
            shortest_path_length = nx.average_shortest_path_length(G_threshold)

        shortest_path_lengths.append(shortest_path_length)

    ax.plot(thresholds, shortest_path_lengths, marker='o', label=f'{ion_class} ({ion_class_counts[ion_class]} nodes)', linewidth=2, color=color)

# Load the JSON data
with open('network_data.json', 'r') as f:
    data = json.load(f)

# Create a dictionary to map node 'id' to the actual node item
node_dict = {node['id']: node for node in data['nodes']}

# Create a NetworkX graph
G = nx.Graph()

# Add nodes to the graph
for node in data['nodes']:
    G.add_node(node['id'], **node)

# Define similarity thresholds
thresholds = list(range(99, 81, -1))

# Get unique ion classes
ion_classes = set(data['nodes'][link['source']]['original_model']['ion_class'] for link in data['links'])
ion_classes.discard(None)
ion_classes.discard('Other')

# Count the number of nodes for each ion class
ion_class_counts = {}
for node in data['nodes']:
    ion_class = node['original_model']['ion_class']
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
    plot_shortest_path_length(G, thresholds, ax, ion_class, ion_class_counts, palette[i], node_dict)

ax.set_xlabel('Similarity Threshold', fontsize=14, labelpad=20)
ax.set_ylabel('Shortest Path Length', fontsize=14, labelpad=20)
ax.set_xticks(thresholds[::2])
ax.set_xticklabels(thresholds[::2])
ax.tick_params(axis='both', labelsize=12)
ax.grid(True, linestyle='--', alpha=0.7)
ax.legend(title='Ion Channel Class', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0., fontsize=12)
ax.set_title("Shortest Path Length vs Similarity Thresholds", fontsize=18, fontweight='bold', pad=20)

# Adjust the layout to make space for the legend and increase margins
plt.tight_layout()
plt.subplots_adjust(right=0.75, top=0.9, bottom=0.15)

fname = 'shortest_path_length.png'
fig.savefig(fname, dpi=300)

if sys.platform == 'darwin':
    os.system(f'open {fname}')
elif sys.platform == 'linux':
    os.system(f'xdg-open {fname}')