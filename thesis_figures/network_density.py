import os
import sys
import json
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

def plot_network_density(G, thresholds, ax, ion_classes, ion_class_counts, palette, node_dict):
    for i, ion_class in enumerate(ion_classes):
        network_densities = []
        for threshold in thresholds:
            G_threshold = nx.Graph()
            for node_id in node_dict:
                if node_dict[node_id]['original_model']['ion_class'] == ion_class:
                    G_threshold.add_node(node_id, **node_dict[node_id])

            for link in data['links']:
                if link['weight'] >= threshold and node_dict[link['source']]['original_model']['ion_class'] == ion_class and node_dict[link['target']]['original_model']['ion_class'] == ion_class:
                    G_threshold.add_edge(link['source'], link['target'], **link)

            network_densities.append(nx.density(G_threshold))
            num_components = nx.number_connected_components(G_threshold)
            num_edges = G_threshold.number_of_edges()
            num_nodes = G_threshold.number_of_nodes()
            print(f"Threshold: {threshold}, Components: {num_components}, Edges: {num_edges}, Nodes: {num_nodes}")

        ax.plot(thresholds, network_densities, marker='o', label=f'{ion_class} ({ion_class_counts[ion_class]} nodes)', linewidth=2, color=palette[i])

    ax.set_xlabel('Similarity Threshold', fontsize=14, labelpad=20)
    ax.set_ylabel('Network Density', fontsize=14, labelpad=20)
    ax.set_xticks(thresholds[::2])
    ax.set_xticklabels(thresholds[::2])
    ax.tick_params(axis='both', labelsize=12)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend(title='Ion Channel Class', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0., fontsize=12)
    ax.set_title("Network Density vs Similarity Thresholds", fontsize=18, fontweight='bold', pad=20)

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
plot_network_density(G, thresholds, ax, sorted_ion_classes, ion_class_counts, palette, node_dict)

# Adjust the layout to make space for the legend and increase margins
plt.tight_layout()
plt.subplots_adjust(right=0.75, top=0.9, bottom=0.15)

fname = 'network_density.png'
fig.savefig(fname, dpi=300)

if sys.platform == 'darwin':
    os.system(f'open {fname}')
elif sys.platform == 'linux':
    os.system(f'xdg-open {fname}')