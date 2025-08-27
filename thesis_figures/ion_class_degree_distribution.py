import os
import sys
import json
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

def plot_degree_distribution(G, threshold, ax):
    degrees = [d for n, d in G.degree()]
    degree_counts = nx.degree_histogram(G)
    degrees = range(len(degree_counts))
    ax.scatter(degrees, degree_counts, s=10, marker='o')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Degree', fontsize=8)
    ax.set_ylabel('Count', fontsize=8)
    ax.set_title(f'Threshold: {threshold}%', fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=8)

    # Fit a power-law distribution
    x = np.array(degrees[1:])
    y = np.array(degree_counts[1:])
    mask = y > 0
    x = x[mask]
    y = y[mask]
    coefficients = np.polyfit(np.log(x), np.log(y), 1)
    exponent = coefficients[0]
    ax.plot(x, np.exp(coefficients[1]) * x**exponent, 'r-', label=f'Power-law (exponent: {exponent:.2f})')
    ax.legend(fontsize=6)

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

# Create subplots for each ion class and threshold
fig1, axes1 = plt.subplots(len(sorted_ion_classes), len(thresholds), figsize=(15, 5*len(sorted_ion_classes)))
fig1.subplots_adjust(hspace=0.6, wspace=0.6)

# Add supertitle at the top of the figure
supertitle = "Degree Distribution vs Similarity Thresholds"
fig1.suptitle(supertitle, fontsize=20, y=0.95, fontweight='bold')

# Determine the minimum and maximum values of degrees and degree counts for each ion class
ion_class_limits = {}
for ion_class in sorted_ion_classes:
    min_degree = float('inf')
    max_degree = -float('inf')
    min_count = float('inf')
    max_count = -float('inf')
    for threshold in thresholds:
        G_threshold = nx.Graph()
        for node_id in node_dict:
            if node_dict[node_id]['original_model']['ion_class'] == ion_class:
                G_threshold.add_node(node_id, **node_dict[node_id])
        for link in data['links']:
            if link['weight'] >= threshold and node_dict[link['source']]['original_model']['ion_class'] == ion_class and node_dict[link['target']]['original_model']['ion_class'] == ion_class:
                G_threshold.add_edge(link['source'], link['target'], **link)
        degrees = [d for n, d in G_threshold.degree()]
        degree_counts = nx.degree_histogram(G_threshold)
        min_degree = min(min_degree, min(degrees))
        max_degree = max(max_degree, max(degrees))
        min_count = min(min_count, min(degree_counts))
        max_count = max(max_count, max(degree_counts))
    ion_class_limits[ion_class] = (min_degree, max_degree, min_count, max_count)

panel_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']

for i, ion_class in enumerate(sorted_ion_classes):
    # Get the top position of the current row of subplots
    top = axes1[i, 0].get_position().y1

    # Add a centered title above each row of subplots with the number of nodes
    fig1.text(0.5, top + 0.0215, f'Ion Channel Class: {ion_class} — {ion_class_counts[ion_class]} nodes', fontsize=17, ha='center')

    # Add panel labels on the left side of each row, aligned with the ion class name text
    fig1.text(0.05, top + 0.0215, panel_labels[i], fontsize=17, ha='left', va='center', fontweight='bold')

    for j, threshold in enumerate(thresholds):
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

        ax_degree = axes1[i, j] if len(sorted_ion_classes) > 1 else axes1[j]
        plot_degree_distribution(G_threshold, threshold, ax_degree)
        ax_degree.set_xlim(left=ion_class_limits[ion_class][0], right=ion_class_limits[ion_class][1])
        ax_degree.set_ylim(bottom=ion_class_limits[ion_class][2], top=ion_class_limits[ion_class][3])

fname = 'ion_class_degree_distribution.png'
fig1.savefig(fname, dpi=600)

if sys.platform == 'darwin':
    os.system(f'open {fname}')
elif sys.platform == 'linux':
    os.system(f'xdg-open {fname}')