import json
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import numpy as np
import os
import sys
from community import best_partition

def plot_louvain_modularity(G, thresholds, node_dict):
    modularities = []
    num_communities = []
    for threshold in thresholds:
        G_threshold = nx.Graph()
        for node_id in node_dict:
            G_threshold.add_node(node_id, **node_dict[node_id])

        for link in data['links']:
            if link['weight'] >= threshold:
                G_threshold.add_edge(link['source'], link['target'], **link)

        num_components = nx.number_connected_components(G_threshold)
        num_edges = G_threshold.number_of_edges()
        num_nodes = G_threshold.number_of_nodes()
        print(f"Threshold: {threshold}, Components: {num_components}, Edges: {num_edges}, Nodes: {num_nodes}")

        partition = best_partition(G_threshold)
        communities = [set() for _ in range(max(partition.values()) + 1)]
        for node, community in partition.items():
            communities[community].add(node)

        modularities.append(nx.algorithms.community.modularity(G_threshold, communities))
        num_communities.append(len(communities))
    
    fig, ax1 = plt.subplots(figsize=(10, 8))
    
    color1 = 'royalblue'
    ax1.plot(thresholds, modularities, marker='o', linewidth=2, color=color1)
    ax1.set_xlabel('Similarity Threshold', fontsize=14, labelpad=20)
    ax1.set_ylabel('Modularity', fontsize=14, labelpad=20, color=color1)
    ax1.tick_params(axis='both', labelsize=12)
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    ax2 = ax1.twinx()
    color2 = 'darkorange'
    ax2.plot(thresholds, num_communities, marker='o', linewidth=2, color=color2)
    ax2.set_ylabel('Number of Communities', fontsize=14, labelpad=20, color=color2)
    ax2.tick_params(axis='y', labelcolor=color2, labelsize=12)
    
    # Set x-tick labels to display every other threshold value
    xticks = range(min(thresholds), max(thresholds) + 1, 2)
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xticks)
    
    # Adjust the spacing between the title and the plot
    plt.subplots_adjust(top=0.9)
    
    # Add a title with increased padding
    title = 'Louvain Modularity and Number of Communities\nvs Similarity Thresholds'
    fig.suptitle(title, fontsize=18, fontweight='bold')
    
    fname = 'louvain_modularity.png'
    fig.savefig(fname, dpi=300, bbox_inches='tight')
    
    if sys.platform == 'darwin':
        os.system(f'open {fname}')
    elif sys.platform == 'linux':
        os.system(f'xdg-open {fname}')

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

# Set up the plot using Seaborn
sns.set(style="whitegrid", font_scale=1.2)

plot_louvain_modularity(G, thresholds, node_dict)