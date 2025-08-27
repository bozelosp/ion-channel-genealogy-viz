import json
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import numpy as np
import os
import sys
from community import best_partition

def plot_louvain_modularity(G, thresholds):
    modularities = []
    num_communities = []
    
    for threshold in thresholds:
        G_threshold = nx.Graph()
        for link in data['links']:
            if link['weight'] >= threshold:
                G_threshold.add_edge(link['source'], link['target'], **link)
        
        partition = best_partition(G_threshold)
        communities = [set() for _ in range(max(partition.values()) + 1)]
        for node, community in partition.items():
            communities[community].add(node)
        
        modularities.append(nx.algorithms.community.modularity(G_threshold, communities))
        num_communities.append(len(communities))
    
    # Plot modularity
    fig1, ax1 = plt.subplots(figsize=(10, 8))
    color1 = 'royalblue'
    ax1.plot(thresholds, modularities, marker='o', linewidth=2, color=color1)
    ax1.set_xlabel('Similarity Threshold', fontsize=14, labelpad=20)
    ax1.set_ylabel('Modularity', fontsize=14, labelpad=20, color=color1)
    ax1.tick_params(axis='both', labelsize=12)
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    xticks = range(min(thresholds), max(thresholds) + 1, 2)
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xticks)
    
    plt.subplots_adjust(top=0.9)
    title1 = 'Louvain Modularity vs Similarity Thresholds'
    fig1.suptitle(title1, fontsize=18, fontweight='bold')
    
    fname1 = 'louvain_modularity.png'
    fig1.savefig(fname1, dpi=300, bbox_inches='tight')
    
    # Plot number of communities
    fig2, ax2 = plt.subplots(figsize=(10, 8))
    color2 = 'darkorange'
    ax2.plot(thresholds, num_communities, marker='o', linewidth=2, color=color2)
    ax2.set_xlabel('Similarity Threshold', fontsize=14, labelpad=20)
    ax2.set_ylabel('Number of Communities', fontsize=14, labelpad=20, color=color2)
    ax2.tick_params(axis='both', labelsize=12)
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticks)
    
    plt.subplots_adjust(top=0.9)
    title2 = 'Number of Communities vs Similarity Thresholds'
    fig2.suptitle(title2, fontsize=18, fontweight='bold')
    
    fname2 = 'louvain_num_communities.png'
    fig2.savefig(fname2, dpi=300, bbox_inches='tight')
    
    if sys.platform == 'darwin':
        os.system(f'open {fname1}')
        os.system(f'open {fname2}')
    elif sys.platform == 'linux':
        os.system(f'xdg-open {fname1}')
        os.system(f'xdg-open {fname2}')

# Load the JSON data
with open('network_data.json', 'r') as f:
    data = json.load(f)

# Create a NetworkX graph
G = nx.Graph()

# Add nodes to the graph
for node in data['nodes']:
    G.add_node(node['id'], **node)

# Define similarity thresholds
thresholds = list(range(99, 82, -1))

# Set up the plot using Seaborn
sns.set(style="whitegrid", font_scale=1.2)

plot_louvain_modularity(G, thresholds)