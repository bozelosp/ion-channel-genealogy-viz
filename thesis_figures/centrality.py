import json
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import numpy as np
import os
import sys

def plot_centrality_measure(G, thresholds, ax, ion_class, ion_class_counts, centrality_func, centrality_name, color):
    centrality_values = []
    for threshold in thresholds:
        G_threshold = G.copy()
        for link in data['links']:
            if link['weight'] >= threshold and data['nodes'][link['source']]['original_model']['ion_class'] == ion_class:
                G_threshold.add_edge(link['source'], link['target'], **link)
        try:
            centrality = centrality_func(G_threshold, max_iter=1000)
        except nx.PowerIterationFailedConvergence:
            centrality = {}
        centrality_values.append(np.mean(list(centrality.values())) if centrality else 0)
    ax.plot(thresholds, centrality_values, marker='o', label=f'{ion_class} ({ion_class_counts[ion_class]} nodes)', linewidth=2, color=color)

# Load the JSON data
with open('network_data.json', 'r') as f:
    data = json.load(f)

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

# Define centrality measures and their corresponding functions
centrality_measures = [
    ('Betweenness Centrality', nx.betweenness_centrality),
    ('Closeness Centrality', nx.closeness_centrality),
    ('Eigenvector Centrality', nx.eigenvector_centrality)
]

# Set up the plot using Seaborn
sns.set(style="whitegrid", font_scale=1.2)

# Define color palette
palette = sns.color_palette("pastel", len(sorted_ion_classes))

for centrality_name, centrality_func in centrality_measures:
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot lines for each ion class
    for i, ion_class in enumerate(sorted_ion_classes):
        plot_centrality_measure(G, thresholds, ax, ion_class, ion_class_counts, centrality_func, centrality_name, palette[i])
    
    ax.set_xlabel('Similarity Threshold', fontsize=14, labelpad=20)
    ax.set_ylabel(centrality_name, fontsize=14, labelpad=20)
    ax.set_xticks(thresholds[::2])
    ax.set_xticklabels(thresholds[::2])
    ax.tick_params(axis='both', labelsize=12)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend(title='Ion Channel Class', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0., fontsize=12)
    ax.set_title(f"{centrality_name} vs Similarity Thresholds", fontsize=18, fontweight='bold', pad=20)
    
    # Adjust the layout to make space for the legend and increase margins
    plt.tight_layout()
    plt.subplots_adjust(right=0.75, top=0.9, bottom=0.15)
    
    fname = f"{centrality_name.lower().replace(' ', '_')}.png"
    fig.savefig(fname, dpi=300)
    
    if sys.platform == 'darwin':
        os.system(f'open {fname}')
    elif sys.platform == 'linux':
        os.system(f'xdg-open {fname}')