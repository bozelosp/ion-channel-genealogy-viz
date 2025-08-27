import json
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

def plot_degree_distribution(G, threshold, ax, ion_class):
    degrees = [d for n, d in G.degree()]
    degree_counts = nx.degree_histogram(G)
    degrees = range(len(degree_counts))
    ax.scatter(degrees, degree_counts, s=10, marker='o')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Degree')
    ax.set_ylabel('Count')
    ax.set_title(f'Degree Distribution (Threshold: {threshold}%, Ion Class: {ion_class})')

    # Fit a power-law distribution
    x = np.array(degrees[1:])
    y = np.array(degree_counts[1:])
    mask = y > 0
    x = x[mask]
    y = y[mask]
    coefficients = np.polyfit(np.log(x), np.log(y), 1)
    exponent = coefficients[0]
    ax.plot(x, np.exp(coefficients[1]) * x**exponent, 'r-', label=f'Power-law (exponent: {exponent:.2f})')
    ax.legend()

def plot_network_density(G, threshold, ax, ion_class):
    network_density = nx.density(G)
    ax.bar(0, network_density, width=0.5, color='b')
    ax.set_xlim([-0.5, 0.5])
    ax.set_xlabel('Network Density')
    ax.set_ylabel('Value')
    ax.set_title(f'Network Density (Threshold: {threshold}%, Ion Class: {ion_class})')
    ax.text(0, network_density, f'{network_density:.2f}', ha='center', va='bottom')

def plot_clustering_coefficient(G, threshold, ax, ion_class):
    clustering_coefficient = nx.average_clustering(G)
    ax.bar(0, clustering_coefficient, width=0.5, color='b')
    ax.set_xlim([-0.5, 0.5])
    ax.set_xlabel('Clustering Coefficient')
    ax.set_ylabel('Value')
    ax.set_title(f'Clustering Coefficient (Threshold: {threshold}%, Ion Class: {ion_class})')
    ax.text(0, clustering_coefficient, f'{clustering_coefficient:.2f}', ha='center', va='bottom')

def plot_shortest_path_length(G, threshold, ax, ion_class):
    if not nx.is_connected(G):
        # If the graph is not connected, calculate the average shortest path length for each connected component
        shortest_path_length_sum = 0
        num_components = 0
        for component in nx.connected_components(G):
            subgraph = G.subgraph(component)
            shortest_path_length_sum += nx.average_shortest_path_length(subgraph)
            num_components += 1
        shortest_path_length = shortest_path_length_sum / num_components
    else:
        shortest_path_length = nx.average_shortest_path_length(G)
    ax.bar(0, shortest_path_length, width=0.5, color='b')
    ax.set_xlim([-0.5, 0.5])
    ax.set_xlabel('Shortest Path Length')
    ax.set_ylabel('Value')
    ax.set_title(f'Shortest Path Length (Threshold: {threshold}%, Ion Class: {ion_class})')
    ax.text(0, shortest_path_length, f'{shortest_path_length:.2f}', ha='center', va='bottom')

def plot_centrality_distribution(G, threshold, ax, ion_class):
    betweenness_centrality = nx.betweenness_centrality(G)
    closeness_centrality = nx.closeness_centrality(G)
    try:
        eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=2000)
    except nx.PowerIterationFailedConvergence:
        eigenvector_centrality = {}

    ax.hist(list(betweenness_centrality.values()), bins=20, alpha=0.5, label='Betweenness Centrality')
    ax.hist(list(closeness_centrality.values()), bins=20, alpha=0.5, label='Closeness Centrality')
    if eigenvector_centrality:
        ax.hist(list(eigenvector_centrality.values()), bins=20, alpha=0.5, label='Eigenvector Centrality')
    ax.set_xlabel('Centrality Score')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Centrality Distribution (Threshold: {threshold}%, Ion Class: {ion_class})')
    ax.legend()

# Load the JSON data
with open('network_data.json', 'r') as f:
    data = json.load(f)

# Create a NetworkX graph
G = nx.Graph()

# Add nodes to the graph
for node in data['nodes']:
    G.add_node(node['id'], **node)

# Define similarity thresholds
thresholds = [95, 90, 85, 80, 75]

# Get unique ion classes
ion_classes = set(data['nodes'][link['source']]['original_model']['ion_class'] for link in data['links'])
ion_classes.discard(None)
ion_classes.discard('Other')

# print(ion_classes)
# {'KCa', 'Na', 'K', 'Ih', 'Ca'}

# Create subplots for each ion class and threshold
fig1, axes1 = plt.subplots(len(ion_classes), len(thresholds), figsize=(15, 5*len(ion_classes)))
fig2, axes2 = plt.subplots(len(ion_classes), len(thresholds), figsize=(15, 5*len(ion_classes)))
fig3, axes3 = plt.subplots(len(ion_classes), len(thresholds), figsize=(15, 5*len(ion_classes)))

for i, ion_class in enumerate(ion_classes):
    for j, threshold in enumerate(thresholds):
        G_threshold = G.copy()
        for link in data['links']:
            if link['weight'] >= threshold and data['nodes'][link['source']]['original_model']['ion_class'] == ion_class:
                G_threshold.add_edge(link['source'], link['target'], **link)

        ax_degree = axes1[i, j] if len(ion_classes) > 1 else axes1[j]
        ax_density = ax_degree.twinx()
        plot_degree_distribution(G_threshold, threshold, ax_degree, ion_class)
        plot_network_density(G_threshold, threshold, ax_density, ion_class)

        ax_clustering = axes2[i, j] if len(ion_classes) > 1 else axes2[j]
        ax_shortest_path = ax_clustering.twinx()
        plot_clustering_coefficient(G_threshold, threshold, ax_clustering, ion_class)
        plot_shortest_path_length(G_threshold, threshold, ax_shortest_path, ion_class)

        ax_centrality = axes3[i, j] if len(ion_classes) > 1 else axes3[j]
        plot_centrality_distribution(G_threshold, threshold, ax_centrality, ion_class)

plt.tight_layout()
# plt.show()

# actually save this figure
fig1.savefig('degree_distribution_and_network_density.png')
fig2.savefig('clustering_coefficient_and_shortest_path_length.png')
fig3.savefig('centrality_distribution.png')