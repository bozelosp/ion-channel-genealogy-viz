import json
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

def plot_clustering_coefficient(thresholds, clustering_coefficients, filename):
    plt.figure(figsize=(6, 4))
    plt.plot(thresholds, clustering_coefficients, marker='o')
    plt.xlabel('Similarity Threshold')
    plt.ylabel('Average Clustering Coefficient')
    plt.title('Average Clustering Coefficient vs. Similarity Threshold')
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def plot_network_density(thresholds, network_densities, filename):
    plt.figure(figsize=(6, 4))
    plt.plot(thresholds, network_densities, marker='o')
    plt.xlabel('Similarity Threshold')
    plt.ylabel('Network Density')
    plt.title('Network Density vs. Similarity Threshold')
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def plot_shortest_path_length(thresholds, shortest_path_lengths, filename):
    plt.figure(figsize=(6, 4))
    plt.plot(thresholds, shortest_path_lengths, marker='o')
    plt.xlabel('Similarity Threshold')
    plt.ylabel('Average Shortest Path Length')
    plt.title('Average Shortest Path Length vs. Similarity Threshold')
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def plot_community_detection(G, threshold, filename):
    communities = nx.community.louvain_communities(G)
    num_communities = len(communities)
    avg_community_size = np.mean([len(c) for c in communities])
    
    pos = nx.spring_layout(G)
    plt.figure(figsize=(8, 6))
    
    # Create a dictionary to map node IDs to indices
    node_index = {node: i for i, node in enumerate(G.nodes)}
    
    # Assign colors to nodes based on their community membership
    node_colors = np.zeros(len(G.nodes))
    for i, community in enumerate(communities):
        for node in community:
            node_colors[node_index[node]] = i
    
    nx.draw_networkx_nodes(G, pos, node_size=50, cmap=plt.cm.viridis, node_color=node_colors)
    nx.draw_networkx_edges(G, pos, alpha=0.3)
    plt.title(f'Community Detection (Threshold: {threshold}%)\nNumber of Communities: {num_communities}, Average Community Size: {avg_community_size:.2f}')
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def plot_centrality_distribution(G, threshold, filename):
    betweenness_centrality = nx.betweenness_centrality(G)
    closeness_centrality = nx.closeness_centrality(G)
    eigenvector_centrality = nx.eigenvector_centrality(G)
    
    plt.figure(figsize=(8, 6))
    plt.hist(list(betweenness_centrality.values()), bins=20, alpha=0.5, label='Betweenness Centrality')
    plt.hist(list(closeness_centrality.values()), bins=20, alpha=0.5, label='Closeness Centrality')
    plt.hist(list(eigenvector_centrality.values()), bins=20, alpha=0.5, label='Eigenvector Centrality')
    plt.xlabel('Centrality Score')
    plt.ylabel('Frequency')
    plt.title(f'Centrality Distribution (Threshold: {threshold}%)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def plot_assortativity(thresholds, assortativity_coefficients, filename):
    plt.figure(figsize=(6, 4))
    plt.plot(thresholds, assortativity_coefficients, marker='o')
    plt.xlabel('Similarity Threshold')
    plt.ylabel('Assortativity Coefficient')
    plt.title('Assortativity Coefficient vs. Similarity Threshold')
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# Load the JSON data
with open('network_data.json', 'r') as f:
    data = json.load(f)

# Create a NetworkX graph
G = nx.Graph()

# Add nodes to the graph
for node in data['nodes']:
    G.add_node(node['id'], **node)

# Define similarity thresholds
thresholds = [95, 90, 85]

# Plot clustering coefficient
clustering_coefficients = []
for threshold in thresholds:
    G_threshold = G.copy()
    for link in data['links']:
        if link['weight'] >= threshold:
            G_threshold.add_edge(link['source'], link['target'], **link)
    clustering_coefficient = nx.average_clustering(G_threshold)
    clustering_coefficients.append(clustering_coefficient)
plot_clustering_coefficient(thresholds, clustering_coefficients, 'clustering_coefficient.png')

# Plot network density
network_densities = []
for threshold in thresholds:
    G_threshold = G.copy()
    for link in data['links']:
        if link['weight'] >= threshold:
            G_threshold.add_edge(link['source'], link['target'], **link)
    network_density = nx.density(G_threshold)
    network_densities.append(network_density)
plot_network_density(thresholds, network_densities, 'network_density.png')

# Plot shortest path length
shortest_path_lengths = []
for threshold in thresholds:
    G_threshold = G.copy()
    for link in data['links']:
        if link['weight'] >= threshold:
            G_threshold.add_edge(link['source'], link['target'], **link)
    
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

plot_shortest_path_length(thresholds, shortest_path_lengths, 'shortest_path_length.png')

# Plot community detection
for threshold in thresholds:
    G_threshold = G.copy()
    for link in data['links']:
        if link['weight'] >= threshold:
            G_threshold.add_edge(link['source'], link['target'], **link)
    plot_community_detection(G_threshold, threshold, f'community_detection_{threshold}.png')

# Plot centrality distribution
for threshold in thresholds:
    G_threshold = G.copy()
    for link in data['links']:
        if link['weight'] >= threshold:
            G_threshold.add_edge(link['source'], link['target'], **link)
    plot_centrality_distribution(G_threshold, threshold, f'centrality_distribution_{threshold}.png')

# Plot assortativity
assortativity_coefficients = []
for threshold in thresholds:
    G_threshold = G.copy()
    for link in data['links']:
        if link['weight'] >= threshold:
            G_threshold.add_edge(link['source'], link['target'], **link)
    assortativity_coefficient = nx.degree_assortativity_coefficient(G_threshold)
    assortativity_coefficients.append(assortativity_coefficient)
plot_assortativity(thresholds, assortativity_coefficients, 'assortativity.png')