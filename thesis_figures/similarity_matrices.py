import json
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import numpy as np

# Load the JSON data
with open('network_data.json', 'r') as f:
    data = json.load(f)

# Create a dictionary to map node 'id' to the actual node item
node_dict = {node['id']: node for node in data['nodes']}

# Create a new link_dict where source and target are sets
link_dict = {}
for link in data['links']:
    source = link['source']
    target = link['target']
    weight = link['weight']
    link_set = frozenset([source, target])
    if link_set in link_dict:
        link_dict[link_set]['weight'] = max(link_dict[link_set]['weight'], weight)
    else:
        link_dict[link_set] = {'source': source, 'target': target, 'weight': weight}

# Get unique ion classes
ion_classes = set(node_dict[node]['original_model']['ion_class'] for node in link_dict.values())
ion_classes.discard(None)
ion_classes.discard('Other')

# Sort the ion classes alphabetically
sorted_ion_classes = sorted(ion_classes)

# Set the similarity thresholds
thresholds = [99, 95, 91, 87, 83, 79]

# Set up the plot using Seaborn
sns.set(style="whitegrid", font_scale=1.2)
num_cols = 3
num_rows = len(thresholds) // num_cols + (len(thresholds) % num_cols > 0)
fig, axs = plt.subplots(num_rows, num_cols, figsize=(6 * num_cols, 6 * num_rows))

# Flatten the axs array for easier indexing
axs = axs.flatten()

# Create a list to store all similarity matrices
all_similarity_matrices = []

# Create similarity matrices for each threshold
for i, threshold in enumerate(thresholds):
    # Initialize the similarity matrix with zeros
    similarity_matrix = np.zeros((len(sorted_ion_classes), len(sorted_ion_classes)), dtype=int)
    
    # Fill the similarity matrix
    for link in link_dict.values():
        source_class = node_dict[link['source']]['original_model']['ion_class']
        target_class = node_dict[link['target']]['original_model']['ion_class']
        if source_class in sorted_ion_classes and target_class in sorted_ion_classes and link['weight'] >= threshold:
            source_index = sorted_ion_classes.index(source_class)
            target_index = sorted_ion_classes.index(target_class)
            if source_index != target_index:  # Exclude the diagonal
                similarity_matrix[source_index, target_index] += 1
                similarity_matrix[target_index, source_index] += 1  # Add the symmetric entry
    
    # Append the similarity matrix to the list
    all_similarity_matrices.append(similarity_matrix)

# Calculate the global min and max values across all similarity matrices
global_min = min(np.min(matrix) for matrix in all_similarity_matrices)
global_max = max(np.max(matrix) for matrix in all_similarity_matrices)

# Plot the heatmaps for each threshold
for i, threshold in enumerate(thresholds):
    ax = axs[i]
    sns.heatmap(all_similarity_matrices[i], annot=True, fmt='d', cmap='YlGnBu', vmin=global_min, vmax=global_max,
                cbar_kws={'label': 'Number of Similar Models'}, ax=ax, mask=np.eye(len(sorted_ion_classes), dtype=bool))
    ax.set_title(f"Similarity Threshold: {threshold}%")
    ax.set_xticklabels(sorted_ion_classes, rotation=45, ha='right')
    ax.set_yticklabels(sorted_ion_classes, rotation=0)
    ax.set_xlabel("Target Ion Channel Class")
    ax.set_ylabel("Source Ion Channel Class")

# Remove empty subplots
for j in range(i + 1, len(axs)):
    fig.delaxes(axs[j])

# Adjust the spacing between subplots
plt.tight_layout(pad=3.0)  # Increase the padding between subplots

# Add a super title to the figure
fig.suptitle("Cross-Ion Channel Class Similarity Links vs. Similarity Threshold", fontsize=16, fontweight='bold')

# Adjust the spacing at the top to accommodate the super title
fig.subplots_adjust(top=0.92)

fname = 'similarity_matrices.png'
fig.savefig(fname, dpi=300)

if sys.platform == 'darwin':
    os.system(f'open {fname}')
elif sys.platform == 'linux':
    os.system(f'xdg-open {fname}')