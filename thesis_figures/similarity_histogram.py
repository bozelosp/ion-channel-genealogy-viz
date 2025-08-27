import json
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import math

# Load the JSON data
with open('network_data.json', 'r') as f:
    data = json.load(f)

# Create a dictionary to map node 'id' to the actual node item
node_dict = {node['id']: node for node in data['nodes']}

# Get unique ion classes
ion_classes = set(data['nodes'][link['source']]['original_model']['ion_class'] for link in data['links'])
ion_classes.discard(None)
ion_classes.discard('Other')

# Sort the ion classes alphabetically
sorted_ion_classes = sorted(ion_classes)

# Set up the plot using Seaborn
sns.set(style="whitegrid", font_scale=1.2)
num_rows = math.ceil(len(sorted_ion_classes) / 2)
fig, axs = plt.subplots(num_rows, 2, figsize=(16, 5 * num_rows), sharex=True)

# Define color palette
palette = sns.color_palette("pastel", len(sorted_ion_classes))

# Set the fixed threshold
threshold = 82.5

# Plot histograms and KDE of weights for each ion class
for i, ion_class in enumerate(sorted_ion_classes):
    weights = [link['weight'] for link in data['links'] if link['weight'] >= threshold and node_dict[link['source']]['original_model']['ion_class'] == ion_class and node_dict[link['target']]['original_model']['ion_class'] == ion_class]
    row = i // 2
    col = i % 2
    sns.histplot(weights, kde=True, color=palette[i], ax=axs[row, col])  # , kde_kws={'linewidth': 2})
    axs[row, col].set_title(f"{ion_class}", fontsize=16)
    axs[row, col].set_xlabel("Levenshtein Similarity (%)", fontsize=14, labelpad=15)
    axs[row, col].set_ylabel("Frequency", fontsize=14, labelpad=15)
    axs[row, col].tick_params(axis='both', labelsize=12)
    axs[row, col].grid(True, linestyle='--', alpha=0.7)
    
    # Add x-ticks
    xticks = list(range(int(min(weights)), int(max(weights)) + 1, 2))
    axs[row, col].set_xticks(xticks)
    axs[row, col].set_xticklabels(xticks)

# Remove empty subplots
if len(sorted_ion_classes) % 2 != 0:
    fig.delaxes(axs[num_rows - 1, 1])

# Adjust the spacing between subplots
fig.subplots_adjust(hspace=0.4, wspace=0.2)

# Add an overall title
fig.suptitle(f"Histogram of Levenshtein Similarity by Ion Channel Class (Threshold = {threshold}%)", fontsize=20, fontweight='bold', y=0.95)

fname = 'similarity_histogram.png'
fig.savefig(fname, dpi=300)

if sys.platform == 'darwin':
    os.system(f'open {fname}')
elif sys.platform == 'linux':
    os.system(f'xdg-open {fname}')

"""import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import math

# Load the JSON data
with open('network_data.json', 'r') as f:
    data = json.load(f)

# Create a dictionary to map node 'id' to the actual node item
node_dict = {node['id']: node for node in data['nodes']}

# Get unique ion classes
ion_classes = set(data['nodes'][link['source']]['original_model']['ion_class'] for link in data['links'])
ion_classes.discard(None)
ion_classes.discard('Other')

# Sort the ion classes alphabetically
sorted_ion_classes = sorted(ion_classes)

# Set up the plot using Seaborn
sns.set(style="whitegrid", font_scale=1.2)
num_rows = math.ceil(len(sorted_ion_classes) / 2)
fig, axs = plt.subplots(num_rows, 2, figsize=(16, 5 * num_rows), sharex=True)

# Define color palette
palette = sns.color_palette("pastel", len(sorted_ion_classes))

# Set the fixed threshold
threshold = 82.5

# Find the maximum y-value across all ion classes
max_y = 0
for ion_class in sorted_ion_classes:
    weights = [link['weight'] for link in data['links'] if link['weight'] >= threshold and node_dict[link['source']]['original_model']['ion_class'] == ion_class and node_dict[link['target']]['original_model']['ion_class'] == ion_class]
    hist, _ = np.histogram(weights, bins=range(int(min(weights)), int(max(weights)) + 2))
    max_y = max(max_y, max(hist)) + 32.5

# Plot histograms and KDE of weights for each ion class
for i, ion_class in enumerate(sorted_ion_classes):
    weights = [link['weight'] for link in data['links'] if link['weight'] >= threshold and node_dict[link['source']]['original_model']['ion_class'] == ion_class and node_dict[link['target']]['original_model']['ion_class'] == ion_class]
    row = i // 2
    col = i % 2
    sns.histplot(weights, kde=True, color=palette[i], ax=axs[row, col])
    axs[row, col].set_title(f"{ion_class}", fontsize=16)
    axs[row, col].set_xlabel("Levenshtein Similarity (%)", fontsize=14, labelpad=15)
    axs[row, col].set_ylabel("Frequency", fontsize=14, labelpad=15)
    axs[row, col].tick_params(axis='both', labelsize=12)
    axs[row, col].grid(True, linestyle='--', alpha=0.7)
    axs[row, col].set_ylim(top=max_y * 1.1)  # Set the y-axis maximum to max_y + 10%
    
    # Add x-ticks
    xticks = list(range(int(min(weights)), int(max(weights)) + 1, 2))
    axs[row, col].set_xticks(xticks)
    axs[row, col].set_xticklabels(xticks)

# Remove empty subplots
if len(sorted_ion_classes) % 2 != 0:
    fig.delaxes(axs[num_rows - 1, 1])

# Adjust the spacing between subplots
fig.subplots_adjust(hspace=0.4, wspace=0.2)

# Add an overall title
fig.suptitle(f"Histogram of Levenshtein Similarity by Ion Channel Class (Threshold = {threshold}%)", fontsize=20, fontweight='bold', y=0.95)

fname = 'similarity_histogram.png'
fig.savefig(fname, dpi=300)

if sys.platform == 'darwin':
    os.system(f'open {fname}')
elif sys.platform == 'linux':
    os.system(f'xdg-open {fname}')"""