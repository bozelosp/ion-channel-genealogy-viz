import sys
import json
import matplotlib.pyplot as plt
import seaborn as sns
from pprint import pprint

fname = "network_data.json"

# Load the JSON data
with open(fname, 'r') as f:
    data = json.load(f)

# Extract the weights (similarity scores) and ion classes from the data
weights = [link['weight'] for link in data['links']]
ion_classes = [data['nodes'][link['source']]['original_model']['ion_class'] for link in data['links']]

# Create a dictionary to store the weights for each ion class
class_weights = {}
for weight, ion_class in zip(weights, ion_classes):
    if ion_class is not None and ion_class != 'Other':
        if ion_class not in class_weights:
            class_weights[ion_class] = []
        class_weights[ion_class].append(weight)

# Set the Seaborn style to "darkgrid"
sns.set_style("darkgrid")

# Determine the number of subplots needed
num_classes = len(class_weights)
num_cols = 2
num_rows = (num_classes + num_cols - 1) // num_cols

# Create a figure and subplots for each ion class
fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 8))
axes = axes.flatten()

# Set the color palette
colors = sns.color_palette("husl", num_classes)

for i, (ion_class, weights) in enumerate(class_weights.items()):
    ax = axes[i]
    if len(weights) > 0:
        sns.histplot(weights, stat='density', ax=ax, color=colors[i], alpha=0.7, bins=20)
        sns.kdeplot(weights, ax=ax, color=colors[i], linewidth=2)
        ax.set_xlabel('Similarity Score')
        ax.set_ylabel('Density')
        ax.set_title(f'Ion Class: {ion_class}')
        ax.set_xlim(70, 105)  # Set the x-axis limits
        ax.set_xticks(range(70, 105, 5))  # Set the x-axis tick positions
        ax.grid(True, axis='y', linestyle='--', alpha=0.7)  # Add horizontal grid lines
    else:
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center')
        ax.set_xlabel('Similarity Score')
        ax.set_ylabel('Density')
        ax.set_title(f'Ion Class: {ion_class}')

# Remove unused subplots
for i in range(num_classes, len(axes)):
    fig.delaxes(axes[i])

# Adjust the spacing between subplots
plt.tight_layout(pad=2.0)

plt.show()