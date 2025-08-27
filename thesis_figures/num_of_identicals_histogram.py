import json
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import math

# Load the JSON data
with open('network_data.json', 'r') as f:
    data = json.load(f)

# Get unique ion classes
ion_classes = set(node['original_model']['ion_class'] for node in data['nodes'])
ion_classes.discard(None)
ion_classes.discard('Other')

# Sort the ion classes alphabetically
sorted_ion_classes = sorted(ion_classes)

# Find the maximum frequency across all ion channel classes
max_frequency = 0
for ion_class in sorted_ion_classes:
    num_of_identicals = [node['num_of_identicals'] for node in data['nodes'] if node['original_model']['ion_class'] == ion_class]
    max_frequency = max(max_frequency, max(plt.hist(num_of_identicals, bins=range(max(num_of_identicals) + 2))[0]))
    # print ion class and max frequency
    print(f"{ion_class}: {max_frequency}")

# Set up the plot using Seaborn
sns.set(style="whitegrid", font_scale=1.2)
num_rows = math.ceil(len(sorted_ion_classes) / 2)
fig, axs = plt.subplots(num_rows, 2, figsize=(16, 5 * num_rows), sharex=True)

# Define color palette
palette = sns.color_palette("pastel", len(sorted_ion_classes))

# Plot histograms for each ion class
for i, ion_class in enumerate(sorted_ion_classes):
    num_of_identicals = [node['num_of_identicals'] for node in data['nodes'] if node['original_model']['ion_class'] == ion_class]
    row = i // 2
    col = i % 2
    sns.histplot(num_of_identicals, kde=False, bins=range(max(num_of_identicals) + 2), color=palette[i], ax=axs[row, col])
    axs[row, col].set_title(f"{ion_class}", fontsize=16)
    axs[row, col].set_xlabel("Number of Identicals", fontsize=14, labelpad=15)
    axs[row, col].set_ylabel("Frequency", fontsize=14, labelpad=15)
    axs[row, col].tick_params(axis='both', labelsize=12)
    axs[row, col].grid(True, linestyle='--', alpha=0.7)
    axs[row, col].set_yscale('log')
    axs[row, col].set_ylim(top=max_frequency * 1.1)  # Set the y-axis limit to max_frequency + 10%

# Remove empty subplots
if len(sorted_ion_classes) % 2 != 0:
    fig.delaxes(axs[num_rows - 1, 1])

# Adjust the layout and spacing between subplots
fig.subplots_adjust(hspace=0.4, wspace=0.2)

# Add an overall title
fig.suptitle("Histogram of Number of Identicals by Ion Channel Class", fontsize=20, fontweight='bold', y=0.95)

fname = 'num_of_identicals_histogram.png'
fig.savefig(fname, dpi=300)

if sys.platform == 'darwin':
    os.system(f'open {fname}')
elif sys.platform == 'linux':
    os.system(f'xdg-open {fname}')