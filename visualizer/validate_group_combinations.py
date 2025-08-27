from itertools import combinations, product

def get_combinations(elements):
    all_combinations = []
    for r in range(1, len(elements) + 1):
        for combo in combinations(elements, r):
            all_combinations.append(combo)
    return all_combinations

# Define active filters and ion classes
active_filters = ['Supermodel 1', 'Supermodel 2', 'ICG entry']
active_ion_classes = ['K', 'Na', 'Ca', 'Ih' ]

# Generate all combinations of active filters
filter_combinations = get_combinations(active_filters)

# Initialize the list to hold the final combinations
final_combinations = []

# If 'All' is active, generate combinations only from the active filters
all_active = False

if all_active:
    final_combinations = [('All',)]
else:
    # Combine each filter combination with each ion class
    for ion_class in active_ion_classes:
        for filter_combo in filter_combinations:
            final_combinations.append((ion_class,) + filter_combo)

# Print the final combinations
for combo in final_combinations:
    print(combo)