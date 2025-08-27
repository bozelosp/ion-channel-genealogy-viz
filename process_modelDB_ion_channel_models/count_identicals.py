import pickle
from collections import defaultdict, OrderedDict

# Load the updated dictionary
with open("mod_source_code_dict_stripped_comments_whitespaces.pkl", "rb") as f:
    mod_source_code_dict = pickle.load(f)

# Dictionary to hold counts of occurrences of stripped_comments_whitespaces_source_code
count_dict = defaultdict(int)

# Populate the count_dict with counts
for mod_name, mod_data in mod_source_code_dict.items():
    stripped_code = mod_data["stripped_comments_whitespaces_source_code"]
    count_dict[stripped_code] += 1

# Sort the count_dict by count in ascending order
sorted_count_dict = OrderedDict(sorted(count_dict.items(), key=lambda x: x[1]))

# Print the counts
for stripped_code, count in sorted_count_dict.items():
    print(f"Occurrences: {count}\nCode: {stripped_code[:10]}... (truncated for brevity)\n{'-'*80}")

# total number of unique stripped_comments_whitespaces_source_code
print(f"Total number of unique stripped_comments_whitespaces_source_code: {len(sorted_count_dict)}")