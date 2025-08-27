import pickle
from collections import defaultdict

# Load the updated dictionary
with open("modelDB_processed_mod_files.pkl", "rb") as f:
    mod_source_code_dict = pickle.load(f)

# Dictionary to group identical models based on stripped_comments_whitespaces_source_code
identical_pools = defaultdict(set)

# Populate the identical_pools dictionary
for mod_name, mod_data in mod_source_code_dict.items():
    stripped_code = mod_data["stripped_comments_whitespaces_source_code"]
    identical_pools[stripped_code].add(mod_name)

# Assign identifiers to each pool and print the sets
pool_id = 1
pools_of_identicals = {}
for stripped_code, pool in identical_pools.items():
    pools_of_identicals[pool_id] = list(pool)
    pool_id += 1

# Save the identified pools as a pickle file
with open("pools_of_identicals.pkl", "wb") as f:
    pickle.dump(pools_of_identicals, f)

print(f"Identified {len(pools_of_identicals)} pools of identical ion channel models.")


''' RATIONALE

Grouping identical ion channel models based on their fully stripped-down source code serves a crucial function in optimizing our subsequent Levenshtein distance calculations. By categorizing these models into distinct pools, we essentially create sets of ion channel models that are functionally identical, as evidenced by their identical source code after stripping away comments and whitespaces. This pre-emptive categorization allows us to significantly reduce the computational burden of the Levenshtein analysis. Specifically, models within the same pool would have a Levenshtein distance of zero, obviating the need for pairwise comparison within these groups. This is particularly advantageous given that the number of unique handshakes, or pairwise comparisons, scales quadratically with the number of unique stripped-down source codes. As we have identified around 5267 unique stripped-down source codes, the number of unique handshakes would be in the order of millions. By grouping identical models beforehand, we can substantially reduce this number, thereby accelerating the Levenshtein analysis. Moreover, these identified pools will serve as singular nodes in our forthcoming graphical representation, allowing us to focus on the meaningful differences between truly distinct ion channel models. This not only simplifies the graph but also makes it more informative, as we eliminate redundancy and highlight functional diversity. Therefore, this step is not just a computational optimization but also an analytical refinement, enabling a more efficient and insightful exploration of ion channel model genealogy.

'''