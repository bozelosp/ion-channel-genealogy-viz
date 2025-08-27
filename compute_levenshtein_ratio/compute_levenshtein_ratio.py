import pickle
from rapidfuzz import fuzz
from tqdm import tqdm

# Load the identified pools of identical .mod files
with open("pools_of_identicals.pkl", "rb") as f:
    pools_of_identicals = pickle.load(f)

# Load the dictionary with stripped-down source codes
with open("modelDB_processed_mod_files.pkl", "rb") as f:
    modelDB_processed_mod_files = pickle.load(f)

# Dictionary to hold Levenshtein distances
levenshtein_distances = {}

# Initialize tqdm progress bar
total_comparisons = len(pools_of_identicals) * (len(pools_of_identicals) - 1) // 2
progress_bar = tqdm(total=total_comparisons, desc="Calculating", dynamic_ncols=True)

# Loop through each pool to get a unique representative model
for pool_key, pool_models in pools_of_identicals.items():
    unique_model_a = pool_models[0]
    code_a = modelDB_processed_mod_files[unique_model_a]['stripped_comments_whitespaces_source_code']
    
    # Calculate Levenshtein distances with other pools
    for other_pool_key, other_pool_models in pools_of_identicals.items():
        if other_pool_key <= pool_key:
            continue  # Skip duplicate and self-comparisons
            
        unique_model_b = other_pool_models[0]
        code_b = modelDB_processed_mod_files[unique_model_b]['stripped_comments_whitespaces_source_code']
        
        try:
            # Using fuzz.ratio to get the Levenshtein Distance ratio
            distance_ratio = fuzz.ratio(code_a, code_b)
            if distance_ratio >= 75:

                # Store the distance ratio in the dictionary
                levenshtein_distances[(pool_key, other_pool_key)] = distance_ratio
            
            #print(f"Levenshtein distance between {unique_model_a} and {unique_model_b}: {distance_ratio}")  
        except Exception as e:
            print(f"Could not calculate distance between {unique_model_a} and {unique_model_b} due to {e}")

        # Update progress bar
        progress_bar.update(1)

# Close progress bar
progress_bar.close()

# Save the levenshtein_distances dictionary
with open("levenshtein_ratios.pkl", "wb") as f:
    pickle.dump(levenshtein_distances, f)

print("Finished calculating Levenshtein ratios.")