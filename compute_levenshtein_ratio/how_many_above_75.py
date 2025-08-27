import pickle

# now let's read it again and print the number of entries in the dictionary
with open("levenshtein_ratios.pkl", "rb") as f:
    levenshtein_distances = pickle.load(f)

print(f"Number of entries in the levenshtein_distances dictionary: {len(levenshtein_distances)}")