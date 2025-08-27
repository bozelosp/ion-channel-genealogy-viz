from bs4 import BeautifulSoup
import os
import re

# Initialize an empty dictionary to hold the mapping
class_to_model_name_to_ICG_detailed_page = {}

# List of HTML files in the 'icg_links_to_detailed_pages' directory
html_files = [
    "ICGenealogy - Class - K.html",
    "ICGenealogy - Class - Na.html",
    "ICGenealogy - Class - Ca.html",
    "ICGenealogy - Class - IH.html",
    "ICGenealogy - Class - KCa.html",
    "ICGenealogy - Class - Other.html"  # Make sure the name is correct
]

# Directory path
dir_path = "icg_links_to_detailed_pages"

for html_file in html_files:
    # Extract ion channel class from the file name
    ion_channel_class = re.search(r'Class - (.+).html', html_file).group(1)

    # Initialize the inner dictionary for this ion_channel_class
    class_to_model_name_to_ICG_detailed_page[ion_channel_class] = {}

    with open(os.path.join(dir_path, html_file), 'r', encoding='utf-8') as f:
        soup = BeautifulSoup(f, 'lxml')
        
        # Find all list items
        list_items = soup.find_all('li')
        
        # Loop through the list items and populate the dictionary
        for item in list_items:
            anchor = item.find('a')
            if anchor:
                href = anchor['href']
                model_name = anchor.text
                
                # valid model names begin with \d+-
                if not re.match(r'\d+-', model_name):
                    continue
                
                class_to_model_name_to_ICG_detailed_page[ion_channel_class][model_name] = 'https://icg.neurotheory.ox.ac.uk' + href

import pickle

# Save the dictionary to a pickle file
with open('class_to_model_name_to_ICG_detailed_page.pkl', 'wb') as f:
    pickle.dump(class_to_model_name_to_ICG_detailed_page, f)

import json

# Save the dictionary to a json file
with open('class_to_model_name_to_ICG_detailed_page.json', 'w') as f:
    json.dump(class_to_model_name_to_ICG_detailed_page, f, indent=4)