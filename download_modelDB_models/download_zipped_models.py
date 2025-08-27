import requests
import os
import pandas as pd
from tqdm import tqdm
import time
import random

# Function to download a single file
def download_file(model_id, file_path, download_url):
    try:
        # Download the file with a 50-second timeout
        response = requests.get(download_url, timeout=45)
        
        # Check if the request was successful
        if response.status_code == 200:
            # Save the file
            with open(file_path, "wb") as f:
                f.write(response.content)
            return True
        else:
            return False
    except requests.exceptions.RequestException:
        return False

# Create the directory if it doesn't exist
directory_name = "modelDB_zipped_models"
if not os.path.exists(directory_name):
    os.mkdir(directory_name)

# Load the CSV file into a Pandas DataFrame
df = pd.read_csv("model_data.csv")

# Initialize tqdm bar
pbar = tqdm(total=len(df), desc="Downloading Models", dynamic_ncols=True)

# Download the zipped files
for index, row in df.iterrows():
    model_id = row['model_id']
    file_path = f"{directory_name}/{model_id}.zip"
    
    # Skip if the model has already been downloaded
    if os.path.exists(file_path):
        pbar.write(f"{model_id}.zip already downloaded, skipping.")
        pbar.update(1)
        continue

    download_url = f"https://modeldb.science/download/{model_id}"

    # Try to download the file up to 3 times
    for attempt in range(2):
        if download_file(model_id, file_path, download_url):
            pbar.write(f"Downloaded {model_id}.zip")
            break
        else:
            pbar.write(f"Attempt {attempt + 1} failed for {model_id}.zip. Retrying...")
            # Add random delay before retrying
            #time.sleep(random.uniform(2, 10))

    # Update tqdm bar
    pbar.update(1)

    # Add random delay
    #time.sleep(random.uniform(2, 10))

# Close tqdm bar
pbar.close()

print("Download complete!")
