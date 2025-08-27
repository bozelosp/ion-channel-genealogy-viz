from bs4 import BeautifulSoup
import pandas as pd
import re

# Read the saved HTML file
with open("modeldb.science_ListByModelName.html", "r") as file:
    html_content = file.read()

# Initialize BeautifulSoup
soup = BeautifulSoup(html_content, 'html.parser')

# Initialize list to store model data
model_data = []

# Find all the model list items
model_list_items = soup.find_all("li")

# Extract href, model ID, name, and year for each model
for item in model_list_items:
    anchor = item.find("a")
    href = anchor['href']

    model_id_match = re.search(r"model=(\d+)", href)
    if model_id_match is None:
        print(f"Skipping item with href {href} as it does not include a model ID.")
        continue
    model_id = model_id_match.group(1)
    
    full_name = anchor.get_text()
    year_match = re.search(r"\((.*?)(\d{4})(.*?)\)", full_name)
    year = year_match.group(2) if year_match else "None"

    model_data.append({
        "href": href,
        "model_id": model_id,
        "full_name": full_name,
        "year": year
    })

# Create a Pandas DataFrame
df = pd.DataFrame(model_data)

# Save the DataFrame in JSON and CSV formats
df.to_json("model_data.json", orient="records", lines=True)
df.to_csv("model_data.csv", index=False)