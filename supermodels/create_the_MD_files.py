import pickle
import os

# Define the directory where you want to save the markdown files
output_directory = 'output_markdown_files'

# Create the output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Load the supermodel_data.pkl file
with open("supermodel_data.pkl", "rb") as f:
    supermodels = pickle.load(f)

# Function to create markdown content for each entry
def create_markdown_content(entry_name, entry_data):
    markdown_content = f"# {entry_name}\n\n"  # Start with a title
    for key, value in entry_data.items():
        markdown_content += f"## {key}\n\n"
        
        # Convert different types of values to string for addition to markdown content
        if isinstance(value, dict):
            for sub_key, sub_value in value.items():
                markdown_content += f"- **{sub_key}**: {sub_value}\n"
        elif isinstance(value, list):
            for item in value:
                markdown_content += f"- {item}\n"
        elif isinstance(value, bool) or value is None:
            markdown_content += f"- {value}\n"
        else:
            markdown_content += f"{value}\n"
        
        markdown_content += "\n"  # Add a newline for readability
    
    return markdown_content

# Loop through the nested dictionary and create a markdown file for each entry
for ion_class, models in supermodels.items():
    ion_class_directory = os.path.join(output_directory, ion_class)
    os.makedirs(ion_class_directory, exist_ok=True)
    
    for model_name, model_data in models.items():
        markdown_content = create_markdown_content(model_name, model_data)
        file_name = f"{ion_class_directory}/{model_name}.md"
        
        with open(file_name, 'w') as f:
            f.write(markdown_content)

print(f"Markdown files have been written to {output_directory}/")

