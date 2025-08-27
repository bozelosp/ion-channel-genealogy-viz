import pickle
import re

def check_all_comment_types_exist(source_code):
    # Check for single-line comments starting with "//"
    has_double_slash_comments = re.search(r'//', source_code)
    
    # Check for single-line comments starting with ":"
    has_colon_comments = re.search(r':', source_code)
    
    # Check for multi-line comments starting with "/*" and ending with "*/"
    has_star_comments = re.search(r'/\*.*?\*/', source_code, re.DOTALL)
    
    # Check for multi-line comments starting with "COMMENT" and ending with "ENDCOMMENT"
    has_comment_block = re.search(r'COMMENT.*?ENDCOMMENT', source_code, re.DOTALL)
    
    return all([has_double_slash_comments, has_colon_comments, has_star_comments, has_comment_block])

# Load the dictionary containing the source code of .mod files
with open('mod_source_code_dict.pkl', 'rb') as f:
    mod_source_code_dict = pickle.load(f)

# Initialize a list to store the names of .mod files containing all types of comments
files_with_all_comments = []

# Check each .mod file
for mod_file_name, mod_data in mod_source_code_dict.items():
    source_code = mod_data.get('source_code', '')
    if check_all_comment_types_exist(source_code):
        files_with_all_comments.append(mod_file_name)
        
        # Stop when we find 10 such files
        if len(files_with_all_comments) >= 10:
            break

# Save the names of the .mod files containing all types of comments
with open('files_with_all_comments.pkl', 'wb') as f:
    pickle.dump(files_with_all_comments, f)

print(f"Found {len(files_with_all_comments)} files with all types of comments. Saved as 'files_with_all_comments.pkl'.")