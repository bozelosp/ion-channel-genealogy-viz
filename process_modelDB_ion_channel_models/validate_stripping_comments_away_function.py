import pickle
import re
import sys

def strip_nmodl_comments(source_code):

    source_code = re.sub(r'^TITLE.*$', '', source_code, flags=re.MULTILINE)

    # Step 1: Identify VERBATIM blocks
    verbatim_blocks = re.findall(r'VERBATIM(.*?)ENDVERBATIM', source_code, flags=re.DOTALL)
    
    # Step 2: Replace VERBATIM blocks with placeholders
    placeholders = []
    for i, block in enumerate(verbatim_blocks):
        placeholder = f"__PLACEHOLDER__{i}__"
        placeholders.append(placeholder)
        source_code = re.sub(f"VERBATIM{re.escape(block)}ENDVERBATIM", placeholder, source_code, flags=re.DOTALL)
        
    # Step 3: Strip NMODL-specific comments
    # Remove single-line comments starting with :
    source_code = re.sub(r':.*$', '', source_code, flags=re.MULTILINE)
    # Remove multi-line comments (COMMENT ... ENDCOMMENT)
    source_code = re.sub(r'COMMENT.*?ENDCOMMENT', '', source_code, flags=re.DOTALL)
    
    # Step 4: Restore VERBATIM blocks and remove C-style comments within them
    for i, (block, placeholder) in enumerate(zip(verbatim_blocks, placeholders)):
        # Remove C-style single-line comments (// ...)
        block = re.sub(r'//.*$', '', block, flags=re.MULTILINE)
        # Remove C-style multi-line comments (/* ... */)
        block = re.sub(r'/\*.*?\*/', '', block, flags=re.DOTALL)
        # Restore the VERBATIM block
        source_code = source_code.replace(placeholder, f"VERBATIM{block}ENDVERBATIM")
    
    source_code = source_code.strip()

    return source_code

# Load the dictionary containing the source code of .mod files
with open('mod_source_code_dict.pkl', 'rb') as f:
    mod_source_code_dict = pickle.load(f)

# Load the list of .mod files with all types of comments
with open('files_with_all_comments.pkl', 'rb') as f:
    files_with_all_comments = pickle.load(f)

# Process each .mod file
for mod_file_name in files_with_all_comments:
    print(f"Processing {mod_file_name}...")
    
    # Get the original source code
    mod_data = mod_source_code_dict.get(mod_file_name, {})
    original_source_code = mod_data.get('source_code', '')
    
    # Print the original source code
    print("=== Original Source Code ===")
    print(original_source_code)
    
    # Strip comments and print the modified source code
    modified_source_code = strip_nmodl_comments(original_source_code)
    
    print("\n=== Modified Source Code (Comments Stripped) ===")
    print(modified_source_code)
    print("\n" + "="*80 + "\n")
