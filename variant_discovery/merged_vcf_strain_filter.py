import pandas as pd
import sys

def filter_variants(input_file, target_strains_list):
    try:
        # Load the data - assuming tab-separated
        df = pd.read_csv(input_file, sep='\t')
    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
        sys.exit(1)

    # Parse target strains
    targets = [s.strip() for s in target_strains_list.split(',')]
    
    # Identify strain columns, has to start with BL
    all_strain_cols = [col for col in df.columns if col.startswith('BL')]
    other_strains = [s for s in all_strain_cols if s not in targets]

    # Check for missing targets
    missing = [s for s in targets if s not in df.columns]
    if missing:
        print(f"Warning: These strains were not found in the table: {missing}")
        targets = [s for s in targets if s in df.columns]

    if not targets:
        print("Error: No valid target strains found. Exiting.")
        sys.exit(1)

    # Filtering Logic:
    # 1. Target strains must NOT be '0/0'
    target_condition = (df[targets] != "0/0").all(axis=1)
    # 2. All other 'BL' strains MUST be '0/0'
    other_condition = (df[other_strains] == "0/0").all(axis=1)

    # Combine and Filter
    filtered_df = df[target_condition & other_condition]
    
    # Output results
    output_name = f"filtered_{input_file}"
    filtered_df.to_csv(output_name, sep='\t', index=False)
    print(f"Success! Found {len(filtered_df)} matching rows. Saved to: {output_name}")

if __name__ == "__main__":
    # Check if the right number of arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python filter_strains.py <filename> <strain1,strain2,...>")
        print("Example: python filter_strains.py table.txt BL156,BL200,BL300")
        sys.exit(1)

    # Assign arguments
    file_arg = sys.argv[1]
    strains_arg = sys.argv[2]
    
    filter_variants(file_arg, strains_arg)