import pandas as pd
import sys
import glob
import os

def parse_hmmer_output(file_path):
    """
    Parse HMMER output files and convert to a pandas DataFrame.
    
    Args:
        file_path (str): Path to the HMMER output file
        
    Returns:
        pandas.DataFrame: Parsed data in tabular format
    """
    data = []
    
    with open(file_path, 'r') as f:
        # Get the KEGG ID from the filename
        kegg_id = os.path.basename(file_path).split('.')[0]
        
        for line in f:
            # Skip comments that don't contain column headers
            if line.startswith('#') and 'target name' not in line:
                continue
                
            # Skip empty lines
            if not line.strip():
                continue
                
            # Split the line into columns
            columns = line.strip().split()
            
            # Skip the header line
            if columns[0] == '#':
                continue
                
            # Extract the relevant fields
            row = {
                'kegg_id': kegg_id,
                'target_name': columns[0],
                'target_accession': columns[1],
                'query_name': columns[2],
                'query_accession': columns[3],
                'e_value': float(columns[4]),
                'score': float(columns[5]),
                'bias': float(columns[6]),
                'best_domain_e_value': float(columns[7]),
                'best_domain_score': float(columns[8]),
                'best_domain_bias': float(columns[9])
            }
            
            data.append(row)
    
    # Convert to DataFrame
    df = pd.DataFrame(data)
    return df

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_pattern> <output_file>")
        print("Example: python script.py '*.hmm_plasmid.txt' combined_results.csv")
        sys.exit(1)
        
    input_pattern = sys.argv[1]
    output_file = sys.argv[2]
    
    # Get list of all input files
    input_files = glob.glob(input_pattern)
    
    if not input_files:
        print(f"No files found matching pattern: {input_pattern}")
        sys.exit(1)
    
    # Parse and combine all files
    dfs = []
    for file in input_files:
        try:
            df = parse_hmmer_output(file)
            dfs.append(df)
            print(f"Processed: {file}")
        except Exception as e:
            print(f"Error processing {file}: {str(e)}")
    
    # Combine all DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Save based on file extension
    if output_file.endswith('.csv'):
        combined_df.to_csv(output_file, index=False)
    elif output_file.endswith('.xlsx'):
        combined_df.to_excel(output_file, index=False)
    elif output_file.endswith('.tsv'):
        combined_df.to_csv(output_file, sep='\t', index=False)
    else:
        combined_df.to_csv(output_file, sep='\t', index=False)  # Default to TSV
    
    print(f"\nProcessed {len(input_files)} files")
    print(f"Total rows in combined table: {len(combined_df)}")
    print(f"Output saved to: {output_file}")

if __name__ == "__main__":
    main()
