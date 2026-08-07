import pandas as pd
import glob
import os
from pathlib import Path
import re

def extract_significant_rows():
    """
    Extract rows where qval, FDR_theoretical, or bio_FDR_theoretical < 0.05
    from all matching files and merge into a single CSV file.
    """
    
    # Define file pattern and output path
    root = os.environ.get("FDRREG_RESULTS_DIR", "data/pipeline")
    base_path = os.path.join(root, "adhd2019", "07.metaxcan_fdrreg", "01.fdrreg_results")
    file_pattern = os.path.join(base_path, "*/*.gene.bio.fdrreg.txt")
    
    # Output directory and file
    output_dir = os.path.join(root, "01.extra.analysis", "11.summary_table")
    output_file = os.path.join(output_dir, "significant_fdrreg_results.csv")
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Find all matching files
    all_files = glob.glob(file_pattern)
    
    if not all_files:
        print(f"No files found matching pattern: {file_pattern}")
        return None
    
    print(f"Found {len(all_files)} files to process")
    
    # List to store dataframes from each file
    dataframes = []
    
    # Process each file
    for file_path in all_files:
        try:
            # Extract tissue/region name from file path
            tissue_name = os.path.basename(os.path.dirname(file_path))
            
            # Read the file with flexible whitespace separator
            # Using regex pattern for one or more whitespace characters
            df = pd.read_csv(file_path, sep=r'\s+', engine='python')
            
            # Clean column names - remove any leading/trailing whitespace
            df.columns = df.columns.str.strip()
            
            # Check if required columns exist
            required_columns = ['qval', 'FDR_theoretical', 'bio_FDR_theoretical']
            missing_columns = [col for col in required_columns if col not in df.columns]
            
            if missing_columns:
                print(f"Warning: Missing columns in {tissue_name}: {missing_columns}")
                print(f"Available columns: {df.columns.tolist()}")
                continue
            
            # Convert columns to numeric, handling errors
            for col in required_columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
            
            # Apply filters: qval < 0.05 OR FDR_theoretical < 0.05 OR bio_FDR_theoretical < 0.05
            mask = (df['qval'] < 0.05) | (df['FDR_theoretical'] < 0.05) | (df['bio_FDR_theoretical'] < 0.05)
            filtered_df = df[mask].copy()
            
            # Add tissue information
            filtered_df['tissue'] = tissue_name
            
            # Add source file information
            filtered_df['source_file'] = os.path.basename(file_path)
            
            if not filtered_df.empty:
                dataframes.append(filtered_df)
                print(f"Processed {tissue_name}: {len(filtered_df)} significant rows")
            
        except Exception as e:
            print(f"Error processing file {file_path}:")
            print(f"  Error: {str(e)}")
            print(f"  File exists: {os.path.exists(file_path)}")
            if os.path.exists(file_path):
                # Try to read the file and show its structure
                try:
                    with open(file_path, 'r') as f:
                        first_line = f.readline().strip()
                    print(f"  First line: {first_line}")
                except:
                    pass
            continue
    
    # Merge all dataframes
    if not dataframes:
        print("No significant rows found in any file")
        return None
    
    merged_df = pd.concat(dataframes, ignore_index=True)
    
    # Reorder columns to have tissue information first
    cols = ['tissue', 'ENSEMBL_GENE_ID', 'gene_name', 'zscore', 'pvalue', 
            'qval', 'FDR_theoretical', 'FDR_empirical', 'bio_FDR_theoretical', 
            'bio_FDR_empirical', 'lasso_FDR_theoretical', 'lasso_FDR_empirical', 
            'source_file']
    
    # Only include columns that exist in the dataframe
    existing_cols = [col for col in cols if col in merged_df.columns]
    merged_df = merged_df[existing_cols]
    
    # Save to CSV
    merged_df.to_csv(output_file, index=False)
    
    print(f"\nSummary:")
    print(f"Total files processed: {len(all_files)}")
    print(f"Total significant rows: {len(merged_df)}")
    print(f"Unique genes: {merged_df['ENSEMBL_GENE_ID'].nunique()}")
    print(f"Unique tissues: {merged_df['tissue'].nunique()}")
    print(f"Results saved to: {output_file}")
    
    return merged_df

if __name__ == "__main__":
    # Execute the function
    result_df = extract_significant_rows()
    
    # Display sample of results
    if result_df is not None:
        print("\nSample of extracted data:")
        print(result_df.head(10))
        
        # Save summary statistics
        summary_file = os.path.join(output_dir, "summary_statistics.txt")
        with open(summary_file, 'w') as f:
            f.write("FDRREG Results Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Total significant rows: {len(result_df)}\n")
            f.write(f"Unique genes: {result_df['ENSEMBL_GENE_ID'].nunique()}\n")
            f.write(f"Unique tissues: {result_df['tissue'].nunique()}\n\n")
            
            f.write("Distribution by tissue:\n")
            tissue_counts = result_df['tissue'].value_counts()
            for tissue, count in tissue_counts.items():
                f.write(f"{tissue}: {count} rows\n")
            
            f.write("\nFilter criteria applied:\n")
            f.write("qval < 0.05 OR FDR_theoretical < 0.05 OR bio_FDR_theoretical < 0.05\n")
        
        print(f"Summary statistics saved to: {summary_file}")
