import os
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm  # For progress bars (install with pip install tqdm)

# Configuration
INPUT_DIR = "Pathogenic_Missmatch/code55_Primer_Results"
OUTPUT_FILE = "Species_Without_Primer_Pairs.xlsx"
DEBUG_FILE = "Primer_Coverage_Debug_Info.xlsx"
MAX_WORKERS = 4  # Adjust based on your CPU cores

def load_primer_file(file_path):
    """Load and validate a single primer file"""
    try:
        df = pd.read_excel(file_path)
        # Validate required columns
        required_cols = {'Species_Name', 'Sequence_Length', 'Primer_Pair_Present'}
        if not required_cols.issubset(df.columns):
            missing = required_cols - set(df.columns)
            print(f"âš ï¸ Missing columns in {os.path.basename(file_path)}: {missing}")
            return None
        return df
    except Exception as e:
        print(f"âš ï¸ Error loading {file_path}: {str(e)}")
        return None

def analyze_primer_coverage():
    """Optimized primer coverage analysis"""
    # Get all primer directories
    try:
        primer_dirs = [d for d in os.listdir(INPUT_DIR) 
                     if os.path.isdir(os.path.join(INPUT_DIR, d))]
    except FileNotFoundError:
        print(f"Error: Input directory '{INPUT_DIR}' not found")
        return pd.DataFrame(), pd.DataFrame()

    if not primer_dirs:
        print("âš ï¸ No primer directories found in", INPUT_DIR)
        return pd.DataFrame(), pd.DataFrame()

    # Prepare file paths for parallel processing
    file_paths = [os.path.join(INPUT_DIR, p, f"{p}_results.xlsx") for p in primer_dirs]
    
    # Load all files in parallel with progress bar
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        results = list(tqdm(executor.map(load_primer_file, file_paths), 
                      total=len(file_paths), desc="Loading files"))

    # Filter out None results (failed loads)
    valid_dfs = [(primer_dirs[i], df) for i, df in enumerate(results) if df is not None]
    
    if not valid_dfs:
        print("âš ï¸ No valid primer files found")
        return pd.DataFrame(), pd.DataFrame()

    # Get all unique species from the first valid file
    all_species = valid_dfs[0][1]['Species_Name'].unique()
    
    # Initialize coverage tracking with numpy arrays for speed
    coverage_counts = {species: 0 for species in all_species}
    total_primers = len(valid_dfs)
    
    # Process each primer file
    for primer_name, df in tqdm(valid_dfs, desc="Analyzing primers"):
        # Vectorized operation to count present primers
        present_species = df.loc[df['Primer_Pair_Present'] == 'Yes', 'Species_Name']
        for species in present_species.unique():
            if species in coverage_counts:
                coverage_counts[species] += 1
    
    # Create debug DataFrame efficiently
    debug_data = []
    for species in sorted(coverage_counts.keys()):
        present = coverage_counts[species]
        debug_data.append({
            'Species': species,
            'Primers Present': present,
            'Total Primers': total_primers,
            'Coverage Percentage': (present / total_primers) * 100 if total_primers > 0 else 0,
            'Completely Missing': present == 0
        })
    
    debug_df = pd.DataFrame(debug_data)
    
    # Get sequence lengths from first valid file
    first_df = valid_dfs[0][1]
    seq_lengths = first_df.set_index('Species_Name')['Sequence_Length'].to_dict()
    
    # Filter for missing species and add sequence lengths
    missing_species = debug_df[debug_df['Completely Missing']].copy()
    if not missing_species.empty:
        missing_species['Sequence Length'] = missing_species['Species'].map(seq_lengths)
        missing_species = missing_species[['Species', 'Sequence Length',
                                         'Primers Present', 'Total Primers',
                                         'Coverage Percentage']]
    else:
        missing_species = pd.DataFrame(columns=['Species', 'Sequence Length',
                                              'Primers Present', 'Total Primers',
                                              'Coverage Percentage'])
    
    return missing_species, debug_df

def main():
    print("ðŸš€ Starting primer coverage analysis...")
    
    # Analyze primer coverage
    result_df, debug_df = analyze_primer_coverage()
    
    # Save results
    try:
        if not result_df.empty:
            result_df.to_excel(OUTPUT_FILE, index=False)
            print(f"\nâœ… Found {len(result_df)} species without any primer pairs")
            print("Saved results to:", OUTPUT_FILE)
        else:
            print("\nðŸŽ‰ All species have at least one primer pair present!")
            # Create empty file for consistency
            pd.DataFrame().to_excel(OUTPUT_FILE, index=False)
        
        # Save debug info
        debug_df.to_excel(DEBUG_FILE, index=False)
        print("ðŸ“Š Detailed coverage information saved to:", DEBUG_FILE)
        
        # Show coverage statistics
        if not debug_df.empty:
            print("\nðŸ“ˆ Coverage Statistics:")
            print(f"â€¢ Mean coverage: {debug_df['Coverage Percentage'].mean():.1f}%")
            print(f"â€¢ Minimum coverage: {debug_df['Coverage Percentage'].min():.1f}%")
            print(f"â€¢ Species with <50% coverage: {len(debug_df[debug_df['Coverage Percentage'] < 50])}")
            print(f"â€¢ Species with <25% coverage: {len(debug_df[debug_df['Coverage Percentage'] < 25])}")
    except Exception as e:
        print(f"\nâš ï¸ Error saving results: {str(e)}")

if __name__ == "__main__":
    main()
