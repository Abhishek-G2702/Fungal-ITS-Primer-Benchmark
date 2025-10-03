import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ====================== CONFIGURATION ======================
INPUT_DIR = "code1_Primer_Results"  # Directory from the first analysis
OUTPUT_FILE = "code2_Primer_Pair_Presence_Summary.xlsx"
GRAPH_FILE = "Primer_Pair_Presence.png"

# ====================== CORE FUNCTIONS ======================
def calculate_primer_presence():
    """Calculate primer pair presence statistics across all primer pairs"""
    primer_stats = []
    
    # Get all primer directories
    primer_dirs = [d for d in os.listdir(INPUT_DIR) 
                  if os.path.isdir(os.path.join(INPUT_DIR, d))]
    
    for primer_name in primer_dirs:
        # Find the Excel file (matches naming convention from first code)
        excel_file = os.path.join(INPUT_DIR, primer_name, f"{primer_name}_results.xlsx")
        
        if not os.path.exists(excel_file):
            continue
            
        df = pd.read_excel(excel_file)
        
        # Calculate statistics
        total_species = len(df)
        present_species = df[df['Primer_Pair_Present'] == 'Yes'].shape[0]
        percentage = (present_species / total_species) * 100 if total_species > 0 else 0
        
        primer_stats.append({
            'Primer Pair': primer_name,
            'Species with Primer Pair': present_species,
            'Total Species': total_species,
            'Percentage': round(percentage, 2)
        })
    
    # Sort by percentage (descending)
    primer_stats.sort(key=lambda x: x['Percentage'], reverse=True)
    
    return pd.DataFrame(primer_stats)

def create_bar_chart(df):
    """Create a bar chart showing primer pair presence percentages"""
    plt.figure(figsize=(12, 6))
    
    # Create bars with color gradient
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(df)))
    bars = plt.bar(df['Primer Pair'], df['Percentage'], color=colors)
    
    # Add percentage labels on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{height:.1f}%',
                 ha='center', va='bottom')
    
    # Customize the plot
    plt.title('Percentage of Species with Each Primer Pair (<2 Mismatches)', fontsize=14)
    plt.xlabel('Primer Pair', fontsize=12)
    plt.ylabel('Percentage (%)', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.ylim(0, 100)  # Percentage scale
    
    # Add grid lines
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save and close
    plt.savefig(GRAPH_FILE, dpi=300)
    plt.close()

# ====================== MAIN EXECUTION ======================
def main():
    # Calculate primer presence statistics
    primer_df = calculate_primer_presence()
    
    # Save to Excel
    primer_df.to_excel(OUTPUT_FILE, index=False)
    print(f"Primer pair statistics saved to {OUTPUT_FILE}")
    
    # Create visualization
    create_bar_chart(primer_df)
    print(f"Visualization saved to {GRAPH_FILE}")
    
    # Print summary
    print("\nPrimer Pair Performance Summary:")
    print(f"- Best performing: {primer_df.iloc[0]['Primer Pair']} ({primer_df.iloc[0]['Percentage']}%)")
    print(f"- Worst performing: {primer_df.iloc[-1]['Primer Pair']} ({primer_df.iloc[-1]['Percentage']}%)")
    print(f"- Average presence: {primer_df['Percentage'].mean():.1f}% across all primer pairs")

if __name__ == "__main__":
    main()
