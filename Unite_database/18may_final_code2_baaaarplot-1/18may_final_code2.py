import pandas as pd
import matplotlib.pyplot as plt
import os

# Constants
TOTAL_ALL_SPECIES = 168030
TOTAL_PATHOGENIC_SPECIES = 1640

# Output folder to save plots and explanation
output_folder = "fff18may_code2_baaaarplot-1"
os.makedirs(output_folder, exist_ok=True)

def calculate_percentages(file_path, total_species_count):
    """
    Calculate percentages and counts from primer search results.
    
    Args:
        file_path (str): Path to the Excel file with primer results
        total_species_count (int): Total number of species in the database
    
    Returns:
        pd.DataFrame: DataFrame with primer pairs, counts and percentages
    """
    df = pd.read_excel(file_path)
    primer_pairs = [col for col in df.columns if " + " in col and "Fwd_Pos" not in col and "Rev_Pos" not in col]
    
    results = []
    for pair in primer_pairs:
        count = df[pair].fillna("").str.contains("Found").sum()
        percentage = (count / total_species_count) * 100
        results.append({
            "Primer Pair": pair,
            "Count": count,
            "Percentage": round(percentage, 2)
        })
    
    return pd.DataFrame(results)

def plot_bar(df_percent, title, output_name):
    """
    Create bar plot with counts and percentages displayed
    
    Args:
        df_percent (pd.DataFrame): DataFrame with results
        title (str): Plot title
        output_name (str): Output filename
    """
    plt.figure(figsize=(18, 7))  # Increased width from 14 to 18
    bars = plt.bar(df_percent["Primer Pair"], df_percent["Percentage"], 
                  color="red", width=0.6)  # Increased bar width from 0.3 to 0.6
    
    # Add labels with both count and percentage
    for bar in bars:
        height = bar.get_height()
        idx = bars.index(bar)
        count = df_percent.loc[idx, "Count"]
        plt.text(bar.get_x() + bar.get_width()/2, height + 0.5,
                 f"{count}\n({height:.2f}%)", 
                 ha="center", va="bottom", fontsize=8)
    
    plt.xlabel("Primer Pair Name", fontsize=12)
    plt.ylabel("Percentage of Species (%)", fontsize=12)
    plt.title(title, fontsize=14, pad=20)
    plt.ylim(0, 100)  # Set y-axis limit to 20
    plt.yticks([0, 5, 10, 15, 20])  # Set specific y-axis ticks
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, output_name), dpi=300)
    plt.close()

# File paths
all_file = "/home/abhishek/ITS/UNITE-REFINED_PRIMERS/primer_search_results_code1/primer_pair_results_all.xlsx"
pathogens_file = "/home/abhishek/ITS/UNITE-REFINED_PRIMERS/primer_search_results_code1/primer_pair_results_pathogens.xlsx"

# Calculate percentages and counts
df_all = calculate_percentages(all_file, TOTAL_ALL_SPECIES)
df_pathogens = calculate_percentages(pathogens_file, TOTAL_PATHOGENIC_SPECIES)

# Save results to Excel
df_all.to_excel(os.path.join(output_folder, "all_species_results.xlsx"), index=False)
df_pathogens.to_excel(os.path.join(output_folder, "pathogens_results.xlsx"), index=False)

# Plot graphs
plot_bar(df_all, "Primer Pair Coverage in All Fungal Species", "all_species_coverage.png")
plot_bar(df_pathogens, "Primer Pair Coverage in Pathogenic Fungi", "pathogens_coverage.png")

# Save detailed explanation
explanation_file = os.path.join(output_folder, "calculation_explanation.txt")
with open(explanation_file, "w") as f:
    f.write("Primer Pair Analysis - Detailed Explanation\n")
    f.write("="*50 + "\n\n")
    
    f.write("1. COUNT CALCULATION:\n")
    f.write("   - For each primer pair column, count how many cells contain 'Found'\n")
    f.write("   - Missing/empty cells are not counted\n\n")
    
    f.write("2. PERCENTAGE CALCULATION:\n")
    f.write("   Formula: (Count / Total Species) × 100\n\n")
    
    f.write("3. SPECIFIC VALUES USED:\n")
    f.write(f"   - Total fungal species: {TOTAL_ALL_SPECIES}\n")
    f.write(f"   - Total pathogenic fungi: {TOTAL_PATHOGENIC_SPECIES}\n\n")
    
    f.write("4. OUTPUT FILES:\n")
    f.write("   - Excel files with raw counts and percentages\n")
    f.write("   - Bar plots showing visual comparison\n")
    f.write("   - This explanation file\n")

# Print summary to console
print("\n" + "="*50)
print("ANALYSIS SUMMARY")
print("="*50)
print(f"\nAll Species Results (Total={TOTAL_ALL_SPECIES}):")
print(df_all[["Primer Pair", "Count", "Percentage"]].to_string(index=False))
print(f"\nPathogenic Fungi Results (Total={TOTAL_PATHOGENIC_SPECIES}):")
print(df_pathogens[["Primer Pair", "Count", "Percentage"]].to_string(index=False))
print(f"\nAll results saved to: {os.path.abspath(output_folder)}")
