import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# Define the primer pairs
primer_pairs = {
    "ITS1F (F)-ITS2(R)": {
        "Forward": "CTTGGTCATTTAGAGGAAGTAA",
        "Reverse": "GCTGCGTTCTTCATCGATGC"
    },
    "ITS3 (F)-ITS4(R)": {
        "Forward": "GCATCGATGAAGAACGCAGC",
        "Reverse": "TCCTCCGCTTATTGATATGC"
    },
    "ITS86F(F)-ITS4(R)": {
        "Forward": "GTGAATCATCGAATCTTTGAA",
        "Reverse": "TCCTCCGCTTATTGATATGC"
    },
    "ITS5(F)-ITS4(R)": {
        "Forward": "GGAAGTAAAAGTCGTAACAAGG",
        "Reverse": "TCCTCCGCTTATTGATATGC"
    },
    "NS7(F)-ITS2(R)": {
        "Forward": "GAGGCAATAACAGGTCTGTGATGC",
        "Reverse": "GCTGCGTTCTTCATCGATGC"
    },
    "ITS3(F)-LR3(R)": {
        "Forward": "GCATCGATGAAGAACGCAGC",
        "Reverse": "CCGTGTTTCAAGACGGG"
    },
}

# Define CBHAPCA genera
cbhpaca_genera = ['Cryptococcus', 'Blastomyces', 'Histoplasma',
                  'Aspergillus', 'Candida']

# Number of top abundant primers to consider
TOP_ABUNDANT_PRIMERS = 4

def find_related_genus(species_name):
    """Find if species name contains any of our target genera names"""
    lower_species = species_name.lower()
    for genus in cbhpaca_genera:
        if genus.lower() in lower_species:
            return genus
    return None

def load_and_categorize_sequences(fasta_file):
    """Load sequences from FASTA and categorize by genus with automatic related genus detection"""
    sequences = []
    genera_counts = defaultdict(int)
    cbhpaca_counts = {genus: 0 for genus in cbhpaca_genera}
    total_species = 0
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        full_species = record.description.split('|')[0].strip()
        genus_species = full_species.split('_')[0] + '_' + full_species.split('_')[1] if '_' in full_species else full_species
        
        original_genus = full_species.split('_')[0].strip().capitalize()
        genera_counts[original_genus] += 1
        total_species += 1
        
        genus_found = None
        for genus in cbhpaca_genera:
            if original_genus.lower().startswith(genus.lower()):
                genus_found = genus
                break
        
        if not genus_found:
            genus_found = find_related_genus(genus_species)
        
        if genus_found:
            sequences.append({
                'id': record.id,
                'description': record.description,
                'genus': genus_found,
                'species': full_species.replace('_', ' '),
                'sequence': str(record.seq).upper(),
                'original_genus': original_genus
            })
            cbhpaca_counts[genus_found] += 1
    
    return sequences, dict(genera_counts), cbhpaca_counts, total_species

def find_primer_positions(sequence, primer):
    """Find all positions of a primer in a sequence (forward strand)"""
    seq_upper = sequence.upper()
    primer_upper = primer.upper()
    positions = []
    start = 0
    while True:
        pos = seq_upper.find(primer_upper, start)
        if pos == -1:
            break
        positions.append((pos + 1, pos + len(primer_upper)))
        start = pos + 1
    return positions

def find_reverse_complement_positions(sequence, primer):
    """Find positions of reverse primer by searching reverse complement"""
    rev_complement = str(Seq(primer).reverse_complement())
    return find_primer_positions(sequence, rev_complement)

def analyze_primers_for_sequences(sequences):
    """Analyze all primer pairs for all sequences"""
    results = []
    
    for seq_data in sequences:
        species_results = {
            'id': seq_data['id'],
            'description': seq_data['description'],
            'genus': seq_data['genus'],
            'species': seq_data['species'],
            'original_genus': seq_data['original_genus'],
            'primer_count': 0
        }
        
        for primer_name, primers in primer_pairs.items():
            fwd_positions = find_primer_positions(seq_data['sequence'], primers['Forward'])
            rev_positions = find_reverse_complement_positions(seq_data['sequence'], primers['Reverse'])
            both_present = len(fwd_positions) > 0 and len(rev_positions) > 0
            
            species_results[f'{primer_name}_present'] = 1 if both_present else 0
            species_results[f'{primer_name}_fwd_positions'] = '; '.join([f"{s}-{e}" for s, e in fwd_positions]) if fwd_positions else 'Not found'
            species_results[f'{primer_name}_rev_positions'] = '; '.join([f"{s}-{e}" for s, e in rev_positions]) if rev_positions else 'Not found'
            
            if both_present:
                species_results['primer_count'] += 1
        
        results.append(species_results)
    
    return results

def identify_special_species(all_results, top_abundant_primers):
    """Identify species with no primers and species missing all top abundant primers"""
    no_primer_species = []
    missing_all_top_primers_species = defaultdict(list)
    
    for species in all_results:
        if species['primer_count'] == 0:
            no_primer_species.append({
                'Species': species['species'],
                'Assigned Genus': species['genus'],
                'Original Genus': species['original_genus'],
                'ID': species['id'],
                'Description': species['description']
            })
        
        genus = species['genus']
        if genus in top_abundant_primers:
            top_primers = top_abundant_primers[genus]
            missing_all = all(species[f'{primer}_present'] == 0 for primer in top_primers)
            
            if missing_all and top_primers:
                present_primers = [p for p in primer_pairs if species[f'{p}_present'] == 1]
                missing_all_top_primers_species[genus].append({
                    'Species': species['species'],
                    'Assigned Genus': genus,
                    'Original Genus': species['original_genus'],
                    'ID': species['id'],
                    'Description': species['description'],
                    'Missing Top Primers': ', '.join(top_primers),
                    'Other Primers Present': ', '.join(present_primers) if present_primers else 'None'
                })
    
    return no_primer_species, missing_all_top_primers_species

def create_genus_summary(genus_data, genus_name, total_species):
    """Create summary statistics for a genus"""
    summary = []
    
    for primer_name in primer_pairs.keys():
        both_present = sum(1 for species in genus_data if species[f'{primer_name}_present'] == 1)
        percentage = (both_present / total_species) * 100 if total_species > 0 else 0
        
        summary.append({
            'Primer Name': primer_name,
            'Count of Species with Both Primers': both_present,
            'Total Species in Genus': total_species,
            'Percentage': percentage
        })
    
    return pd.DataFrame(summary)

def determine_top_abundant_primers(output_dir, top_n=TOP_ABUNDANT_PRIMERS):
    """Determine the top N most abundant primers for each genus"""
    top_abundant = {}
    
    for genus in cbhpaca_genera:
        summary_file = os.path.join(output_dir, genus, f"{genus}_summary.xlsx")
        if os.path.exists(summary_file):
            df = pd.read_excel(summary_file)
            if not df.empty:
                top_primers = df.nlargest(top_n, 'Percentage')['Primer Name'].tolist()
                top_abundant[genus] = top_primers
    
    return top_abundant

def save_percentage_calculation(genus_name, summary_df, output_dir):
    """Save percentage calculation details to a text file"""
    filename = os.path.join(output_dir, f"{genus_name}_percentage_calculations.txt")
    with open(filename, 'w') as f:
        f.write(f"Percentage Calculations for {genus_name}\n")
        f.write("="*50 + "\n\n")
        
        for _, row in summary_df.iterrows():
            f.write(f"Primer Pair: {row['Primer Name']}\n")
            f.write(f"Formula: (Number of species with both primers present / Total number of species in genus) * 100\n")
            f.write(f"Calculation: ({row['Count of Species with Both Primers']} / {row['Total Species in Genus']}) * 100 = {row['Percentage']:.2f}%\n\n")

def create_genus_barplot(genus_data, genus_name, output_dir):
    """Create a bar plot showing primer pair abundance in a genus"""
    primer_counts = {primer: 0 for primer in primer_pairs}
    total_species = len(genus_data)
    
    for primer in primer_pairs:
        primer_counts[primer] = sum(1 for species in genus_data if species[f'{primer}_present'] == 1)
    
    if not primer_counts:
        return
    
    df = pd.DataFrame.from_dict(primer_counts, orient='index', columns=['Count'])
    df['Percentage'] = (df['Count'] / total_species) * 100
    
    plt.figure(figsize=(12, 6))
    bars = plt.bar(df.index, df['Percentage'], color='#D8BFD8')
    
    plt.title(f"Primer Pair Presence in {genus_name} Species (n={total_species})")
    plt.xlabel("Primer Pairs")
    plt.ylabel("Percentage of Species with Primer Pair Present")
    plt.xticks(rotation=45, ha='right')
    plt.ylim(0, 100)
    
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{height:.1f}%',
                 ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{genus_name}_primer_abundance.png"), dpi=300)
    plt.close()

def create_genus_primer_heatmap(all_results, cbhpaca_counts, output_dir):
    """Create a heatmap with genera as rows and primer pairs as columns"""
    heatmap_data = []
    
    for genus in cbhpaca_genera:
        genus_species = [s for s in all_results if s['genus'] == genus]
        total_species = cbhpaca_counts[genus]
        
        if total_species == 0:
            continue
            
        for primer_name in primer_pairs:
            count = sum(1 for s in genus_species if s[f'{primer_name}_present'] == 1)
            percentage = (count / total_species) * 100
            
            heatmap_data.append({
                'Genus': genus,
                'Primer Pair': primer_name,
                'Percentage': percentage,
                'Count': count,
                'Total Species': total_species
            })
    
    if not heatmap_data:
        print("No data available for heatmap")
        return
    
    heatmap_df = pd.DataFrame(heatmap_data)
    pivot_table = heatmap_df.pivot(index='Genus', columns='Primer Pair', values='Percentage')
    pivot_table = pivot_table.sort_index()
    pivot_table = pivot_table[sorted(pivot_table.columns)]
    
    plt.figure(figsize=(12, 8))
    
    # Create smooth white to red gradient
    colors = ["#FFFFFF", "#FFCCCC", "#FF9999", "#FF6666", "#FF3333", "#FF0000"]
    cmap = LinearSegmentedColormap.from_list("white_to_red", colors)
    
    ax = sns.heatmap(
        pivot_table, 
        annot=True, 
        fmt=".1f",
        cmap=cmap,
        linewidths=0.5,
        linecolor='lightgray',
        cbar_kws={
            'label': 'Percentage of Species',
            'ticks': [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        },
        vmin=0,
        vmax=100,
        annot_kws={
            'size': 10,
            'color': 'black'
        },
        square=True
    )
    
    # Force 0.0% values to be pure white
    for i in range(len(pivot_table)):
        for j in range(len(pivot_table.columns)):
            if pivot_table.iloc[i, j] == 0:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, color='white', lw=0.5, ec='lightgray'))
                ax.text(j + 0.5, i + 0.5, "0.0", 
                       ha="center", va="center", 
                       color="black", fontsize=10)
    
    cbar = ax.collections[0].colorbar
    cbar.set_label('Percentage of Species', rotation=270, labelpad=25, fontsize=12)
    
    plt.title('Primer Pair Presence by Pathogenic Genus', pad=20, fontsize=16, weight='bold')
    plt.xlabel('Primer Pairs', fontsize=14, labelpad=15)
    plt.ylabel('Pathogenic Genera', fontsize=14, labelpad=15)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(rotation=0, fontsize=12)
    
    plt.tight_layout()
    
    heatmap_path = os.path.join(output_dir, 'genus_primer_heatmap.png')
    plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    heatmap_df.to_excel(os.path.join(output_dir, 'genus_primer_heatmap_data.xlsx'), index=False)
    print(f"Heatmap saved to {heatmap_path}")

def main(fasta_file, output_dir):
    """Main analysis function"""
    os.makedirs(output_dir, exist_ok=True)
    
    sequences, all_genera_counts, cbhpaca_counts, total_species = load_and_categorize_sequences(fasta_file)
    
    pd.DataFrame.from_dict(all_genera_counts, orient='index', columns=['Count']).to_excel(
        os.path.join(output_dir, 'all_genera_counts.xlsx'))
    
    pd.DataFrame.from_dict(cbhpaca_counts, orient='index', columns=['Count']).to_excel(
        os.path.join(output_dir, 'CBHPACA_genus_counts.xlsx'))
    
    mapping_data = []
    for seq in sequences:
        mapping_data.append({
            'Original Genus': seq['original_genus'],
            'Assigned Genus': seq['genus'],
            'Species': seq['species'],
            'ID': seq['id']
        })
    pd.DataFrame(mapping_data).to_excel(
        os.path.join(output_dir, 'genus_mapping.xlsx'), index=False)
    
    with open(os.path.join(output_dir, 'total_species_count.txt'), 'w') as f:
        f.write(f"Total number of pathogenic fungal species analyzed: {total_species}")
    
    all_results = analyze_primers_for_sequences(sequences)
    results_df = pd.DataFrame(all_results)
    
    for genus in cbhpaca_genera:
        genus_data = [s for s in all_results if s['genus'] == genus]
        if not genus_data:
            print(f"No sequences found for genus {genus}")
            continue
        
        genus_dir = os.path.join(output_dir, genus)
        os.makedirs(genus_dir, exist_ok=True)
        
        genus_df = pd.DataFrame(genus_data)
        columns_order = ['id', 'description', 'genus', 'original_genus', 'species', 'primer_count']
        for primer_name in primer_pairs.keys():
            columns_order.extend([
                f'{primer_name}_present',
                f'{primer_name}_fwd_positions',
                f'{primer_name}_rev_positions'
            ])
        
        genus_df = genus_df[columns_order]
        column_rename = {}
        for primer in primer_pairs:
            column_rename[f'{primer}_present'] = f'{primer} (1=present)'
            column_rename[f'{primer}_fwd_positions'] = f'{primer} forward positions'
            column_rename[f'{primer}_rev_positions'] = f'{primer} reverse positions'
        
        genus_df = genus_df.rename(columns=column_rename)
        genus_df.to_excel(os.path.join(genus_dir, f"{genus}_detailed_results.xlsx"), index=False)
        
        summary_df = create_genus_summary(genus_data, genus, cbhpaca_counts[genus])
        summary_df.to_excel(os.path.join(genus_dir, f"{genus}_summary.xlsx"), index=False)
        save_percentage_calculation(genus, summary_df, genus_dir)
        create_genus_barplot(genus_data, genus, genus_dir)
    
    top_abundant_primers = determine_top_abundant_primers(output_dir)
    top_primers_dir = os.path.join(output_dir, "top_abundant_primers")
    os.makedirs(top_primers_dir, exist_ok=True)
    
    with open(os.path.join(top_primers_dir, "top_abundant_primers.txt"), 'w') as f:
        f.write("Top Abundant Primers for Each Genus\n")
        f.write("="*50 + "\n\n")
        for genus, primers in top_abundant_primers.items():
            f.write(f"{genus}: {', '.join(primers)}\n")
    
    no_primer_species, missing_all_top_primers_species = identify_special_species(all_results, top_abundant_primers)
    special_dir = os.path.join(output_dir, "special_species")
    os.makedirs(special_dir, exist_ok=True)
    
    if no_primer_species:
        pd.DataFrame(no_primer_species).to_excel(
            os.path.join(special_dir, "species_with_no_primers.xlsx"), index=False)
    else:
        with open(os.path.join(special_dir, "no_species_without_primers.txt"), 'w') as f:
            f.write("All species had at least one primer pair match.")
    
    for genus, species_list in missing_all_top_primers_species.items():
        if species_list:
            df = pd.DataFrame(species_list)
            df = df[['Species', 'Assigned Genus', 'Original Genus', 'Missing Top Primers', 'Other Primers Present', 'ID', 'Description']]
            df.to_excel(
                os.path.join(special_dir, f"{genus}_species_missing_all_top_primers.xlsx"), index=False)
    
    all_missing_species = []
    for species_list in missing_all_top_primers_species.values():
        all_missing_species.extend(species_list)
    
    if all_missing_species:
        pd.DataFrame(all_missing_species).to_excel(
            os.path.join(special_dir, "all_species_missing_top_primers.xlsx"), index=False)
    else:
        with open(os.path.join(special_dir, "no_species_missing_top_primers.txt"), 'w') as f:
            f.write("All species had at least one of their genus's top primers.")
    
    create_genus_primer_heatmap(all_results, cbhpaca_counts, output_dir)
    print("Analysis completed successfully!")

if __name__ == "__main__":
    input_fasta = "primer_search_results_code1/pathogenic_fungi_sequences.fasta"
    output_directory = "primer_analysis_results555"
    main(input_fasta, output_directory)
