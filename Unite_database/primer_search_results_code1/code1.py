import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# =============================================
# CONFIGURATION
# =============================================
INPUT_FASTA = "UNITE_database_fasta_file.fasta"
PRIMER_EXCEL = "primer_ITS.xlsx"
OUTPUT_FOLDER = "primer_search_results"
TARGET_PATHOGENS = [
    "Cryptococcus", "Blastomyces",
    "Histoplasma", "Aspergillus", "Candida"
]

# =============================================
# INITIALIZATION
# =============================================
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# Load primer pairs
df_primers = pd.read_excel(PRIMER_EXCEL, sheet_name=0, usecols="A:D")
df_primers.columns = ["Forward_Name", "Forward_Seq", "Reverse_Name", "Reverse_Seq"]
df_primers.dropna(inplace=True)
df_primers["Forward_Seq"] = df_primers["Forward_Seq"].str.upper()
df_primers["Reverse_Seq"] = df_primers["Reverse_Seq"].str.upper()

# =============================================
# CORE FUNCTIONS
# =============================================
def find_primer_positions(seq, primer_seq, seq_length):
    """Find primer matches in both strands with correct position reporting"""
    rev_seq = str(Seq(seq).reverse_complement())
    
    # Forward strand matches (1-based)
    fwd_matches = [f"{m.start()+1}-{m.end()}" for m in re.finditer(primer_seq, seq)]
    
    # Reverse complement matches (converted to original coordinates)
    rc_matches = [
        f"{seq_length - m.end() + 1}-{seq_length - m.start()} (RevComp)" 
        for m in re.finditer(primer_seq, rev_seq)
    ]
    
    return fwd_matches + rc_matches

def process_sequence(record):
    """Process a single sequence record"""
    seq = str(record.seq).upper()
    seq_length = len(seq)
    desc = record.description
    
    # Extract metadata
    accession = desc.split("|")[1] if "|" in desc else "Unknown"
    species = desc.split("|")[0].replace(">", "").replace("_", " ").strip()
    is_pathogen = any(p.lower() in species.lower() for p in TARGET_PATHOGENS)
    
    # Initialize result rows
    indiv_row = {"Accession": accession, "Species": species}
    pair_row = {"Accession": accession, "Species": species}
    
    # Track hits for this sequence
    hits_for_sequence = {
        "is_pathogen": is_pathogen,
        "record": record,
        "individual": indiv_row,
        "pair": pair_row,
        "primer_pairs_found": set()
    }
    
    return hits_for_sequence, seq, seq_length

# =============================================
# MAIN PROCESSING
# =============================================
all_results = {
    "individual": [],
    "pair": [],
    "pathogen_sequences": [],
    "primer_pair_counts": {"all": {}, "pathogens": {}},
    "species_counts": {"all": set(), "pathogens": set()}
}

with open(INPUT_FASTA, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        hits, seq, seq_length = process_sequence(record)
        
        for _, row in df_primers.iterrows():
            fname = row["Forward_Name"]
            fseq = row["Forward_Seq"]
            rname = row["Reverse_Name"]
            rseq = row["Reverse_Seq"]
            pair_label = f"{fname} + {rname}"
            
            # Find primer positions
            fwd_positions = find_primer_positions(seq, fseq, seq_length)
            rev_positions = find_primer_positions(seq, rseq, seq_length)
            
            # Record individual primer results
            hits["individual"][fname] = len(fwd_positions)
            hits["individual"][f"{fname}_Positions"] = "; ".join(fwd_positions) or "None"
            hits["individual"][rname] = len(rev_positions)
            hits["individual"][f"{rname}_Positions"] = "; ".join(rev_positions) or "None"
            
            # Check for primer pair
            if fwd_positions and rev_positions:
                hits["pair"][pair_label] = "Found"
                hits["pair"][f"{pair_label}_Fwd_Pos"] = "; ".join(fwd_positions)
                hits["pair"][f"{pair_label}_Rev_Pos"] = "; ".join(rev_positions)
                hits["primer_pairs_found"].add(pair_label)
        
        # Update global results
        all_results["individual"].append(hits["individual"])
        all_results["pair"].append(hits["pair"])
        all_results["species_counts"]["all"].add(hits["individual"]["Species"])
        
        if hits["is_pathogen"]:
            all_results["pathogen_sequences"].append(hits["record"])
            all_results["species_counts"]["pathogens"].add(hits["individual"]["Species"])
            for pair in hits["primer_pairs_found"]:
                all_results["primer_pair_counts"]["pathogens"][pair] = (
                    all_results["primer_pair_counts"]["pathogens"].get(pair, 0) + 1
                )
        
        for pair in hits["primer_pairs_found"]:
            all_results["primer_pair_counts"]["all"][pair] = (
                all_results["primer_pair_counts"]["all"].get(pair, 0) + 1
            )

# =============================================
# SAVE RESULTS
# =============================================
# Save individual and pair results
pd.DataFrame(all_results["individual"]).to_excel(
    os.path.join(OUTPUT_FOLDER, "individual_primer_results_all.xlsx"), index=False
)
pd.DataFrame(all_results["pair"]).to_excel(
    os.path.join(OUTPUT_FOLDER, "primer_pair_results_all.xlsx"), index=False
)

# Save pathogenic sequences
with open(os.path.join(OUTPUT_FOLDER, "pathogenic_fungi_sequences.fasta"), "w") as f:
    SeqIO.write(all_results["pathogen_sequences"], f, "fasta")

# Save pathogen-specific results
if all_results["pathogen_sequences"]:
    pathogen_pair_results = [
        r for r in all_results["pair"] 
        if r["Species"] in all_results["species_counts"]["pathogens"]
    ]
    pd.DataFrame(pathogen_pair_results).to_excel(
        os.path.join(OUTPUT_FOLDER, "primer_pair_results_pathogens.xlsx"), index=False
)

# Generate summary statistics
summary_data = []
total_species = len(all_results["species_counts"]["all"])
total_pathogens = len(all_results["species_counts"]["pathogens"])

for pair in df_primers.apply(
    lambda x: f"{x['Forward_Name']} + {x['Reverse_Name']}", axis=1
):
    percent_all = (
        all_results["primer_pair_counts"]["all"].get(pair, 0) / total_species * 100
        if total_species else 0
    )
    percent_pathogen = (
        all_results["primer_pair_counts"]["pathogens"].get(pair, 0) / total_pathogens * 100
        if total_pathogens else 0
    )
    
    summary_data.append({
        "Primer Pair": pair,
        "% Species Found (All)": round(percent_all, 2),
        "% Species Found (Pathogens)": round(percent_pathogen, 2),
        "Total Hits (All)": all_results["primer_pair_counts"]["all"].get(pair, 0),
        "Total Hits (Pathogens)": all_results["primer_pair_counts"]["pathogens"].get(pair, 0)
    })

pd.DataFrame(summary_data).to_excel(
    os.path.join(OUTPUT_FOLDER, "primer_pair_summary.xlsx"), index=False
)

# =============================================
# FINAL OUTPUT
# =============================================
print(f"""
✅ Analysis complete!
============================
Total sequences processed: {len(all_results["individual"])}
Total species: {total_species}
Pathogenic species: {total_pathogens}

Results saved in: {os.path.abspath(OUTPUT_FOLDER)}
""")
