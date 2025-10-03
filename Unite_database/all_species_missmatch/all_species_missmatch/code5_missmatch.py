import os
from Bio import SeqIO
import pandas as pd
from multiprocessing import Pool, cpu_count
from collections import defaultdict

# ====================== CONFIGURATION ======================
FASTA_FILE = "UNITE_database_fasta_file.fasta"  # Your input file
OUTPUT_DIR = "code55_Primer_Results"       # Main output directory
MAX_MISMATCHES = 2                  # Allowed mismatches
STRICT_3PRIME = True                # Reject 3'-end mismatches
N_PROCESSES = max(1, cpu_count() - 1)

# Primer pairs (add/remove as needed)
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

# ====================== CORE FUNCTIONS ======================
def parse_fasta_header(description):
    """Extract accession number and species name from FASTA header"""
    parts = description.split('|')
    species_name = parts[0].split('>')[-1]  # Remove '>' if present
    accession = parts[1] if len(parts) > 1 else "N/A"
    return accession, species_name

def approximate_search(primer, sequence):
    """Find primer matches with mismatch tolerance"""
    len_primer = len(primer)
    matches = []
    
    for i in range(len(sequence) - len_primer + 1):
        mismatches = 0
        details = []
        bad_3prime = False
        
        for j in range(len_primer):
            if primer[j] != sequence[i + j]:
                mismatches += 1
                details.append(f"{primer[j]}→{sequence[i+j]}@pos{j+1}")
                
                # Check 3' end (last 3 bases)
                if STRICT_3PRIME and j >= (len_primer - 3):
                    bad_3prime = True
                
                if mismatches > MAX_MISMATCHES:
                    break
        
        if mismatches <= MAX_MISMATCHES:
            matches.append({
                'start': i + 1,  # 1-based position
                'end': i + len_primer,
                'mismatches': mismatches,
                'details': "; ".join(details),
                'bad_3prime': bad_3prime,
                'multi_bind': False  # Updated later
            })
    
    # Flag multi-binding primers
    if len(matches) > 1:
        for m in matches:
            m['multi_bind'] = True
    
    return matches

def process_record(record, primer_name, primers):
    """Process one sequence for a specific primer pair"""
    seq = str(record.seq).upper()
    rev_seq = str(record.seq.reverse_complement()).upper()
    
    # Extract accession and species name
    accession, species_name = parse_fasta_header(record.description)
    
    # Forward primer
    fwd_matches = approximate_search(primers["Forward"], seq)
    valid_fwd = [m for m in fwd_matches if not (STRICT_3PRIME and m['bad_3prime'])]
    
    # Reverse primer (with position adjustment)
    rev_matches = approximate_search(primers["Reverse"], rev_seq)
    valid_rev = [m for m in rev_matches if not (STRICT_3PRIME and m['bad_3prime'])]
    adj_positions = [
        f"{len(seq) - m['end'] + 1}-{len(seq) - m['start'] + 1}"
        for m in valid_rev
    ]
    
    # Amplicon size calculation
    if valid_fwd and valid_rev:
        amp_size = (len(seq) - valid_rev[0]['end'] + 1) - valid_fwd[0]['start'] + len(primers["Reverse"])
        pair_present = "Yes"
    else:
        amp_size = "N/A"
        pair_present = "No"
    
    return {
        "Accession": accession,
        "Species_Name": species_name,
        "Sequence_Length": len(seq),
        # Forward primer results
        "Fwd_Present": "Yes" if valid_fwd else "No",
        "Fwd_Positions": "; ".join(f"{m['start']}-{m['end']}" for m in valid_fwd) or "None",
        "Fwd_Mismatches": valid_fwd[0]["mismatches"] if valid_fwd else 0,
        "Fwd_Mismatch_Details": valid_fwd[0]["details"] if valid_fwd else "None",
        "Fwd_3prime_Issue": "Yes" if any(m['bad_3prime'] for m in fwd_matches) else "No",
        "Fwd_Multi_Bind": "Yes" if any(m['multi_bind'] for m in fwd_matches) else "No",
        # Reverse primer results
        "Rev_Present": "Yes" if valid_rev else "No",
        "Rev_Positions": "; ".join(adj_positions) or "None",
        "Rev_Mismatches": valid_rev[0]["mismatches"] if valid_rev else 0,
        "Rev_Mismatch_Details": valid_rev[0]["details"] if valid_rev else "None",
        "Rev_3prime_Issue": "Yes" if any(m['bad_3prime'] for m in rev_matches) else "No",
        "Rev_Multi_Bind": "Yes" if any(m['multi_bind'] for m in rev_matches) else "No",
        # Amplicon info
        "Amplicon_Size": amp_size,
        "Primer_Pair_Present": pair_present,
        "Amplicon_Valid": "Yes" if (valid_fwd and valid_rev) else "No"
    }

# ====================== OUTPUT HANDLING ======================
def save_primer_results(primer_name, results):
    """Save results for one primer pair"""
    primer_dir = os.path.join(OUTPUT_DIR, primer_name)
    os.makedirs(primer_dir, exist_ok=True)
    
    df = pd.DataFrame(results)
    # Reorder columns to make important info first
    cols = ["Accession", "Species_Name", "Sequence_Length", "Primer_Pair_Present"] + \
           [col for col in df.columns if col not in ["Accession", "Species_Name", "Sequence_Length", "Primer_Pair_Present"]]
    df = df[cols]
    
    output_file = os.path.join(primer_dir, f"{primer_name}_results.xlsx")
    df.to_excel(output_file, index=False)

# ====================== MAIN EXECUTION ======================
def main():
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load sequences
    with open(FASTA_FILE) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    
    # Process each primer pair
    for primer_name, primers in primer_pairs.items():
        print(f"🔍 Processing {primer_name}...")
        
        with Pool(N_PROCESSES) as pool:
            results = pool.starmap(
                process_record,
                [(record, primer_name, primers) for record in records]
            )
        
        save_primer_results(primer_name, results)
    
    print(f"""
✅ Analysis complete!
📂 Results organized in '{OUTPUT_DIR}/' with:
   - Separate folder for each primer pair
   - Excel file per primer with all validation metrics
   - Accession numbers and species names in separate columns
   - 'Primer_Pair_Present' column showing if BOTH primers are found
""")

if __name__ == "__main__":
    main()
