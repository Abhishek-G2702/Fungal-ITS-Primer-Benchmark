import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def main():
    try:
        # --- 1. Load ITS sequence from FASTA ---
        fasta_file = "ITS_reference.fasta"
        try:
            record = SeqIO.read(fasta_file, "fasta")
            sequence = str(record.seq).upper()
            seq_length = len(sequence)
            print(f"Loaded reference sequence '{record.id}' with length: {seq_length} bp")
        except Exception as e:
            print(f"Error loading FASTA file: {e}")
            return

        # --- 2. Load primers from Excel ---
        excel_file = "primer_ITS.xlsx"
        try:
            df = pd.read_excel(excel_file)
            print("\nLoaded primer data with columns:", df.columns.tolist())
            
            # Clean up primer names and sequences
            df['Primer'] = df['Primer'].str.strip()
            df["Forward_primer(5'-3')"] = df["Forward_primer(5'-3')"].astype(str).str.strip().str.upper()
            df["Reverse_primer(5'-3')"] = df["Reverse_primer(5'-3')"].astype(str).str.strip().str.upper()
        except Exception as e:
            print(f"Error loading Excel file: {e}")
            return

        # --- 3. Find primer binding positions ---
        positions = []
        primer_pairs = {}

        for i, row in df.iterrows():
            try:
                primer_name = row['Primer']
                fwd_seq = row["Forward_primer(5'-3')"]
                rev_seq = row["Reverse_primer(5'-3')"]
                
                if pd.isna(fwd_seq) or pd.isna(rev_seq) or fwd_seq == "NAN" or rev_seq == "NAN":
                    print(f"Skipping {primer_name} due to missing primers")
                    continue
                
                reverse_primer_name = row['Primer.1'] if 'Primer.1' in df.columns else f"Rev_{primer_name}"
                primer_pairs[primer_name] = reverse_primer_name
                
                fwd_start = sequence.find(fwd_seq)
                if fwd_start != -1:
                    fwd_end = fwd_start + len(fwd_seq)
                    rev_rc = str(Seq(rev_seq).reverse_complement())
                    rev_start = sequence.find(rev_rc)
                    
                    if rev_start != -1:
                        rev_end = rev_start + len(rev_rc)
                        positions.append((primer_name, reverse_primer_name, fwd_start, fwd_end, rev_start, rev_end))
                        print(f"Found positions for {primer_name}/{reverse_primer_name}: Fwd({fwd_start}-{fwd_end}), Rev({rev_start}-{rev_end})")
                    else:
                        print(f"Could not find reverse primer for {primer_name} in sequence")
                else:
                    print(f"Could not find forward primer for {primer_name} in sequence")
            except Exception as e:
                print(f"Error processing {row['Primer']}: {e}")
                continue

        if not positions:
            print("\nNo valid primer binding positions found!")
            return

        # --- 4. Save positions to Excel ---
        output_data = []
        for fwd_name, rev_name, f_start, f_end, r_start, r_end in positions:
            output_data.append([fwd_name, rev_name, f_start, f_end, r_start, r_end, r_end - f_start])

        output_df = pd.DataFrame(output_data, 
                               columns=["Forward Primer", "Reverse Primer", "Forward Start", "Forward End", 
                                       "Reverse Start", "Reverse End", "Amplicon Length"])
        output_df.to_excel("ITS_primer_positions_output.xlsx", index=False)
        print("\nSaved results to ITS_primer_positions_output.xlsx")

        # --- 5. Visualization with ITS regions at bottom ---
        fig, (ax, ax_regions) = plt.subplots(
            2, 1, 
            figsize=(16, 8),
            gridspec_kw={'height_ratios': [3, 1]},  # Primer plot takes 3/4 space, regions take 1/4
            sharex=True
        )
        
        # Define ITS regions (adjust these based on your sequence)
        its1_start = 0
        its1_end = 400  # End of ITS1
        five8s_start = its1_end
        five8s_end = five8s_start + 160  # 5.8S is ~160bp
        its2_start = five8s_end
        its2_end = 800  # End of ITS2
        lsu_start = its2_end  # Start of 28S
        
        regions = {
            "18S": (0, its1_start),
            "ITS1": (its1_start, five8s_start),
            "5.8S": (five8s_start, its2_start),
            "ITS2": (its2_start, its2_end),
            "28S": (its2_end, seq_length)
        }
        
        # Create a color palette for primer pairs
        colors = plt.cm.tab20.colors
        
        # --- Plot primer pairs on top axes ---
        for i, (fwd_name, rev_name, f_start, f_end, r_start, r_end) in enumerate(positions):
            y_pos = i + 1
            color = colors[i % len(colors)]
            
            # Amplicon region
            ax.add_patch(mpatches.Rectangle(
                (f_start, y_pos - 0.3), 
                r_end - f_start, 
                0.6, 
                color=color, 
                alpha=0.4,
                linewidth=1
            ))
            
            # Connecting line
            ax.plot([f_start, r_end], [y_pos, y_pos], color='black', linewidth=1.5, linestyle='-')
            
            # Forward primer (green arrow)
            ax.plot(
                f_start, y_pos, 
                '>', 
                color='green', 
                markersize=14,
                markeredgecolor='black',
                markeredgewidth=0.7,
                label='Forward Primer' if i == 0 else ""
            )
            
            # Reverse primer (blue arrow)
            ax.plot(
                r_end, y_pos, 
                '<', 
                color='blue', 
                markersize=14,
                markeredgecolor='black',
                markeredgewidth=0.7,
                label='Reverse Primer' if i == 0 else ""
            )
            
            # Primer pair label
            ax.text(
                (f_start + r_end)/2, y_pos + 0.25, 
                f"{fwd_name}/{rev_name}", 
                ha='center', 
                va='center',
                fontsize=10,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=2)
            )
            
            # Position markers
            ax.text(
                f_start, y_pos - 0.25, 
                f"{f_start}", 
                ha='center', 
                va='center',
                color='green',
                fontsize=9,
                fontweight='bold'
            )
            
            ax.text(
                r_end, y_pos - 0.25, 
                f"{r_end}", 
                ha='center', 
                va='center',
                color='blue',
                fontsize=9,
                fontweight='bold'
            )
            
            # Amplicon length label
            amp_len = r_end - f_start
            ax.text(
                (f_start + r_end)/2, y_pos - 0.25, 
                f"{amp_len} bp", 
                ha='center', 
                va='center',
                fontsize=9,
                fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1)
            )

        # --- Plot ITS regions on bottom axes ---
        region_colors = {
            "18S": "#FFDDC1",  # Light orange
            "ITS1": "#C1FFC1",  # Light green
            "5.8S": "#FFC1FF",  # Light pink
            "ITS2": "#C1FFFF",  # Light cyan
            "28S": "#E6E6FA"    # Light purple
        }
        
        # Draw region bars
        for region, (start, end) in regions.items():
            if end <= start:  # Skip zero-length regions
                continue
                
            ax_regions.add_patch(mpatches.Rectangle(
                (start, 0), 
                end - start, 
                1, 
                color=region_colors[region], 
                alpha=0.7,
                linewidth=1
            ))
            
            # Add region label - now including all regions with significant length
            if (end - start) > 0:  # Only label if region has length
                ax_regions.text(
                    (start + end)/2, 
                    0.5, 
                    f"{region}\n({end-start} bp)",
                    ha='center',
                    va='center',
                    fontsize=10,
                    fontweight='bold'
                )
            
            # Add region boundary lines
            if region != "18S":  # Skip leftmost boundary
                ax_regions.axvline(x=start, color='black', linestyle=':', alpha=0.7, linewidth=1)

        # --- Configure axes ---
        # Top axes (primers)
        ax.set_xlim(0, seq_length)
        ax.set_ylim(0.5, len(positions) + 1)
        ax.set_ylabel("Primer Pairs", fontsize=12, labelpad=15)
        ax.set_yticks([])
        ax.grid(True, axis='x', linestyle=':', alpha=0.5)
        ax.legend(loc='upper right', fontsize=10)
        
        # Bottom axes (regions)
        ax_regions.set_ylim(0, 1)
        ax_regions.set_yticks([])
        ax_regions.set_ylabel("ITS Regions", fontsize=12, labelpad=15)
        ax_regions.set_xlabel("Nucleotide Position in ITS Reference Sequence", fontsize=12, labelpad=10)
        ax_regions.grid(False)
        
        # Title
        fig.suptitle(
            f"Primer Binding Sites on {record.id}\nReference Sequence Length: {seq_length} bp", 
            fontsize=14, 
            y=1.02,
            fontweight='bold'
        )
        
        plt.tight_layout()
        
        # Save high-resolution image
        plt.savefig("ITS_primer_map_with_bottom_regions.png", dpi=300, bbox_inches='tight')
        print("\nSaved visualization with bottom ITS regions to ITS_primer_map_with_bottom_regions.png")
        plt.show()

    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
