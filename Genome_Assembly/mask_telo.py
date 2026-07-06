import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq

def mask_global_telomeres(input_fasta, output_fasta, motif="TTTAGGG", min_repeats=3):
    rev_comp = str(Seq(motif).reverse_complement())
    
    # Create regex to find tandem arrays anywhere in the sequence
    # e.g., (TTTAGGG){3,} matches 3 or more consecutive motifs
    fwd_pattern = re.compile(f"({motif}){{{min_repeats},}}", re.IGNORECASE)
    rev_pattern = re.compile(f"({rev_comp}){{{min_repeats},}}", re.IGNORECASE)
    
    # Process chromosome by chromosome to maintain low RAM usage
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq)
            original_len = len(seq)
            
            # Use a lambda function to dynamically measure each match 
            # and replace it with the exact same number of Ns.
            # This guarantees coordinates will never shift.
            masked_seq = fwd_pattern.sub(lambda m: "N" * len(m.group(0)), seq)
            masked_seq = rev_pattern.sub(lambda m: "N" * len(m.group(0)), masked_seq)
            
            # Calculate exactly how many bases were newly masked for reporting
            masked_count = masked_seq.count('N') - seq.count('N')
            
            record.seq = Seq(masked_seq)
            SeqIO.write(record, out_f, "fasta")
            
            print(f"{record.id}: Masked {masked_count} bp of telomeric repeats. Total length: {original_len}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python mask_all_telomeres.py <input.fasta> <output.fasta>")
        sys.exit(1)
        
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    mask_global_telomeres(input_file, output_file)
