#!/bin/bash
#SBATCH --job-name=blast_loop
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=2:00:00
#SBATCH --output=blast_loop.%j.out
#SBATCH --error=blast_loop.%j.err
ml BLAST+/2.14.1-gompi-2023a
ml SAMtools/1.23.1-GCC-13.3.0
ml Miniforge3/24.1.2-0

# Create the Python script on the fly
cat << 'EOF' > run_sv_mapping.py
import subprocess
import os
import csv

# ==========================================
PLUG YOUR PATHS HERE
# ==========================================
C2_FASTA = "/scratch/kbb45638/telomere_final_analysis/reference/ABSxc2v2.1.fasta"
W22_BLASTDB = "/scratch/kbb45638/telomere_final_analysis/blast_w22/w22_db"
CHROMOSOME = "CM039153.1_RagTag"
OUTPUT_CSV = "w22_mapped_coordinates.csv"

# Temporary files generated during the run
QUERY_FASTA = "temp_queries.fasta"
BLAST_OUT = "temp_blast_results.tsv"

# Hardcoded data
SV_DATA = [
    ("4a1", "duplication", 114221395),
    ("4a1", "duplication", 114993680),
    ("4b2", "telomere", 150412744),
    ("4b2", "foldback", 134947803),
    ("4a3.d6", "telomere", 171651351),
    ("4a3.d6", "foldback", 172771714),
    ("4a3.d3", "telomere", 114937974),
    ("4a3.d3", "foldback", 114926723),
    ("4a3.d5", "telomere", 162247299),
    ("4a3.d4", "telomere", 149915150),
    ("4a3.d2", "telomere", 152665977),
    ("4a3.d1", "telomere", 141286410),
    ("4a8", "telomere", 151407899),
    ("4a9", "telomere", 163424315)
]

def extract_sequences():
    print("Extracting 20kb windows from C2 genome...")
    with open(QUERY_FASTA, "w") as out_fasta:
        for sample, sv_type, location in SV_DATA:
            start = location - 10000
            end = location + 10000
            query_id = f"{sample}|{sv_type}|{location}"
            region = f"{CHROMOSOME}:{start}-{end}"
            
            cmd = ["samtools", "faidx", C2_FASTA, region]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"Error extracting {region}: {result.stderr}")
                continue
            
            fasta_lines = result.stdout.strip().split('\n')
            if fasta_lines:
                fasta_lines[0] = f">{query_id}"
                out_fasta.write('\n'.join(fasta_lines) + '\n')

def run_blast():
    print("Running BLAST against W22 database...")
    threads = os.environ.get("SLURM_CPUS_PER_TASK", "1")
    cmd = [
        "blastn",
        "-query", QUERY_FASTA,
        "-db", W22_BLASTDB,
        "-outfmt", "6",
        "-out", BLAST_OUT,
        "-num_threads", threads
    ]
    subprocess.run(cmd, check=True)

def parse_blast_and_calculate():
    print("Parsing BLAST results and calculating W22 coordinates...")
    mapped_queries = set()
    results_table = []
    
    with open(BLAST_OUT, "r") as blast_file:
        for line in blast_file:
            cols = line.strip().split('\t')
            if len(cols) < 12: continue
                
            query_id = cols[0]
            if query_id in mapped_queries: continue
            
            sseqid = cols[1]
            qstart, qend = int(cols[6]), int(cols[7])
            sstart, send = int(cols[8]), int(cols[9])
            
            BREAKPOINT_REL_POS = 10001
            
            if qstart <= BREAKPOINT_REL_POS <= qend:
                offset = BREAKPOINT_REL_POS - qstart
                w22_location = (sstart + offset) if (sstart < send) else (sstart - offset)
                
                sample, sv_type, c2_loc = query_id.split('|')
                results_table.append({
                    "samples": sample,
                    "sv_type": sv_type,
                    "c2_location": c2_loc,
                    "w22_chromosome": sseqid,
                    "w22_location": w22_location
                })
                mapped_queries.add(query_id)

    for sample, sv_type, location in SV_DATA:
        query_id = f"{sample}|{sv_type}|{location}"
        if query_id not in mapped_queries:
            results_table.append({
                "samples": sample,
                "sv_type": sv_type,
                "c2_location": str(location),
                "w22_chromosome": "UNMAPPED",
                "w22_location": "UNMAPPED"
            })
            
    return results_table

def write_output(results_table):
    print(f"Writing final mapped coordinates to {OUTPUT_CSV}...")
    fieldnames = ["samples", "sv_type", "c2_location", "w22_chromosome", "w22_location"]
    with open(OUTPUT_CSV, "w", newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results_table:
            writer.writerow(row)
            
    if os.path.exists(QUERY_FASTA): os.remove(QUERY_FASTA)
    if os.path.exists(BLAST_OUT): os.remove(BLAST_OUT)
    print("Grind complete. Data is ready.")

if __name__ == "__main__":
    extract_sequences()
    run_blast()
    results = parse_blast_and_calculate()
    write_output(results)
EOF

# Execute the generated Python script
python run_sv_mapping.py
