import pandas as pd
import pyranges as pr

# ==========================================
# USER CONFIGURATION
# ==========================================
MAPPED_CSV = "w22_mapped_coordinates.csv"
GENE_GFF3 = "/Users/kempton/Desktop/Telomere_Paper/breaksite_enrichment/w22/Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.gff3"
TE_GFF3 = "/Users/kempton/Desktop/Telomere_Paper/breaksite_enrichment/w22/W22.structuralTEv2.disjoined.2018-09-22.gff3"
OUTPUT_CSV = "w22_feature_annotations.csv"

def load_annotations():
    print("Loading gene and TE annotations...")
    
    genes_pr = pr.read_gff3(GENE_GFF3)
    # CRITICAL FIX: Only look at actual genes, ignoring 'chromosome', 'mRNA', 'exon' lines
    genes_pr = genes_pr[genes_pr.Feature == "gene"]
    
    tes_pr = pr.read_gff3(TE_GFF3)
    te_df = tes_pr.df
    
    def fix_chrom(c):
        c_str = str(c)
        if c_str.isdigit():
            return f"chr{c_str}"
        return c_str
        
    te_df["Chromosome"] = te_df["Chromosome"].apply(fix_chrom)
    fixed_tes_pr = pr.PyRanges(te_df)
    
    return genes_pr, fixed_tes_pr

def annotate_variants():
    print("Parsing mapped structural variant coordinates...")
    df_sv = pd.read_csv(MAPPED_CSV)
    
    df_mapped = df_sv[df_sv["w22_location"] != "UNMAPPED"].copy()
    df_mapped["w22_location"] = df_mapped["w22_location"].astype(int)
    
    df_mapped["Start"] = df_mapped["w22_location"] - 1
    df_mapped["End"] = df_mapped["w22_location"]
    
    sv_pr = pr.PyRanges(df_mapped.rename(columns={"w22_chromosome": "Chromosome"}))
    
    genes_pr, tes_pr = load_annotations()
    
    print("Intersecting variants with genes...")
    gene_overlaps = sv_pr.join(genes_pr, suffix="_gene")
    gene_df = gene_overlaps.df
    
    gene_hit_dict = {}
    if not gene_df.empty:
        for _, row in gene_df.iterrows():
            key = (str(row["samples"]), str(row["c2_location"]))
            gene_id = row.get("ID", row.get("Name", "known_gene"))
            gene_hit_dict[key] = f"Gene ({gene_id})"

    print("Intersecting variants with TEs...")
    te_overlaps = sv_pr.join(tes_pr, suffix="_te")
    te_df = te_overlaps.df
    
    te_hit_dict = {}
    if not te_df.empty:
        for _, row in te_df.iterrows():
            key = (str(row["samples"]), str(row["c2_location"]))
            te_family = row.get("sup", row.get("ID", "TE")) 
            
            # Check the intact status column
            intact_val = str(row.get("intact", "FALSE")).upper()
            is_intact = (intact_val == "TRUE")
            
            # Prioritize Intact over Fragmented if a breakpoint overlaps multiple TE annotations
            if key not in te_hit_dict or is_intact:
                if is_intact:
                    te_hit_dict[key] = f"TE Intact ({te_family})"
                elif key not in te_hit_dict:
                    te_hit_dict[key] = f"TE Fragmented ({te_family})"

    print("Classifying genomic context...")
    final_results = []
    for _, row in df_sv.iterrows():
        sample = str(row["samples"])
        c2_loc = str(row["c2_location"])
        key = (sample, c2_loc)
        
        if row["w22_location"] == "UNMAPPED":
            context = "UNMAPPED"
        else:
            if key in gene_hit_dict:
                context = gene_hit_dict[key]
            elif key in te_hit_dict:
                # Merge fragmented TEs into Intergenic space per the biological rationale
                if "Fragmented" in te_hit_dict[key]:
                    context = f"Intergenic (Degraded TE: {te_hit_dict[key].split('(')[1]}"
                else:
                    context = te_hit_dict[key]
            else:
                context = "Intergenic"
                
        row_dict = row.to_dict()
        row_dict["genomic_context"] = context
        final_results.append(row_dict)
        
    print(f"Writing annotated results to {OUTPUT_CSV}...")
    out_df = pd.DataFrame(final_results)
    out_df.to_csv(OUTPUT_CSV, index=False)
    print("Annotation grind complete.")

if __name__ == "__main__":
    annotate_variants()