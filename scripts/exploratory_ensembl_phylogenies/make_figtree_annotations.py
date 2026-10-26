# Usage: python make_figtree_annotations.py path/to/main_folde

import os
import argparse
import pandas as pd
from ete3 import NCBITaxa

# Initialize NCBI taxonomy DB (first run downloads it locally)
ncbi = NCBITaxa()
# ncbi.update_taxonomy_database()  # uncomment if you need to refresh

def get_taxonomy(species_list):
    """Return DataFrame with taxonomy for a list of species names."""
    target_ranks = ["phylum", "class", "order", "family", "genus"]
    results = []

    for sp in species_list:
        try:
            name2taxid = ncbi.get_name_translator([sp])
            if not name2taxid:
                print(f"❌ Not found: {sp}")
                continue
            taxid = name2taxid[sp][0]

            lineage = ncbi.get_lineage(taxid)
            ranks = ncbi.get_rank(lineage)
            names = ncbi.get_taxid_translator(lineage)

            extracted = {rank: None for rank in target_ranks}
            for tid in lineage:
                rank = ranks.get(tid, "")
                if rank in extracted:
                    extracted[rank] = names[tid]

            extracted["species"] = sp
            results.append(extracted)

        except Exception as e:
            print(f"Error processing {sp}: {e}")

    return pd.DataFrame(results, columns=target_ranks + ["species"])


def process_master_table(master_file):
    """Process one master_table.csv and generate annotation file."""
    df_master = pd.read_csv(master_file)

    # normalize species names from e.g. Homo_sapiens_male -> "Homo sapiens"
    df_master['species'] = df_master['species'].apply(
        lambda x: ' '.join(x.split('_')[:2]).capitalize()
    )

    # make annotation subset
    df_annotation = df_master[[
        'protein_id', 'source_id', 'source_name', 'gene_id',
        'gene_name', 'transcript_id', 'perc_id', 'species'
    ]].copy()

    # add taxonomy
    df = get_taxonomy(df_master['species'].unique())
    df_annotation = pd.merge(df_annotation, df, on='species', how='left')

    # save annotation file
    base = os.path.splitext(os.path.basename(master_file))[0]  # without .csv
    out_file = os.path.join(os.path.dirname(master_file), f"{base}_annotation.tsv")
    df_annotation.to_csv(out_file, sep='\t', index=False)

    print(f"✅ Saved annotation: {out_file}")


def main():
    parser = argparse.ArgumentParser(description="Generate FigTree annotation files from master tables.")
    parser.add_argument("folder", help="Input folder to search for master_table.csv files")
    args = parser.parse_args()

    # walk through folder and subfolders
    for root, _, files in os.walk(args.folder):
        for f in files:
            if f.endswith("master_table.csv"):
                master_file = os.path.join(root, f)
                print(f"Processing: {master_file}")
                process_master_table(master_file)


if __name__ == "__main__":
    main()