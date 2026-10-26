#!/usr/bin/env python3
"""
Build Ensembl + NR annotation and FASTA files per gene.

This script:
1. Iterates over gene subfolders in a parent directory
2. Finds Ensembl curated FASTA + annotation tables
3. Filters MMseqs2 NR hits
4. Merges Ensembl and NR annotations
5. Writes per-gene annotation TSV and FASTA files

Designed for transparency and reproducibility.
"""

import argparse
import os
import pandas as pd
from Bio import SeqIO


# -------------------------------------------------------------------
# Taxonomic ranks expected in MMseqs 'taxlineage' column
# (e.g. d_Eukaryota;k_Metazoa;p_Chordata;...)
# -------------------------------------------------------------------
TAXONOMIC_RANKS = [
    'domain',
    'kingdom',
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species'
]


# -------------------------------------------------------------------
# Command-line argument parsing
# -------------------------------------------------------------------
def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Build Ensembl + NR annotation and FASTA files per gene"
    )
    parser.add_argument(
        "parent_gene_dir",
        help="Path to parent directory containing gene subfolders (RSAD2, ZNFX1, etc.)"
    )
    parser.add_argument(
        "mmseqs_file",
        help="Path to MMseqs2 TSV results file"
    )
    return parser.parse_args()


# -------------------------------------------------------------------
# FASTA utilities
# -------------------------------------------------------------------
def load_fasta_to_df(fasta_file):
    """
    Read a FASTA file into a DataFrame with protein_id and protein_sequence.
    """
    print(f"    📥 Reading FASTA: {os.path.basename(fasta_file)}")
    df = pd.DataFrame(
        [
            {
                "protein_id": rec.id,
                "protein_sequence": str(rec.seq)
            }
            for rec in SeqIO.parse(fasta_file, "fasta")
        ]
    )
    print(f"       → {df.shape[0]} sequences loaded")
    return df


# -------------------------------------------------------------------
# MMseqs taxonomy parsing
# -------------------------------------------------------------------
def extract_taxonomy(df_mmseqs):
    """
    Extract taxonomic ranks from MMseqs 'taxlineage' column.
    """
    print("🔬 Extracting taxonomy columns from MMseqs results")

    for rank in TAXONOMIC_RANKS:
        prefix = f"{rank[0]}_"
        df_mmseqs[rank] = df_mmseqs["taxlineage"].apply(
            lambda x: next(
                (i[2:] for i in x.split(";") if i.startswith(prefix)),
                None
            ) if isinstance(x, str) else None
        )

    return df_mmseqs


# -------------------------------------------------------------------
# Per-gene processing
# -------------------------------------------------------------------
def process_gene_folder(gene_dir, df_mmseqs):
    """
    Process a single gene folder:
    - load Ensembl curated FASTA
    - load Ensembl annotation table
    - filter MMseqs NR hits
    - combine annotations
    - write TSV + FASTA outputs
    """
    gene = os.path.basename(gene_dir)

    print("\n" + "=" * 70)
    print(f"🧬 Processing gene: {gene}")
    print("=" * 70)

    # Locate required files
    fasta_files = [
        f for f in os.listdir(gene_dir)
        if f.endswith("ensembl_curated.fasta")
    ]
    annotation_files = [
        f for f in os.listdir(gene_dir)
        if f.endswith("ensembl_master_table_annotation.tsv")
    ]

    if not fasta_files or not annotation_files:
        print(f"⚠️  SKIP: required files not found in {gene}")
        return

    fasta_path = os.path.join(gene_dir, fasta_files[0])
    annotation_path = os.path.join(gene_dir, annotation_files[0])

    print(f"📂 FASTA file      : {fasta_files[0]}")
    print(f"📂 Annotation file : {annotation_files[0]}")

    # ---------------------------------------------------------------
    # Load Ensembl curated FASTA
    # ---------------------------------------------------------------
    df_ensembl_curated = load_fasta_to_df(fasta_path)

    # ---------------------------------------------------------------
    # Load Ensembl annotation table
    # ---------------------------------------------------------------
    print("    📥 Reading Ensembl annotation table")
    df_annotation = pd.read_csv(annotation_path, sep="\t")
    print(f"       → {df_annotation.shape[0]} annotation rows")

    # Normalize gene names
    df_annotation["gene_name"] = df_annotation["gene_name"].str.upper()

    # Keep only curated proteins
    df_annotation = df_annotation[
        df_annotation["protein_id"].isin(df_ensembl_curated["protein_id"])
    ]
    print(f"       → {df_annotation.shape[0]} curated Ensembl proteins")

    # Attach protein sequences
    df_annotation = df_annotation.merge(
        df_ensembl_curated,
        on="protein_id",
        how="left"
    )

    # ---------------------------------------------------------------
    # Filter MMseqs NR hits
    # ---------------------------------------------------------------
    print("    🔎 Filtering MMseqs NR hits")

    df_hits = df_mmseqs[
        (df_mmseqs["query"].isin(df_annotation["protein_id"])) &
        (df_mmseqs["species"].notna()) &
        (df_mmseqs["tseq"].str.count("X") <= 3)
    ].copy()

    print(f"       → {df_hits.shape[0]} hits after initial filters")

    # Keep best hit per unique sequence
    df_hits = (
        df_hits
        .sort_values("evalue")
        .drop_duplicates(subset="tseq")
    )
    print(f"       → {df_hits.shape[0]} non-redundant NR sequences")

    # Remove NR sequences identical to Ensembl
    df_hits = df_hits[
        ~df_hits["tseq"].isin(df_annotation["protein_sequence"])
    ]
    print(f"       → {df_hits.shape[0]} NR-only sequences retained")

    # ---------------------------------------------------------------
    # NR annotation formatting
    # ---------------------------------------------------------------
    cols = [
        'protein_id', 'source', 'source_id', 'source_name',
        'gene_id', 'gene_name', 'transcript_id',
        'perc_id', 'domain', 'kingdom', 'phylum',
        'class', 'order', 'family', 'genus',
        'species', 'protein_sequence'
    ]

    df_nr = df_hits.rename(columns={
        "target": "protein_id",
        "query": "source_id",
        "theader": "gene_name",
        "tseq": "protein_sequence"
    })

    df_nr["perc_id"] = df_nr["fident"] * 100
    df_nr["source"] = "nr"

    for c in cols:
        if c not in df_nr.columns:
            df_nr[c] = None

    df_nr = df_nr[cols]

    # ---------------------------------------------------------------
    # Ensembl annotation formatting
    # ---------------------------------------------------------------
    df_ens = df_annotation.copy()
    df_ens["domain"] = "Eukaryota"
    df_ens["kingdom"] = "Metazoa"
    df_ens["source"] = "ensembl"
    df_ens = df_ens[cols]

    # ---------------------------------------------------------------
    # Combine Ensembl + NR
    # ---------------------------------------------------------------
    df_combined = pd.concat([df_ens, df_nr], ignore_index=True)
    print(f"📊 Total sequences in final set: {df_combined.shape[0]}")

    # ---------------------------------------------------------------
    # Taxonomic distribution summary
    # ---------------------------------------------------------------
    # ---- Domain-level distribution ----
    print("\n🧫 Hits per domain:")
    domain_counts = (
        df_combined["domain"]
        .fillna("Unclassified")
        .value_counts()
    )
    for domain, count in domain_counts.items():
        print(f"    {domain:15s} : {count}")

    # ---- Kingdoms within Eukaryota ----
    df_euk = df_combined[df_combined["domain"] == "Eukaryota"]

    if not df_euk.empty:
        print("\n🌿 Kingdoms within Eukaryota:")
        kingdom_counts = (
            df_euk["kingdom"]
            .fillna("Unclassified")
            .value_counts()
        )
        for kingdom, count in kingdom_counts.items():
            print(f"    {kingdom:15s} : {count}")
    else:
        print("\n🌿 No Eukaryota hits found")

    # ---- Phyla within Metazoa ----
    df_metazoa = df_euk[df_euk["kingdom"] == "Metazoa"]

    if not df_metazoa.empty:
        print("\n🦴 Phyla within Metazoa:")
        phylum_counts = (
            df_metazoa["phylum"]
            .fillna("Unclassified")
            .value_counts()
        )
        for phylum, count in phylum_counts.items():
            print(f"    {phylum:20s} : {count}")
    else:
        print("\n🦴 No Metazoa hits found")

    # ---------------------------------------------------------------
    # Write outputs
    # ---------------------------------------------------------------
    annotation_out = os.path.join(
        gene_dir,
        f"{gene}_ensembl_nr_annotation.tsv"
    )
    df_combined.to_csv(annotation_out, sep="\t", index=False)
    print(f"\n📝 Wrote annotation TSV: {os.path.basename(annotation_out)}")

    fasta_out = os.path.join(
        gene_dir,
        f"{df_combined.shape[0]}_{gene}_ensembl_nr_hits.fasta"
    )

    with open(fasta_out, "w") as f:
        for _, row in df_combined.iterrows():
            f.write(f">{row['protein_id']}\n{row['protein_sequence']}\n")

    print(f"🧬 Wrote FASTA file     : {os.path.basename(fasta_out)}")
    print(f"✅ Finished {gene}")


# -------------------------------------------------------------------
# Main entry point
# -------------------------------------------------------------------
def main():
    args = parse_args()

    print("\n🚀 Starting Ensembl + NR annotation pipeline")
    print(f"📁 Parent gene directory : {args.parent_gene_dir}")
    print(f"📄 MMseqs results file   : {args.mmseqs_file}\n")

    # Load MMseqs results once
    df_mmseqs = pd.read_csv(args.mmseqs_file, sep="\t")
    print(f"📥 Loaded MMseqs table with {df_mmseqs.shape[0]} rows")

    df_mmseqs = extract_taxonomy(df_mmseqs)

    # Process each gene folder
    for entry in sorted(os.listdir(args.parent_gene_dir)):
        gene_dir = os.path.join(args.parent_gene_dir, entry)
        if os.path.isdir(gene_dir):
            process_gene_folder(gene_dir, df_mmseqs)

    print("\n🎉 All genes processed successfully")


if __name__ == "__main__":
    main()