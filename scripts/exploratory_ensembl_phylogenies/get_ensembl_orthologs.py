#!/usr/bin/env python3

# Usage: python get_ensembl_orthologs.py <input.csv>


import pandas as pd
import numpy as np
import requests, sys, os

# ----------------------------
# Helper functions
# ----------------------------

def fetch_endpoint(server, request, content_type):
    """Fetch a single endpoint (GET request)."""
    r = requests.get(server+request, headers={ "Content-Type" : content_type})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    if content_type == 'application/json':
        return r.json()
    else:
        return r.text

def fetch_multiple_endpoints(server, request, content_type, ids):
    """Fetch multiple IDs (POST request)."""
    ids_str = '{"ids":["' + '","'.join(ids) + '"]}'
    r = requests.post(server+request,
                      headers={ "Content-Type" : content_type, "Accept" : content_type},
                      data=ids_str)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    if content_type == 'application/json':
        return r.json()
    else:
        return r.text


# ----------------------------
# Main workflow
# ----------------------------

def main(input_folder):

    input_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith(".csv")]

    for input_file in input_files:

        print(f"[INFO] Reading input file: {input_file}")
        df = pd.read_csv(input_file)

        server = "http://rest.ensembl.org"

        # ---------------- Homology tables ----------------
        print("[INFO] Fetching homology tables from Ensembl...")

        homology_dict = {
            'source_id': [], 'target_id': [], 'perc_id': [], 'perc_pos': [],
            'protein_id': [], 'seq': [], 'species': [], 'taxon_id': [],
            'type': [], 'method_link_type': []
        }

        for index, row in df.iterrows():
            gene_id = row['gene_id']
            ext = f"/homology/id/human/{gene_id}?aligned=0;cigar_line=0"
            try:
                ens_homology = fetch_endpoint(server, ext, "application/json")
            except:
                print(f"[WARNING] Skipping {gene_id} (could not fetch homology info)")
                continue

            for element in ens_homology['data'][0]['homologies']:
                for key in homology_dict.keys():
                    if key in ['type', 'method_link_type']:
                        homology_dict[key].append(element[key])
                    elif key == 'source_id':
                        homology_dict[key].append(element['source']['id'])
                    elif key == 'target_id':
                        homology_dict[key].append(element['target']['id'])
                    else:
                        homology_dict[key].append(element['target'].get(key, np.nan))

        homology_df = pd.DataFrame(homology_dict)
        homology_df_merge = pd.merge(
            homology_df,
            df[['gene_id', 'gene_name']],
            left_on='source_id',
            right_on='gene_id',
            how='inner'
        )
        homology_df_merge.drop(columns=['gene_id'], inplace=True)
        homology_df_merge.rename(columns={'gene_name': 'source_name'}, inplace=True)
        homology_df_merge.rename(columns={'seq': 'protein_seq'}, inplace=True)
        homology_df_merge = homology_df_merge[
            ['source_id', 'source_name', 'target_id', 'species', 'taxon_id',
            'perc_id', 'perc_pos', 'protein_id', 'type', 'method_link_type', 'protein_seq']
        ]

        source_names = sorted(homology_df_merge['source_name'].unique().tolist())
        folder_out = "_".join(source_names)
        os.makedirs(folder_out, exist_ok=True)
        homology_out = "_".join(source_names) + "_ensembl_homology_table.csv"
        homology_df_merge.to_csv(os.path.join(folder_out, homology_out), index=False)
        print(f"[INFO] Homology data written to {os.path.join(folder_out, homology_out)}")

        # ---------------- Extended info ----------------
        print("[INFO] Fetching extended gene/transcript/protein info...")

        df = pd.read_csv(os.path.join(folder_out, homology_out))
        df.drop_duplicates(subset='target_id', inplace=True)
        df.drop(columns=['species', 'taxon_id', 'protein_id', 'protein_seq'], inplace=True)

        isg_dict = {
            'gene_id': [], 'gene_name': [], 'transcript_id': [],
            'position': [], 'species': []
        }

        ids = df['target_id'].dropna().unique().tolist()
        n = 500
        for i in range(0, len(ids), n):
            ids_chunk = ids[i:i+n]
            ext = "/lookup/id"
            ens_lookup = fetch_multiple_endpoints(server, ext, "application/json", ids_chunk)
            for key, value in ens_lookup.items():
                gene_id = key
                gene_name = value.get('display_name', np.nan)
                transcript_id = value.get('canonical_transcript', np.nan)
                species = value.get('species', np.nan)
                position = f"{value['assembly_name']}:{value['seq_region_name']}:{value['start']}:{value['end']}:{value['strand']}"
                isg_dict['gene_id'].append(gene_id)
                isg_dict['gene_name'].append(gene_name)
                isg_dict['transcript_id'].append(transcript_id)
                isg_dict['position'].append(position)
                isg_dict['species'].append(species)

        isg_df = pd.DataFrame(isg_dict)

        # ---------------- Sequences ----------------
        expand_5prime = 1500
        expand_3prime = 1500
        ids = isg_df['gene_id'].tolist()
        n = 50

        print("[INFO] Fetching genomic sequences...")
        for i in range(0, len(ids), n):
            ids_chunk = ids[i:i+n]
            ext = f"/sequence/id?expand_5prime={expand_5prime}&expand_3prime={expand_3prime}"
            ens_upstream = fetch_multiple_endpoints(server, ext, "application/json", ids_chunk)
            for record in ens_upstream:
                gene_id = record['id']
                upstream_seq = record['seq'][:expand_5prime]
                downstream_seq = record['seq'][-expand_3prime:]
                gene_seq = record['seq'][expand_5prime:-expand_3prime]
                isg_df.loc[isg_df['gene_id'] == gene_id, f'upstream_seq_{expand_5prime}bp'] = upstream_seq
                isg_df.loc[isg_df['gene_id'] == gene_id, 'gene_seq'] = gene_seq
                isg_df.loc[isg_df['gene_id'] == gene_id, f'downstream_seq_{expand_3prime}bp'] = downstream_seq

        # CDS
        print("[INFO] Fetching CDS sequences...")
        ids = [i.split('.')[0] for i in isg_df['transcript_id'].dropna().tolist()]
        for i in range(0, len(ids), n):
            ids_chunk = ids[i:i+n]
            ext = "/sequence/id?type=cds"
            ens_cds = fetch_multiple_endpoints(server, ext, "application/json", ids_chunk)
            for record in ens_cds:
                transcript_id = f"{record['id']}.{record['version']}"
                cds = record['seq']
                isg_df.loc[isg_df['transcript_id'] == transcript_id, 'cds'] = cds

        # cDNA
        print("[INFO] Fetching cDNA sequences...")
        ids = [i.split('.')[0] for i in isg_df['transcript_id'].dropna().tolist()]
        for i in range(0, len(ids), n):
            ids_chunk = ids[i:i+n]
            ext = "/sequence/id?type=cdna"
            ens_cdna = fetch_multiple_endpoints(server, ext, "application/json", ids_chunk)
            for record in ens_cdna:
                transcript_id = f"{record['id']}.{record['version']}"
                cdna_seq = record['seq']
                isg_df.loc[isg_df['transcript_id'] == transcript_id, 'cdna_seq'] = cdna_seq

        # Protein
        print("[INFO] Fetching protein sequences...")
        ids = [i.split('.')[0] for i in isg_df['transcript_id'].dropna().tolist()]
        for i in range(0, len(ids), n):
            ids_chunk = ids[i:i+n]
            ext = "/sequence/id?type=protein"
            ens_protein = fetch_multiple_endpoints(server, ext, "application/json", ids_chunk)
            for record in ens_protein:
                transcript_id = record['query']
                protein_id = f"{record['id']}.{record['version']}"
                protein_seq = record['seq']
                isg_df.loc[isg_df['transcript_id'].str.contains(transcript_id), 'protein_id'] = protein_id
                isg_df.loc[isg_df['transcript_id'].str.contains(transcript_id), 'protein_seq'] = protein_seq

        # ---------------- Save final extended file ----------------
        
        isg_df_merge = pd.merge(
            isg_df,
            df,
            left_on='gene_id',
            right_on='target_id',
            how='inner'
        )
        isg_df_merge.drop(columns=['target_id'], inplace=True)
        
        isg_df_merge = isg_df_merge[
            ['source_id', 'source_name', 'gene_id', 'gene_name', 'transcript_id', 'protein_id',
            'species', 'position', 'perc_id', 'perc_pos', 'type', 'method_link_type',
            f'upstream_seq_{expand_5prime}bp', 'gene_seq', f'downstream_seq_{expand_3prime}bp',
            'cdna_seq', 'cds', 'protein_seq']
        ]

        extended_out = "_".join(source_names) + "_ensembl_master_table.csv"
        isg_df_merge.to_csv(os.path.join(folder_out, extended_out), index=False)
        print(f"[INFO] Extended data written to {os.path.join(folder_out, extended_out)}")

        # ---------------- Save final FASTA file ----------------

        isg_df_merge = isg_df_merge[~isg_df_merge['protein_seq'].isna()] # protein sequences can be missing (remove them to avoid an error)
        isg_df_merge['protein_seq_len'] = isg_df_merge['protein_seq'].apply(lambda x: len(x))
        isg_df_merge['protein_X_len'] = isg_df_merge['protein_seq'].apply(lambda x: x.count('X'))
        
        isg_df_merge_sub = isg_df_merge[
            ((isg_df_merge['method_link_type'] == 'ENSEMBL_ORTHOLOGUES') | 
                (
                    (isg_df_merge['method_link_type'] == 'ENSEMBL_PARALOGUES') &
                    (isg_df_merge['gene_name'].isin(source_names))
                )
            ) &
            (isg_df_merge['protein_X_len'] <= 3)
        ].drop_duplicates(subset='protein_seq')


        fasta_file_provisional = os.path.join(folder_out, "_".join(source_names) + "_ensembl.fasta")

        if os.path.exists(fasta_file_provisional):
            os.remove(fasta_file_provisional)

        fasta_seq_number = 0

        for index, row in isg_df_merge_sub.iterrows():

            source_name = row['source_name']
            protein_id = row['protein_id']
            protein_seq = row['protein_seq']

            protein_seq_len = row['protein_seq_len']
            protein_avg_len = isg_df_merge_sub[isg_df_merge_sub['source_name'] == source_name]['protein_seq_len'].mean()

            if protein_seq_len < (protein_avg_len * 0.5):
                continue
            
            with open(fasta_file_provisional, 'a') as f:
                f.write(f">{protein_id}\n{protein_seq}\n")

            fasta_seq_number += 1

        fasta_file_actual = os.path.join(folder_out, str(fasta_seq_number) + "_" + "_".join(source_names) + "_ensembl.fasta")

        if os.path.exists(fasta_file_actual):
            os.remove(fasta_file_actual)

        os.rename(fasta_file_provisional, fasta_file_actual)

        print(f"[INFO] Filtered data written to {fasta_file_actual}")

    
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_folder>")
        sys.exit(1)
    main(sys.argv[1])