# GBF to bacLIFE LAG/DEG mapper

This utility links gene identifiers from a GenBank/GBFF file to bacLIFE cluster annotations and then merges the result with a bacLIFE LAGs table and a DEGs table.

## What it does

The script:

1. extracts protein sequences from the input `--gbf` file using the qualifier selected with `--gbf-tag`
2. aligns those proteins against `--combined-proteins` using DIAMOND
3. keeps the best hit per query protein above the identity threshold
4. maps those hits to bacLIFE clusters using the `descriptions` column in `--megamatrix`
5. merges the mapping with:
   - a bacLIFE LAGs Excel table from the app (`--baclife-LAGs`)
   - a DEGs table (`--DEGs`)

The final merged output only keeps the overlap between:
- mapped bacLIFE clusters present in the LAGs table
- GBF gene identifiers present in the DEGs table

## Inputs

- `--combined-proteins`: protein FASTA produced by the bacLIFE workflow/app
- `--megamatrix`: megamatrix file produced by the bacLIFE workflow/app
- `--baclife-LAGs`: Excel file with LAGs generated from the app when comparing groups
- `--DEGs`: CSV file with differentially expressed genes with the column with gene ids called `ID`
- `--gbf`: GenBank/GBFF annotation file
- `--gbf-tag`: qualifier from the GBFF CDS entries to use as gene ID

## Important note about `--gbf-tag` and the DEGs table

The first column in the DEGs table is expected to contain the gene identifier used for the join.

That identifier must match the qualifier selected with `--gbf-tag`.

For example, if the first column in the DEGs file contains `locus_tag` values, run the script with:

```bash
--gbf-tag locus_tag
```

If the DEGs file instead uses `old_locus_tag` or `gene`, then `--gbf-tag` should be set to that same qualifier.

## Usage

```bash
python gbf_to_baclife_lag_deg_mapper.py \
  --baclife-LAGs plant_pathogen_LAGs.xlsx \
  --megamatrix MEGAMATRIX_renamed.txt \
  --combined-proteins combined_proteins.fasta \
  --DEGs DEGs.csv \
  --gbf genome.gbff \
  --gbf-tag locus_tag \
  --output results \
  --threads 4
```

## Outputs

The script writes:

- `baclife_gene_cluster_mapping.tsv`: mapping between bacLIFE clusters, bacLIFE genes, and GBF gene names
- `baclife_LAGs_DEGs_merged.tsv`: merged table with the overlap between mapping, LAGs, and DEGs

## Suggested location in the repository

```text
src/gbf_to_baclife_lag_deg_mapper/
├── README.md
└── gbf_to_baclife_lag_deg_mapper.py
```
