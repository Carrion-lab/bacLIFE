# GBF to bacLIFE LAG/DEG mapper

This utility links gene identifiers from a GenBank/GBFF file to bacLIFE cluster annotations and then merges the result with a bacLIFE LAGs table and a DEGs table.

Run this script inside the bacLIFE_environment conda environment, since it depends on tools and Python packages available there.

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

The DEGs table is expected to contain the gene identifier used for the join in a column called `ID`.

That identifier must match the qualifier selected with `--gbf-tag`, so you should inspect both the DEG table and the GBFF file to determine which CDS qualifier is the correct one.

For example, if the `ID` column in the DEGs file contains `locus_tag` values, run the script with:

```bash
--gbf-tag locus_tag
