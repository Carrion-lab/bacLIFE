# bacLIFE in pixi — Setup, Usage & Troubleshooting Guide

**Repository:** https://github.com/bharat1912/taxonomy_bundle  
**Pipeline:** [bacLIFE](https://github.com/ClavelLab/bacLIFE) — Lifestyle-Associated Gene identification  
**Environment:** `env-baclife` in pixi multilayer environment  
**Status:** Successfully validated on thermohalophile comparative genomics project (23 genomes, March 2026)

---

## Table of Contents

1. [Overview](#overview)
2. [Environment Setup](#environment-setup)
3. [Input Requirements](#input-requirements)
4. [Genome Download Strategies](#genome-download-strategies)
5. [Genome QC with CheckM2](#genome-qc-with-checkm2)
6. [Data Organisation](#data-organisation)
7. [Running the Pipeline](#running-the-pipeline)
8. [Known Issues and Fixes](#known-issues-and-fixes)
9. [Post-Pipeline Analysis](#post-pipeline-analysis)
10. [Test Data Validation](#test-data-validation)
11. [Shiny App](#shiny-app)

---

## Overview

bacLIFE identifies Lifestyle-Associated Genes (LAGs) by comparing gene presence/absence across genomes grouped by lifestyle (e.g. Thermohalophile, Thermophile, Mesophile). The pipeline runs:

**Prokka** (annotation) → **MMseqs2 + MCL** (clustering) → **Presence/absence matrix** → **Random Forest** (LAG identification) → **antiSMASH + BiG-SCAPE** (BGC analysis) → **PhyloPhlAn** (phylogeny) → **Shiny app** (visualisation)

---

## Environment Setup

bacLIFE runs in the `env-baclife` pixi environment. Key tools included:

- `prokka 1.14.6` — genome annotation
- `mmseqs2` — protein clustering
- `mcl` — Markov clustering
- `antismash 8.0.4` — biosynthetic gene cluster prediction
- `snakemake` — workflow management
- `r-randomforest` — LAG identification
- `ncbi-datasets-cli` — genome downloads

**Available pixi tasks:**

```bash
pixi run -e env-baclife baclife-run        # Run the full pipeline
pixi run -e env-baclife baclife-link-assets # Link vault genomes into data/input/
pixi run -e env-baclife baclife-help        # Show available tasks
```

---

## Input Requirements

### 1. Genome FASTA files — naming convention (CRITICAL)

The Snakefile discovers genomes by globbing `bacLIFE/data/*.fna` using the pattern:

```
{Genus}_{species}_{strain}_{replicon}.fna
```

**Exactly 4 underscore-delimited tokens are required.** The replicon token is conventionally `O` (for chromosome).

Examples of correct naming:
```
Halothermothrix_orenii_DSM18212_O.fna
Bacillus_subtilis_168_O.fna
Petrotoga_mobilis_SJ95_O.fna
```

Files must be placed directly in `bacLIFE/data/` — **not** in `bacLIFE/data/input/` (the `baclife-link-assets` task creates symlinks there but the Snakefile does not read from that subdirectory).

### 2. mapping_file.txt

Tab-separated, two columns, placed at `bacLIFE/mapping_file.txt`:

```
Sample	Lifestyle
Halothermothrix_orenii	Thermohalophile
Bacillus_subtilis	Mesophile
```

The `Sample` column must match your full organism name (as used in `names_equivalence.txt`), **not** the 4-part bacLIFE internal name.

**Important:** The pipeline's `rename_MEGAMATRIX` step overwrites `mapping_file.txt` with `Lifestyle = Unknown` for all samples. Back up your annotated version and restore it after the run:

```bash
# After baclife-run completes:
cp /path/to/your/mapping_file_annotated.txt bacLIFE/mapping_file.txt
```

### 3. names_equivalence.txt

Maps bacLIFE's internal 4-part names to your sample names. Place at `bacLIFE/names_equivalence.txt`. The internal name is derived from the filename with `_O` stripped:

```
bacLIFE_name	Full_name
Halothermothrix_orenii_DSM18212	Halothermothrix_orenii
Bacillus_subtilis_168	Bacillus_subtilis
```

**How to get the exact internal names** — run the pipeline first, then extract from MEGAMATRIX.txt:

```bash
head -1 bacLIFE/MEGAMATRIX.txt | tr ' ' '\n' | sed 's/"//g' | \
  awk '/^clusters$/{found=1; next} /^completeness$/{exit} found{print}' | sort
```

---

## Genome Download Strategies

Three strategies are required depending on accession type. All downloads should go to `$EXTERNAL_VAULT/ncbi_genomes/`. Add your NCBI API key to `~/.bashrc` to avoid rate limits.

### Strategy 1 — GCF/GCA assemblies (ncbi-datasets-cli)

```bash
pixi run -e env-baclife datasets download genome accession GCF_000371785.1 \
  --include genome --no-progressbar --filename out.zip
unzip out.zip -d ncbi_download
find ncbi_download -name "*.fna" -exec cp {} $EXTERNAL_VAULT/ncbi_genomes/MyGenome.fna \;
```

Note: Downloaded filenames use verbose NCBI names — rename manually to your `Sample_ID.fna` convention.

### Strategy 2 — NZ_ WGS master accessions

NZ_ accessions are WGS master records with no sequence themselves. Use esearch + wget (not entrez-direct, which fails with OpenSSL 3.6.1 inside pixi):

```bash
# Step 1: get WebEnv/QueryKey
wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?\
db=nuccore&term=NZ_MCIB01000000[ACCN]&usehistory=y&api_key=${NCBI_API_KEY}" \
  -O /tmp/search.xml

# Step 2: fetch all contigs (adjust retmax to number of contigs)
wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?\
db=nuccore&query_key=1&WebEnv=PASTE_WEBENV_HERE&rettype=fasta&retmax=41&api_key=${NCBI_API_KEY}" \
  -O $EXTERNAL_VAULT/ncbi_genomes/MyGenome.fna
```

### Strategy 3 — Legacy accessions (CP/FP/AE/AL/BX/U prefixes)

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?\
db=nuccore&id=CP001098&rettype=fasta&retmode=text&api_key=${NCBI_API_KEY}" \
  -O $EXTERNAL_VAULT/ncbi_genomes/MyGenome.fna
sleep 0.4  # respect NCBI rate limit
```

**Important:** Run wget from a plain terminal, not inside `pixi run` — OpenSSL 3.6.1 in pixi environments is incompatible with NCBI's TLS and produces `error:0A00006E: bad extension`.

---

## Genome QC with CheckM2

```bash
# Register database (path to .dmnd FILE, not directory)
pixi run -e env-checkm2 checkm2 database \
  --setdblocation $EXTERNAL_VAULT/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd

# Run QC
pixi run -e env-checkm2 checkm2 predict --threads 16 \
  --input bacLIFE/data/input \
  --output-directory bacLIFE/quality_audit \
  --extension fna --force
```

**Thresholds:** Completeness ≥ 90%, Contamination ≤ 5%.

**Note for archaeal genomes:** CheckM2 is trained on bacteria. Archaeal completeness estimates may be lower (~93%) even for complete genomes — treat with caution rather than excluding.

---

## Data Organisation

```
bacLIFE/
├── data/                          ← Snakefile globs HERE for *.fna
│   ├── Halothermothrix_orenii_DSM18212_O.fna  ← 4-part symlinks
│   ├── Bacillus_subtilis_168_O.fna
│   └── ...
├── data/input/                    ← Created by baclife-link-assets (NOT read by Snakefile)
├── pseudomonas_test_data/         ← Test data stored OUTSIDE data/ to avoid mixing
├── mapping_file.txt               ← Sample → Lifestyle
├── names_equivalence.txt          ← bacLIFE internal name → Sample_ID
└── config.json                    ← Pipeline parameters
```

**Critical rule:** Never put test data and project data in `data/` simultaneously. The Snakefile picks up ALL `*.fna` files in `data/` and its subdirectories. Store test data in a separate directory outside `data/`.

### Project switching

```bash
# Switch to test data
rm bacLIFE/data/*.fna
ln -sfn bacLIFE/pseudomonas_test_data/*.fna bacLIFE/data/

# Switch to project data
rm bacLIFE/data/*.fna
ln -sfn $EXTERNAL_VAULT/ncbi_genomes/baclife_ready/*.fna bacLIFE/data/
```

---

## Running the Pipeline

Always run in tmux to survive terminal disconnection:

```bash
tmux new -s baclife
cd ~/software/taxonomy_bundle
pixi run -e env-baclife baclife-run
# Ctrl+B then D to detach
# tmux attach -t baclife to reattach
```

**Monitor Prokka progress** (in a separate terminal):

```bash
# Count completed annotations
grep -l "Annotation finished" \
  bacLIFE/intermediate_files/annot/*/*.log 2>/dev/null | wc -l

# See which genomes are done
grep -l "Annotation finished" \
  bacLIFE/intermediate_files/annot/*/*.log 2>/dev/null
```

**Resume after failure:** Snakemake automatically resumes from the last completed step. If the pipeline exits due to a failed job, fix the issue and rerun:

```bash
cd ~/software/taxonomy_bundle
pixi run -e env-baclife baclife-run
```

---

## Known Issues and Fixes

### Issue 1 — MMseqs2 stale database files

**Symptom:**
```
intermediate_files/clustering/mmseqDB_clu.dbtype exists already!
CalledProcessError: mmseqs cluster returned non-zero exit status 1
```

**Cause:** MMseqs2 refuses to overwrite its own database files from a previous run.

**Fix:**
```bash
rm -rf bacLIFE/intermediate_files/clustering/mmseqDB*
pixi run -e env-baclife baclife-run
```

---

### Issue 2 — BiG-SCAPE Network_Annotations_Full.tsv path mismatch

**Symptom:**
```
Error in rule extract_binary_table_GCF
Missing output files: .../hybrids_glocal/Network_Annotations_Full.tsv
```

**Cause:** The Snakefile's `bigscape_exe` rule runs BiG-SCAPE and then attempts to reorganise its output directory with:
```bash
rm -r .../hybrids_glocal
mv .../*hybrids_glocal .../hybrids_glocal
```
If BiG-SCAPE has already run and `hybrids_glocal` exists from a previous run, Snakemake skips `bigscape_exe` but the `mv` never executes, leaving `Network_Annotations_Full.tsv` in the wrong location.

**Fix:** Delete the BiG-SCAPE output entirely and rerun:
```bash
rm -rf bacLIFE/intermediate_files/BiG-SCAPE/bigscape_output/
pixi run -e env-baclife baclife-run
```

**Permanent fix (Snakefile patch):** See [GitHub fix section](#snakefile-fix-for-bigscape-path-issue) below.

---

### Issue 3 — FASTA file leading newline causing Prokka failure

**Symptom:**
```
MSG: The sequence does not appear to be FASTA format (lacks a descriptor line '>')
STACK Bio::SeqIO::fasta::next_seq
```

**Cause:** Some NCBI eutils downloads prepend a blank line (`\n`) before the first `>` header. BioPerl's FASTA parser requires the file to begin with `>` as the very first byte.

**Diagnosis:**
```bash
xxd your_genome.fna | head -2
# If first byte is 0a (not 3e), you have this problem
```

**Real example — Thermohalobacter berrensis (NZ_MCIB01000000):**
```
00000000: 0a3e 4e5a 5f4d 4349 42...   ← 0a = newline BEFORE >
```

**Fix — strip blank lines and simplify headers:**
```bash
sed '/^$/d' your_genome.fna | \
  awk '/^>/{split($1,a,">"); print ">"a[2]} !/^>/{print}' \
  > your_genome_clean.fna

# Verify first byte is now 3e (>)
xxd your_genome_clean.fna | head -2

# Update symlink to point to clean file
ln -sfn /full/path/your_genome_clean.fna \
  bacLIFE/data/YourGenome_strain_O.fna
```

---

### Issue 4 — names_equivalence.txt mismatch causing rename_MEGAMATRIX failure

**Symptom:**
```
Error in colnames(matrix)[2:n_samples] <- M$Full_name :
  replacement has length zero
```

**Cause:** The `bacLIFE_name` column in `names_equivalence.txt` does not match the internal names bacLIFE generated in `MEGAMATRIX.txt`. This happens when the file is generated with incorrect suffixes (e.g. `_O_O.fna` instead of the correct internal name).

**Fix:** Extract the exact internal names from MEGAMATRIX and rebuild the file:

```bash
# Get exact internal names
head -1 bacLIFE/MEGAMATRIX.txt | tr ' ' '\n' | sed 's/"//g' | \
  awk '/^clusters$/{found=1; next} /^completeness$/{exit} found{print}' | sort

# Then rebuild names_equivalence.txt manually with these exact names
```

---

### Issue 5 — Snakemake lock from killed sessions

**Symptom:**
```
LockException: Directory cannot be locked
```

**Fix:**
```bash
cd bacLIFE
pixi run -e env-baclife snakemake --unlock --cores 16 -s Snakefile --configfile config.json
```

Note: Run the unlock command from inside `bacLIFE/` where `config.json` is located.

---

### Issue 6 — Zombie cmscan/prokka processes from killed sessions

**Symptom:** After killing tmux sessions, old Prokka/cmscan processes continue running in background consuming cores.

**Fix:**
```bash
pkill -9 -f cmscan
pkill -9 -f prokka
pkill -9 -f snakemake
sleep 3
pgrep -la cmscan || echo "clear"
pgrep -la prokka || echo "clear"
```

---

## Snakefile Fix for BiG-SCAPE Path Issue

The root cause of Issue 2 is in the `bigscape_exe` rule shell command (Snakefile line ~306). The current command:

```bash
python ./databases/BiG-SCAPE/bigscape.py ... ; \
rm -r intermediate_files/BiG-SCAPE/bigscape_output/network_files/hybrids_glocal; \
mv intermediate_files/BiG-SCAPE/bigscape_output/network_files/*hybrids_glocal \
   intermediate_files/BiG-SCAPE/bigscape_output/network_files/hybrids_glocal
```

**The fix** — add a pre-run cleanup step so stale files never block the mv:

```bash
# Replace the shell command in the bigscape_exe rule with:
shell:
    """
    rm -rf intermediate_files/BiG-SCAPE/bigscape_output/network_files/hybrids_glocal
    python ./databases/BiG-SCAPE/bigscape.py -i {params.indir} -o {params.outdir} \
      --pfam_dir databases/PFAM/ --mode glocal --mibig \
      --cutoffs 0.3 0.7 --include_singletons --cores {params.threads} --mix
    mv intermediate_files/BiG-SCAPE/bigscape_output/network_files/*hybrids_glocal \
       intermediate_files/BiG-SCAPE/bigscape_output/network_files/hybrids_glocal
    """
```

This fix can be committed to the repository as a one-line change to the Snakefile. To do so:

```bash
cd ~/software/taxonomy_bundle
git diff bacLIFE/Snakefile   # review the change
git add bacLIFE/Snakefile
git commit -m "fix: pre-clean hybrids_glocal before BiG-SCAPE mv to prevent path mismatch on rerun"
git push
```

---

## Post-Pipeline Analysis

After the pipeline completes, restore your annotated mapping file:

```bash
cp your_mapping_file_annotated.txt bacLIFE/mapping_file.txt
```

Key output files:

| File | Description |
|---|---|
| `MEGAMATRIX_renamed.txt` | Gene presence/absence matrix with sample names |
| `intermediate_files/BiG-SCAPE/big_scape_binary_table_renamed.txt` | BGC presence/absence matrix |
| `intermediate_files/phylophlan/output_phylophlan/RAxML_bestTree.input_refined.tre` | Phylogenetic tree |
| `mapping_file.txt` | Sample → Lifestyle (restore your annotated version) |

### LAG extraction (R script)

Use `extract_TH_LAGs_v2.R` (available in project directory) to extract lifestyle-associated gene clusters from MEGAMATRIX. The script uses fast `rowMeans()` matrix operations and produces:

- `TH_LAGs_enriched.csv` — all enriched clusters
- `TH_LAGs_unique.csv` — clusters unique to target lifestyle
- `TH_LAGs_temperature_signal.csv` — shared with thermophiles
- `TH_LAGs_salt_signal.csv` — shared with halophiles
- `all_clusters_ranked.csv` — full ranked list

```bash
cd bacLIFE
pixi run -e env-baclife Rscript extract_TH_LAGs_v2.R
```

---

## Test Data Validation

The bacLIFE repository includes 21 *Pseudomonas* test genomes for pipeline validation. These are stored in `bacLIFE/pseudomonas_test_data/`.

**To run the validation:**

```bash
# 1. Clear current data symlinks
rm bacLIFE/data/*.fna

# 2. Link test data
ln -sfn ~/software/taxonomy_bundle/bacLIFE/pseudomonas_test_data/*.fna \
  ~/software/taxonomy_bundle/bacLIFE/data/

# 3. Wipe intermediate files from previous runs
rm -rf bacLIFE/intermediate_files/

# 4. Unlock if needed
cd bacLIFE
pixi run -e env-baclife snakemake --unlock --cores 16 -s Snakefile --configfile config.json

# 5. Run in tmux
tmux new -s baclife_test
cd ~/software/taxonomy_bundle
pixi run -e env-baclife baclife-run
```

**Known issue with test data:** The original test run (March 2026) failed at `extract_binary_table_GCF` due to Issue 2 (BiG-SCAPE path mismatch). Apply the fix in Issue 2 before running.

---

## Shiny App

The Shiny app reads from `bacLIFE/Shiny_app/input/`. Copy pipeline outputs before launching:

```bash
cp bacLIFE/MEGAMATRIX_renamed.txt bacLIFE/Shiny_app/input/
cp bacLIFE/mapping_file.txt bacLIFE/Shiny_app/input/
cp bacLIFE/intermediate_files/phylophlan/output_phylophlan/RAxML_bestTree.input_refined.tre \
   bacLIFE/Shiny_app/input/
cp bacLIFE/intermediate_files/BiG-SCAPE/big_scape_binary_table_renamed.txt \
   bacLIFE/Shiny_app/input/

# Delete cached data to force regeneration from your files
rm -f bacLIFE/Shiny_app/data.Rdata

# Launch
cd bacLIFE
pixi run -e env-baclife Rscript -e "shiny::runApp('Shiny_app')"
```

The app will be available at `http://127.0.0.1:PORT` (port shown in terminal output). Navigation warnings about `bslib::nav_panel()` are cosmetic and do not affect functionality.

---

*Guide compiled from bacLIFE thermohalophile project run, March 2026.*
*Pipeline version: bacLIFE (taxonomy_bundle), pixi env-baclife.*
