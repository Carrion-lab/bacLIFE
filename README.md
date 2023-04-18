![workflow](https://user-images.githubusercontent.com/69348873/231155358-7fbebb3c-f6f6-406a-989b-9d273b83aa1e.png)


# Download MicroLife

```
 git clone https://github.com/Carrion-lab/microLife.git microLife
```


# Download genomes and reduce redundancy
### - Install conda environment for download using the file ENVS/MicroLife_download.yml and activate it

```
mamba env create -f ENVS/MicroLife_download.yml
conda activate MicroLife_download
```

### - Download the accession names 
Get a list of accesion ids we are interested (i.e All Mycobacterium RefSeq genomes) by extracting the SQL query from NCBI website search and copying it in the following command

![screenshot](https://user-images.githubusercontent.com/69348873/231155408-8ccbc10e-ce6f-4e24-bf2a-81c6e8f10aad.png)

```
esearch -db assembly -query '("Mycobacterium"[Organism] OR Mycobacterium[All Fields]) AND ("latest refseq"[filter] AND all[filter] NOT anomalous[filter])' | esummary | xtract -pattern         DocumentSummary -element AssemblyAccession > NCBI_accs.txt
```
### - Download the fasta files

```
mkdir genomes/
cd genomes/
bit-dl-ncbi-assemblies -w ../NCBI_accs.txt -f fasta -j 10
```

### - Uncompress the fasta files

```
gzip -d ./*.gz
```

### - Run main script to download metadata and remove redundancy at >99% Average Nucleotide Identity (ANI) with mash

```
cd ../
python main_script_download.py
```

### - Outputs:
- "genomes_renamed/": Folder where the non-redundant genomes with the filenames ready for MicroLife are stored
- "METADATA_MERGED.txt": File with the following columns 'scientific_name', 'NCBIid', 'cluster_membership', 'representative' and 'MicroLife_name'

# Run MicroLife

## Clustering module
### - Install conda environments for MicroLife

```
mamba env create -f ENVS/MicroLife_environment.yml
mamba env create -f ENVS/antismash.yml
mamba env create -f ENVS/bigscape.yml
conda activate MicroLife_environment
```

### - Prepare input
The input must be fasta assemblies stored in the data/ directory following format "Genus_species_strain_O.fna"*. If genomes were downloaded as described above in the section 'Download genomes and reduce redundancy' the folder download/genomes_renamed/ contains the genomes in the correct format and they only have to be moved to the data/ folder.

*In the genome format name user must avoid the presence of any symbol that is not a alphabet letter, number , dash or dot. This applies specially to the strain name which also should be of less than 10 characters long. Users may use the script src/rename_genomes.R to change the input genome strain names into a MicroLife friendly format in the following way:

```
Rscript src/rename_genomes.R data/ names_equivalence.txt
```
This script changes the name of the input genome files by replacing the strain name into a barcode to avoid MicroLife crashes \
I.e Pseudomonas_fluorescens_FDAARGOS-1088.fna --> Pseudomonas_fluorescens_X0001.fna \
The file "names_equivalence.txt" contains the correspondent name matching.

### - Run MicroLife clustering module
MicroLife is written using the snakemake workflow manager and it can be executed using the following command from the main directory

```
snakemake -j 24
```

### - Create your own metadata

To be able to run the lifestyle prediction and the shiny app you have to make a "mapping_file.txt" for your genomes.
You can find an example of the mapping file in the main folder. It consists of a two column tab-separated file. First column specifies the fasta file name (as established in the `data/` folder) and second column specifies the Lifestyle of that genome fasta file given by the user. Genomes which not known lifestyle must be annotated as "Unknown"

### - MicroLife outputs

MEGAMATRIX.txt: This file is the main output of MicroLife and consists in a matrix where each row represent one gene cluster and each column represent the presence of these clusters in the different genomes and its annotation in different databases

# microLife Lifestyle prediction module
User can test the predictability of their metadata with the machine learning model random forest using the script 'src/classifier.R'

```
Rscript src/classifier.R mapping_file.txt MEGAMATRIX_renamed.txt
```
ROC plots showing the model evaluation for the different classes present in the metadata are generated and stored in classifier/

If good accuracy is accomplished, user can use the model for predictions of genomes labeled as 'Unknown' in the mapping file 
Metadata augmented is stored in the new file "mapping_file_augmented.txt" 

# MicroLife App

### - Prepare input 
In order to initiate the shiny app, the following snakemake output files need to be droped in the 'Shiny_app/input/' directory.

- MEGAMATRIX_renamed.txt
- mapping_file.txt or mapping_file_augmented.txt
- intermediate_files/BiG-SCAPE/big_scape_binary_table_renamed.txt
- intermediate_files/BiG-SCAPE/annotation.txt
- intermediate_files/combined_proteins/combined_proteins.fasta

After input preparation users can open the script "app.R" in Rstudio and click in the upper right 'run app' bottom to initiate the MicroLife app. 
