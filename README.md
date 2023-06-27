# microLife: an automated genome mining tool for identification of lifestyle associated genes

microLife is a streamlined computational workflow that annotates bacterial genomes and performs large-scale comparative genomics to predict bacterial lifestyles and to pinpoint candidate genes, denominated  **lifestyle-associated genes (LAGs)**, and biosynthetic gene clusters associated with each lifestyle detected. This whole process is divided into different modules:

- **Clustering module**
	Predicts, clusters and annotates the genes of every input genome
- **Lifestyle prediction**
	Employs a machine learning model to forecast bacterial lifestyle or other specified metadata
- **Analitical module (Shiny app)**
	Results from the previous modules are embedded in a user-friendly interface for comprehensive and interactive comparative genomics


![workflow](https://user-images.githubusercontent.com/69348873/231155358-7fbebb3c-f6f6-406a-989b-9d273b83aa1e.png)

## How to start using microLife
### Download microLife and install dependencies

```
 git clone https://github.com/Carrion-lab/microLife.git microLife
 cd microLife/
```
Install the conda environment necessary to download files from NCBI using the file `ENVS/MicroLife_download.yml` and activate it:

```
mamba env create -f ENVS/MicroLife_download.yml
conda activate MicroLife_download
```

### Download genomes and reduce redundancy



#### Download the accession names 
Retrieve a list of accesion IDs you are interested in processing with microLife (e.g All Mycobacterium RefSeq genomes) by extracting the SQL query (image below) from NCBI website search and integrating it in the following command:

```
esearch -db assembly -query '("Mycobacterium"[Organism] OR Mycobacterium[All Fields]) AND ("latest refseq"[filter] AND all[filter] NOT anomalous[filter])' | esummary | xtract -pattern         DocumentSummary -element AssemblyAccession > NCBI_accs.txt
```

![screenshot](https://user-images.githubusercontent.com/69348873/231155408-8ccbc10e-ce6f-4e24-bf2a-81c6e8f10aad.png)


#### Download and uncompress the FASTA files

```
mkdir genomes/
cd genomes/
bit-dl-ncbi-assemblies -w ../NCBI_accs.txt -f fasta -j 10
gzip -d ./*.gz
```

#### Download metadata and remove redundancy at >99% Average Nucleotide Identity (ANI)

```
cd ../
python main_script_download.py
```

#### Making sure everything is ready
Before starting microLife, you should have the following in your main directory (`microLife/download/`):
-   "genomes_renamed/": Directory where the non-redundant genomes with the filenames ready for microLife are stored.
-   "METADATA_MERGED.txt": File with the following columns: 'scientific_name', 'NCBIid', 'cluster_membership', 'representative' and 'MicroLife_name'

(Remember to deactivate the conda environment `MicroLife_download` if you want to jump to the next section (**"Executing microLife"**) right after this one!

## Executing clustering module
### Install dependencies necessary for MicroLife

```
mamba env create -f ENVS/MicroLife_environment.yml
mamba env create -f ENVS/antismash.yml
mamba env create -f ENVS/bigscape.yml
conda activate MicroLife_environment
```

### Prepare input files
The input must be FASTA assemblies stored in the `data/` directory, all with name that should follow the format "*Genus_species_strain_O.fna*". *
If genomes were downloaded as described in the previous section **'Download genomes and reduce redundancy'**, the folder directory `download/genomes_renamed/` already contains the genome files in the correct name format and they only have to be moved to the `data/` folder.
```
mv download/genomes_renamed/* data/
```

*In the genome format name, the user must avoid adding any special characters, that is not a alphabet letter, number , dash or dot. This applies specially to the strain name which also should be of less than 10 characters long. Users may use the script *src/rename_genomes.R* to change the input genome strain names into a MicroLife friendly format in the following way:

```
Rscript src/rename_genomes.R data/ names_equivalence.txt
```
This script changes the name of the input genome files by replacing the strain name into a barcode to avoid microLife crashes \
I.e Pseudomonas_fluorescens_FDAARGOS-1088.fna --> Pseudomonas_fluorescens_X0001.fna \
The file "names_equivalence.txt" contains the correspondent name matching.

### Executing Snakemake
microLife is written using the Snakemake workflow manager and it can be executed using the following command from the main directory

```
snakemake -j 24
```

> (`-j` specifies the number of CPU cores to use for the whole process. If this flag is ommitted, Snakemake will use all the available CPU cores in the machine)

#### Create your own metadata

To be able to run the Lifestyle prediction module and the Shiny app (next modules of the microLife pipeline) you have to make a `mapping_file.txt` for your genomes. You can find an example of the mapping file in the main folder {INSERT PATH}. It consists of a two column tab-separated file where the first column specifies the FASTA file name (as established in the `data/` directory) and the second column specifies the Lifestyle of that genome given by the user(WHEN DID THE USER DO THIS??). Genomes with no lifestyle information must be annotated as "Unknown".

#### microLife output

The main (and very big) output of microLife is a file called `MEGAMATRIX_renamed.txt`, a matrix where each row represents one gene cluster and each column represents the presence of these clusters in the different genomes and their annotation in different databases


## Lifestyle prediction module
The user can test the predictability of their metadata (WHY WOULD THE USER WANT THIS?) with the machine learning model random forest using the script `src/classifier.R` and specifying the location of `mapping_file.txt` and `MEGAMATRIX_renamed.txt` in addition to the column name of the metadata variable that the user wants to predict

```
Rscript classifier_src/classifier.R mapping_file.txt MEGAMATRIX_renamed.txt Lifestyle
```
ROC plots showing the model evaluation for the different classes present in the metadata are generated and stored in `classifier/`

If good accuracy is accomplished, the user can use the model for predictions of genomes labeled as 'Unknown' in the mapping file. The resulting augmented metadata is stored in the new file `mapping_file_augmented.txt`


## microLife App (Analytical module)

In order to initiate the shiny app, the following snakemake output files need to be droped in the 'Shiny_app/input/' directory.


-   `MEGAMATRIX_renamed.txt`
-   `mapping_file.txt` or `mapping_file_augmented.txt`
-   `intermediate_files/BiG-SCAPE/big_scape_binary_table_renamed.txt`
-   `intermediate_files/BiG-SCAPE/annotation.txt`
-   `intermediate_files/combined_proteins/combined_proteins.fasta`
-   `names_equivalence.txt`

After having prepared these input files, users can open the script `app.R` in Rstudio and click in the upper right '**run app**' button to initiate the microLife app and enjoy visualizing all their results!

An example of the app with a demo dataset (genomes present in (`data`) is available at [http://178.128.251.24:3838/microLife_linux]
