#!/bin/bash

#This is a modified version of entrypoint which for picasso HPC for personell UMA and IHSM

# Check for required input directories and files
if [ ! -d "/workflow/data" ]; then
    echo "Error: /workflow/data directory not found. Please mount it using the -v option."
    exit 1
fi

if [ ! -f "/workflow/config.json" ]; then
    echo "Error: config.json file not found. Please mount it using the -v option."
    exit 1
fi

# Check the command option
if [ "$1" == "download" ]; then
    # Step 1: Download dependencies only
    eval "$(conda shell.bash hook)"
    echo "Activating bacLIFE_environment..."
    conda activate bacLIFE_environment

    echo "Downloading general databases..."
    python /workflow/src/download_dependencies.py

    echo "Activating antismash_bacLIFE..."
    conda activate antismash_bacLIFE

    echo "Downloading AntiSMASH databases..."
    download-antismash-databases

elif [ "$1" == "run-snakemake" ]; then
    # Step 2: Run Snakemake only
    eval "$(conda shell.bash hook)"
    echo "Activating bacLIFE_environment..."
    conda activate bacLIFE_environment
    echo "Launching Snakemake..."
    snakemake -j ${2:-1} --use-conda
    echo "Preparing app inputs"
    python src/prepare_app_inputs.py

elif [ "$1" == "rename-genomes" ]; then
    echo "Activating bacLIFE_environment..."
    conda activate bacLIFE_environment
    echo "Rename genomes ..."
    PATH=$PATH:/opt/conda/envs/bacLIFE_environment/lib/R/bin
    Rscript src/rename_genomes.R data/ names_equivalence.txt
	
else
    echo "Invalid option. Use 'download' to download databases, 'rename-genomes' or 'run-snakemake' to run Snakemake."
    exit 1
fi
