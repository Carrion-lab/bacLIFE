FROM continuumio/miniconda3:latest

# Set working directory inside the container
WORKDIR /workflow

# Copy the entire workflow (including Snakefile, scripts, and config)
COPY . /workflow

# Create and activate the main conda environment (bacLIFE_environment)
RUN conda env create -f ENVS/bacLIFE_environment.yml && \
    conda clean -a

# Create and activate the additional two environments used in specific rules
RUN conda env create -f ENVS/antismash.yml && \
    conda env create -f ENVS/bigscape.yml && \
    conda clean -a

# Ensure the environments are activated within the container
SHELL ["conda", "run", "-n", "bacLIFE_environment", "/bin/bash", "-c"]

# Define an entrypoint script that downloads databases and runs Snakemake
COPY entrypoint.sh /workflow/entrypoint.sh
RUN chmod +x /workflow/entrypoint.sh

# Expose a directory for data input/output (to be mounted as a volume)
VOLUME ["/data"]

# Set the entry point to the script
ENTRYPOINT ["/workflow/entrypoint.sh"]