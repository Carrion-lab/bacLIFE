import os
import argparse

def download_genome_metadata():
  os.system('python src/download_NCBI_metadata.py')

def rename_genomes():
  os.system('Rscript src/obtain_names.R')
  
def mash_redundancy():
  os.system('mkdir -p mash')
  os.system('mkdir -p mash/results')
  os.system('python src/mash_redundance.py -i genomes/')
  os.system("awk -F'\t' '$3<0.01' mash/combined_results > mash/combined_results_filtered")
 
def clean_mash_output():
  os.system('mkdir -p genomes_renamed')
  os.system('Rscript src/clean_mash_output.R ')


def rename_genomes2():
  os.system('Rscript src/rename_genomes.R')

def join_metadata():
  os.system('Rscript src/join_metadata.R')

if __name__=='__main__':
    download_genome_metadata()
    rename_genomes()
    mash_redundancy()
    clean_mash_output()
    rename_genomes2()
    join_metadata()
