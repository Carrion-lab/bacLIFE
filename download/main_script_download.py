import os
import argparse

def download_genome_metadata():
  os.system('python src/download_NCBI_metadata.py')

def rename_genomes():
  os.system('Rscript src/obtain_names.R')
  
def mash_redundancy(threshold=0.01):
  os.system('mkdir -p mash')
  os.system('mkdir -p mash/results')
  os.system('python src/mash_redundance.py -i genomes/')
  os.system(f"awk -F'\t' '$3<{threshold}' mash/combined_results > mash/combined_results_filtered")
 
def clean_mash_output(threshold=0.01):
  os.system('mkdir -p genomes_renamed')
  os.system(f'Rscript src/clean_mash_output.R {threshold}')


def rename_genomes2():
  os.system('Rscript src/rename_genomes.R')

def join_metadata():
  os.system('Rscript src/join_metadata.R')

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Download and process genomes')
    parser.add_argument('--threshold', type=float, default=0.01, 
                        help='Mash distance threshold for redundancy filtering (default: 0.01)')
    args = parser.parse_args()
    
    download_genome_metadata()
    rename_genomes()
    mash_redundancy(args.threshold)
    clean_mash_output(args.threshold)
    rename_genomes2()
    join_metadata()
