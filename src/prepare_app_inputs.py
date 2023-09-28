##Script to move all the outputs into the input/ app folder
import os

os.system('mv MEGAMATRIX_renamed.txt Shiny_app/input/')

os.system('mv mapping_file.txt Shiny_app/input/')

os.system('mv names_equivalence.txt Shiny_app/input/')

os.system('mv intermediate_files/BiG-SCAPE/big_scape_binary_table_renamed.txt Shiny_app/input/')

os.system('mv intermediate_files/BiG-SCAPE/annotation.txt Shiny_app/input/')

os.system('mv intermediate_files/combined_proteins/combined_proteins.fasta Shiny_app/input/')

os.system('mv intermediate_files/BiG-SCAPE/bigscape_output/BGC_descriptions.txt Shiny_app/input/')
