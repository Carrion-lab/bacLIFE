import os
import argparse


def run_mash(input_folder):
  fileDir = input_folder
  fileExt = r".fa"
  
  #Prepare mash database
  os.system(str('mash sketch -o mash/reference ' + str(input_folder) + '*.fa'))
  
  lst = [_ for _ in os.listdir(fileDir) if _.endswith(fileExt)]
  
  for i in lst:
    os.system(str('mash dist mash/reference.msh '+ str(input_folder) + i + ' > mash/results/results_' + i))
  
  os.system('cat mash/results/* > mash/combined_results')

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument('-d', '--diamond_input', dest = 'diamond_input', help = 'diamond results', required = True)
    parser.add_argument('-i', '--folder_input', dest = 'folderinput', help = 'folder where fasta files are stored', required = True)
    args = parser.parse_args()
    run_mash(args.folderinput)
  
##Line to filter huge table 
#awk -F"\t" '$3<0.03' combined_results > combined_results_filtered
