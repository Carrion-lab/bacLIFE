
import json
import os
import argparse
import glob

file_id = []
for name in glob.glob('genomes/*.fa'):
    file_id.append(name)

print(file_id)

ncbi_ids = []
sci_names = []
str_names = []

os.system('mkdir JSON')
for x in file_id:
  json_name = x.strip()
  json_name = json_name.replace('.fa', '')
  json_name = json_name.replace('genomes/', '')
  json_file = str('JSON/') + str(json_name)
  print(json_file)
  os.system('datasets summary genome accession %s > %s' % (json_name, json_file))

  # Opening JSON file
  f = open(json_file)
    
  # returns JSON object as 
  # a dictionary
  data = json.load(f)
  
  key1= 'reports'
  if key1 in data.keys():
    to_print = data['reports']
    to_print= dict(to_print[0])
    #print(to_print.keys())

    to_print= to_print['organism']
    #print(to_print.keys())
    to_print_name= to_print['organism_name']
    #print(to_print.keys())
    sci_names.append(to_print_name)
    
    
    key = 'infraspecific_names'
    if key in to_print.keys():
      to_print = to_print['infraspecific_names']
      if 'strain' in to_print.keys():
        to_print_strain= to_print['strain']
      else:
        if 'isolate' in to_print.keys():
          to_print_strain= to_print['isolate']
        else:
          to_print_strain= 'unknown'
    else:
      to_print_strain= 'unknown'
    
    str_names.append(to_print_strain)
    ncbi_ids.append(json_name)
  else:
    continue
###Write file


output = open("scientific_names.txt", "w")

for i in range(len(sci_names)):
  output.write(ncbi_ids[i])
  output.write('\t')
  output.write(sci_names[i])
  output.write('\t')
  output.write(str_names[i])
  output.write('\n')


