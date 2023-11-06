##Script to download databases and other dependencies
import os


#Create folders
if not os.path.exists('intermediate_files'):
   os.makedirs('intermediate_files')
if not os.path.exists('intermediate_files/mapper_data'):
   os.makedirs('intermediate_files/mapper_data')
if not os.path.exists('intermediate_files/PFAM'):
   os.makedirs('intermediate_files/PFAM')
if not os.path.exists('intermediate_files/DBCAN'):
   os.makedirs('intermediate_files/DBCAN')   
if not os.path.exists('intermediate_files/BiG-SCAPE'):
   os.makedirs('intermediate_files/BiG-SCAPE')     


#Download PFAM

os.system('wget -P ./intermediate_files/PFAM/ ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.gz')
os.system('gunzip intermediate_files/PFAM/Pfam-A.hmm.gz' )
os.system('hmmpress intermediate_files/PFAM/Pfam-A.hmm')

#Download EGGNOG

os.system('download_eggnog_data.py -y --data_dir intermediate_files/mapper_data')


#Download DBCAN

os.system('wget -P ./intermediate_files/DBCAN/ http://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V9.txt')

#BiG-SCAPE

os.system('git clone https://github.com/medema-group/BiG-SCAPE.git')
os.system('mv BiG-SCAPE intermediate_files/')




##Phylophlan

import json
with open('config.json', 'r') as config_file:
    config_data = json.load(config_file)

phylo_database_value = config_data["phylo_database"]

cmd = "wget http://cmprod1.cibio.unitn.it/databases/PhyloPhlAn/%s.tar" % phylo_database_value
os.system(cmd)
cmd = "wget http://cmprod1.cibio.unitn.it/databases/PhyloPhlAn/%s.md5" % phylo_database_value
os.system(cmd)


os.system("mkdir -p src/phylophlan_db/")

cmd = "mv %s.tar src/phylophlan_db/%s.tar" % (phylo_database_value, phylo_database_value)
os.system(cmd)


cmd = "mv %s.md5 src/phylophlan_db/%s.md5" % (phylo_database_value, phylo_database_value)
os.system(cmd)

cmd = "phylophlan_write_default_configs.sh"
os.system(cmd)

cmd = "mv super* src/"







