##Script to download databases and other dependencies
import os


#Create folders
if not os.path.exists('intermediate_files'):
   os.makedirs('intermediate_files')
if not os.path.exists('databases'):
   os.makedirs('databases')
if not os.path.exists('databases/mapper_data'):
   os.makedirs('databases/mapper_data')
if not os.path.exists('databases/PFAM'):
   os.makedirs('databases/PFAM')
if not os.path.exists('databases/DBCAN'):
   os.makedirs('databases/DBCAN')   
if not os.path.exists('intermediate_files/BiG-SCAPE'):
   os.makedirs('intermediate_files/BiG-SCAPE')     
if not os.path.exists('databases/BiG-SCAPE'):
   os.makedirs('databases/BiG-SCAPE')     
if not os.path.exists('databases/BAKTA'):
   os.makedirs('databases/BAKTA')   


#Download PFAM

os.system('wget -P ./databases/PFAM/ ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/Pfam-A.hmm.gz')
os.system('gunzip databases/PFAM/Pfam-A.hmm.gz' )
os.system('hmmpress databases/PFAM/Pfam-A.hmm')

#Download EGGNOG

#os.system('download_eggnog_data.py -y --data_dir databases/mapper_data') broken for now

os.system('wget -P ./databases/mapper_data/ http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz')

os.system('wget -P ./databases/mapper_data/ http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz')

os.system('wget -P ./databases/mapper_data/ http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz')

os.system('gunzip ./databases/mapper_data/*.gz')

os.system('tar -xvf ./databases/mapper_data/eggnog.taxa.tar -C ./databases/mapper_data/')

#Download Bakta

os.system ('conda run -n bacLIFE_environment_BAKTA bakta_db download --output ./databases/mapper_data/BAKTA --type full')

#Download DBCAN

os.system('wget https://pro.unl.edu/dbCAN2/download_file.php?file=dbCAN-HMMdb-V14.txt --no-check-certificate -O ./databases/DBCAN/dbCAN-HMMdb-V14.txt')

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
os.system(cmd)






