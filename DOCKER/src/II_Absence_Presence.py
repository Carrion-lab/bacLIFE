############################ 
### John Baker Hernandez ###
### johnbakerh@gmail.com ###
###        2020          ###
###  Leiden University   ###
### Bioinformatics Group ###
############################
'Python Script used to create an absence presence table of GCFs per Strain from BiG-SCAPE file'

#Import Libraries
import pandas as pd
import sys
import os
from string import digits

#Importing DataFrame
df1 = pd.read_csv(sys.argv[1], sep="\t")
df1.columns = ['BGC Name', 'GCF No', 'Organism', 'BGC Class']
#print(df1.isnull())
#print(df1.head())

'Method for Creating Absene Presence Table'
#Dataframe Copies
df2 = df1.copy() #Copy Original DataFrame
df2.columns = ['BGC Name', 'GCF No', 'Organism', 'BGC Class']
#DataFrame Manipulation
df2['GCF No'] = df1['BGC Class'].astype(str) + 'GCF' + df2['GCF No'].astype(str)
df1[['BGC','BGC2']] = df1['BGC Name'].str.split('_',n=1,expand=True)
df2['Genome'] = df1['Organism'].astype(str) + df1['BGC Name'].astype(str)
df2.set_index('GCF No')
df2 = df2[['GCF No','Genome']]
df2['Genome'] = df2['Genome'].str.replace(' ', '_').str.replace('Unclassified','_').str.replace('__','_').str.replace('.','')
path = os.path.join(os.path.abspath(os.getcwd()), 'intermediate_files/BiG-SCAPE')
df2.to_csv(os.path.join(path,r'abs_pres_table.csv'))


# 'Alternative Method'
# #Initializing the new_df (df2) dataframe
# import pandas as pd
# new_columns =  ['GCF', 'Genome']
# new_df  = pd.DataFrame(columns = new_columns)
# #new_df
# #df1.columns
# #creating the column GCF of the new dataframe
# GCF = [str(df1['BGC Class'][i])+ '_GCF' + str(df1['GCF No'][i]) for i in range(len(df1))]
# #creating the new column 'Genome' of the new dataframe
# organism = [df1['Organism'][i].replace(' ', '_') for i in range(len(df1))] #substituting ' ' by '_'
# BGC_name = [list(df1['BGC Name'].str.split('_'))[i][0] for i in range(len(df1))]
# genome = [organism[i] + '_' + BGC_name[i] for i in range(len(df1))]
# new_df['GCF'] = GCF
# new_df['Genome'] = genome
