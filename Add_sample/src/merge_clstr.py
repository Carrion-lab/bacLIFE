import re
import argparse
import os

##change name to same as in clstr

def change_name(name):
    elements = name.split('_')
    new_name = str(elements[0][2:] + '_' + elements[1] + elements[2])
    return(new_name)

##Arguments

parser = argparse.ArgumentParser(description = 'Merge two .clstr files')
parser.add_argument('-i1', '--input1', dest='input1',  type = str, help= 'filename of original .clstr file')
parser.add_argument('-i2', '--input2',dest='input2', type = str, help= 'filename of new .clstr file')
parser.add_argument('-o', '--output', dest='output',type = str, help= 'filename of merged .clstr file')
args = parser.parse_args()


###Open files

file2 = open(args.input1)
file1 = open(args.input2)

sample_name = os.listdir('data/')

f =open(args.output, 'w')




### Read new.clstr file, extract clusters and file lines to paste in old .clstr file

list_clusters_id= []
list_new_lines = []
cluster_id = 'none'
cluster_id2 = 'none'


for line in file1:
    if line[0] == '>':
        cluster_id = line
    for name in sample_name:
        new_name = change_name(name)
        if re.search(new_name, line):
            list_clusters_id.append(cluster_id)
            list_new_lines.append(line)
print(list_clusters_id[:15])
print(list_new_lines[:15])





#### Write new .clstr file merged

index_cluster_id_list = 0
seq_number = 'none'

for line in file2:
    f.write(line)
    if line[0] == '>':
        #new_line[1] = seq_number
        cluster_id2 = line
        if cluster_id2 in list_clusters_id:
            #index_cluster_id_list = list_clusters_id.index(cluster_id2)
            index_cluster_id_list = [i for i, x in enumerate(list_clusters_id) if x == cluster_id2]
            if len(index_cluster_id_list) > 1:
                for i in index_cluster_id_list:
                    f.write(list_new_lines[i])
            else:
                new_line = list_new_lines[index_cluster_id_list[0]]
                f.write(new_line)
    else:
        seq_number = line[0]







