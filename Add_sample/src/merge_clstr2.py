import argparse
import os

parser = argparse.ArgumentParser(description = 'Merge two .clstr files')
parser.add_argument('-i1', '--input1', dest='input1',  type = str, help= 'filename of temp .clstr file')
parser.add_argument('-i2', '--input2',dest='input2', type = str, help= 'filename of new sequences .clstr file')
parser.add_argument('-o', '--output', dest='output',type = str, help= 'filename of merged .clstr file')
args = parser.parse_args()

file1 = open(args.input1)
file2 = open(args.input2)
f =open(args.output, 'w')

###Detect max id of clusters
cluster_id = 'none'

for line in file1:
    f.write(line)
    if line[0] == '>':
        cluster_id = line

cluster_max = cluster_id.split()
cluster_max = int(cluster_max[1])


print( 'N clusters original clstr file: ' + cluster_id)

cluster_id = 'none'
new_cluster_id = 'none'

list_clusters_id= []
n = 1
for line in file2:
    if line[0] == '>':
        cluster_id = line
        list_clusters_id.append(cluster_id)
        cluster_id_list = cluster_id.split()
        new_cluster_id = str(cluster_id_list[0] + ' ' + str(cluster_max + n))
        n += 1
        f.write(new_cluster_id)
        f.write('\n')
    else:
        f.write(line)
