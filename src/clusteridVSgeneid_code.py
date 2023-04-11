import argparse
import re
parser = argparse.ArgumentParser(description = 'Create a file matching cluster number and geneid')
parser.add_argument('-i', '--input', dest='input', type=str, 
                        help='filename of the .clstr produced by cd-hit')
parser.add_argument('-o', '--output', dest='output', type=str, 
                        help='filename of the output of this script')
args = parser.parse_args()



file = open(args.input)
f =open(args.output, 'w')

geneids = []
clusterids = []
for line in file:
	line = line.strip()
	if 'Cluster' in line:
		print(line)
		id = line[1:]
		clusterids.append(id)
	if '*' in line:
		id = re.search('>(.*).....', line)
		geneids.append(id.group(1))
	



for i in range(len(geneids)):
	f.write(clusterids[i])
	f.write("\t")
	f.write(geneids[i])
	f.write("\n")