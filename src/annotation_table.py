import argparse
parser = argparse.ArgumentParser(description = 'Clean annotatons to import to R')
parser.add_argument('-i', '--input', dest='input', type=str, 
                        help='filename of the output of rule keeg/cog')
parser.add_argument('-o', '--output', dest='output', type=str, 
                        help='filename of the output of rthis script')
args = parser.parse_args()


file = open(args.input)
f =open(args.output, 'w')
list_genes = []
list_id = []
list_description = []
for line in file:
    line.strip()
    a = line.split("\t")
    joined_string = ",".join(a[3:])
    list_description.append(joined_string)
    geneid = a[1]
    list_genes.append(geneid)
    id = a[2]
    list_id.append(id)

for i in range(len(list_id)):
    f.write(list_genes[i])
    f.write("\t")
    f.write(list_id[i])
    f.write("\t")
    f.write(list_description[i])
    f.write("\n")