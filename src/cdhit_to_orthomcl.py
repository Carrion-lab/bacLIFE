#!/usr/bin/env python

# Author: Victor de Jager
# Contact: v.dejager@nioo.knaw.nl

# Create orthomcl like output from the cluster output from CD-hit

# command: cdhit_parser.py -f [cluster .fasta file] -c [.clstr cluster file]

### rewrite
# the output of the CD-hit_filter.py script is a fasta file [input cluster .fasta file with the min # seqs in the file name]
# with the OTU sequences for the clusters with more than or equal number of sequences to the user specified minimum.

# import modules used by the script
import argparse, os, itertools, re, numpy
from Bio import SeqIO

class convert_cdhit_to_orthomcl():
    cluster_dict = {}
    margin = 0.05
    def read_clstr(self, infile):
        
        # parse through the .clstr file and create a dictionary
        # with the sequences per cluster
        m = re.compile(r'(\d+)\s+(\d+).+>(.+)\.\.\.\s(.+)')
        # open the cluster file and set the output dictionary
        cluster_file = open(infile,'r')
                # parse through the cluster file and store the cluster name + sequences in the dictionary
        cluster_groups = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))
        for cluster in cluster_groups:
            name = next(cluster)
            name = name.strip().strip('>').replace(' ','_')
            count = int(name.split('_')[1])
            name = "cluster_%05d"%count
            #seqs = [ seq.split('>')[1].split('...')[0] for seq in cluster_groups.next() ]
            
            cluster_groups_members = list( next(cluster_groups) )
            
            seqs = [ m.search(seq).group(3) for seq in  cluster_groups_members]
            sizes = [ int(m.search(seq).group(2)) for seq in cluster_groups_members ]
            
            self.cluster_dict[name] = {'members' : seqs, 'sizes':sizes, 'complete': (self.margin > (numpy.array(sizes).std()/numpy.array(sizes).mean() ) )}
    
    def get_strains(self,genes):
        strains = [ s.split('|')[0] for s in genes]
        return strains
    
    def cluster_stat(self):
        
        self.strainset = {}
        clusters = self.cluster_dict.keys()
        for k in clusters:
            genes = self.cluster_dict[ k ]['members']
            #print (k,self.cluster_dic[ k ]['complete'], len(genes))
            
            strains = self.get_strains(genes)
            for s in strains:
                if s in self.strainset:
                    self.strainset[s].add(k)
                else:
                    self.strainset[s]=set([k])

        core = set(clusters)    
        for k in self.strainset:
            core = core.intersection(self.strainset[k])
            print (k, len(self.strainset[k]), len(core))
        
        singletons = 0
        for k in self.cluster_dict.keys():
            genes = self.cluster_dict[ k ]['members']
            if len(genes)<2:
                singletons=singletons+1
        
        
        print("Clusters: %s"%len(clusters))
        print("Clusters with two or more proteins: %s"%(len(clusters)-singletons))
        print("Singleton clusters: %s"%(singletons))
        
    def write_og_file(self, ogfile, ogpath):
        
        f =open(ogfile,'w')
        
        if not (os.path.exists(ogpath)):
            try:
                os.makedirs(ogpath)
            except:
                print("path %s cannot be created"%(ogpath))
                pass
        
        for k in self.cluster_dict.keys():
            ogseqs = []
            data = []
            genes = self.cluster_dict[ k ]['members']
            status = self.cluster_dict[ k ]['complete']
            strains = strains = self.get_strains(genes)
            strain_count = len(set(strains))
            descriptions = [ self.seqdict[ s ].description for s in genes ]
            descriptions = set(descriptions)
            
            data.append(k)
            data.append(';'.join(genes))
            data.append(';'.join(descriptions))
            
            f.write("\t".join(data))
            f.write("\n")
            
            ogseqs = [ self.seqdict[s] for s in genes ]
            if (os.path.exists(ogpath)):
                fastafile = '%s/%s.%s.%s.fasta'%(ogpath,k,str(status),strain_count)
                SeqIO.write(ogseqs, fastafile, 'fasta' )
            else:
                print("Cannot write %s"%(fastafile))
        f.close()
        
    def write_abspres_matrix(self,abspres_matrix):
        f =open(abspres_matrix,'w')
        
        strains = sorted(self.strainset.keys())
        
        f.write("clusters\t")
        f.write("\t".join(strains))
        f.write("\tcompleteness\tdescriptions")
        f.write("\n")
        
        clusters = sorted(self.cluster_dict.keys())
        for k in clusters:
            
            f.write(k)
            f.write("\t")
            genes = self.cluster_dict[ k ]['members']
            clusterstrains = self.get_strains(genes)
            for s in strains:
                if s in clusterstrains:
                    f.write("%s\t"%(clusterstrains.count(s)))
                else:
                    f.write("0\t")
                    
            f.write("%s\t"%(self.cluster_dict[ k ]['complete']))
            descriptions = [ self.seqdict[sq].description for sq in genes ]
            descriptions = set(descriptions)
            f.write(';'.join(descriptions))
            f.write("\n")
        f.close()
           
            
        
    def parse_fasta(self,fastafile):
        # parse through the cluster fasta sequences and placed them in a dict
        self.seqdict = SeqIO.to_dict(SeqIO.parse(fastafile,'fasta'))
        


def main():
    # set argument parser
    parser = argparse.ArgumentParser(description = 'Create orthomcl like output from the cluster output from CD-hit')
    parser.add_argument('-f', '--fasta', metavar='.fasta file', dest='fasta', type=str, 
                        help='The .fasta file containing the clusters produced by CD-hit.')
    parser.add_argument('-c', '--clusters', metavar='.clstr file', dest='clusters', type=str,
                        help='The .clstr file producec by CD-hit that contains the cluster information.')
    parser.add_argument('-o', '--ogfile', metavar='ogfile', dest='ogfile', type=str,
                        help='The file that will contain the orthologous groups')

    parser.add_argument('-p', '--path', metavar='path', dest='path', type=str, default='ogs',
                        help='path for ogfiles')
    parser.add_argument('-a', '--absence_presence_matrix', metavar='absence_presence_matrix', dest='absence_presence_matrix', type=str, default='absence_presence.txt',
                        help='absence_presence matrix file')
    
    args = parser.parse_args()

    converter = convert_cdhit_to_orthomcl()
    converter.read_clstr(args.clusters)
    converter.parse_fasta(args.fasta)
    converter.cluster_stat()
    converter.write_og_file(args.ogfile, args.path)
    converter.write_abspres_matrix(args.absence_presence_matrix)
    
    
if __name__ == '__main__':

    main()
