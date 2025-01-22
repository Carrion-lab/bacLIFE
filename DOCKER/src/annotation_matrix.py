#!/usr/bin/env python

#Makes absence presence matrix with the KO  and COG annotations


#Author: Dennis Jansen / Victor de Jager
#Contact: dennisj4995@gmail.com // D.Jansen@nioo.knaw.nl
#Latest edit:15-12-2016


import argparse, os, itertools, re
from collections import defaultdict
from Bio import SeqIO
class KO2matrix():
    species = defaultdict(list)
    KO = defaultdict(list)
    COG = defaultdict(list)##
    species1=defaultdict(list)##
    
    def parse_fasta(self,fastafile):
    # parse through the cluster fasta sequences and placed them in a dict
        self.seqdict = SeqIO.to_dict(SeqIO.parse(fastafile,'fasta'))
            
    def parse_uproc(self, csvfile):
        
        fh = open(csvfile,'r')
        for line in fh:
            data = line.strip().split(',')
            sp = data[1].split('|')[0]
            ko=data[3]
            self.species[sp].append(data[1])
            self.KO[ko].append(data[1])
    def parse_cog(self,cog):
        cog = open(cog,'r')
        for line in cog:
            data= line.strip().split(',')
            splt= data[1].split('|')[0]
            cog=data[3]
            self.species1[splt].append(data[1])
            self.COG[cog].append(data[1])
        
    def write_abspres_matrix(self, abspres_matrix):
        f =open(abspres_matrix,'w')
        
        specieslist = sorted(self.species.keys())
        f.write('KO\t')
        f.write("\t".join(specieslist))
        f.write("\n")
        
        kos = sorted(self.KO.keys())
        for k in kos:
            
            f.write(k)
            
            genes = self.KO[ k ]
            species = [ g.split('|')[0] for g in genes]
            
            for s in specieslist:
                if s in species:
                    f.write("\t%s"%('1'))
                else:
                    f.write("\t0")
                    
            f.write("\n")
        f.close()
    def write_abspres_cog(self, abspres_cog):
        fc =open(abspres_cog,'w')
        
        species1list=sorted(self.species1.keys())
        fc.write('COG\t')
        fc.write('\t'.join(species1list))
        fc.write('\n')
        cogs = sorted(self.COG.keys())
        for c in cogs:
            fc.write(c)
            genes = self.COG[ c ]
            species= [ g.split('|')[0] for g in genes]
            for s in species1list:
                if s in species:
                    fc.write('\t%s'%('1'))
                else:
                    fc.write('\t0')
            fc.write('\n')
        fc.close()
def main():
    # set argument parser
    parser = argparse.ArgumentParser(description = 'Create orthomcl like output from the KO annotation of combined protein files')
    parser.add_argument('-f', '--fasta', metavar='.fasta file', dest='fasta', type=str, 
                        help='The .fasta file containing the combined proteins')
    parser.add_argument('-c', '--csv', metavar='.csv file from uproc KO annotation', dest='csv', type=str,
                        help='The .clstr file produced by CD-hit that contains the cluster information.')
    parser.add_argument('-o', '--kofile', metavar='kofile', dest='kofile', type=str,
                        help='The file that will contain the KO matrix')
    parser.add_argument('-x', '-cogfile', metavar='.csv file from uproc cog annotation', dest='cog',type=str,
                        help='The file that will contain the COG matrix')
    parser.add_argument('-y', '--absence_cog_matrix',metavar='absence_cog_matrix', dest='absence_cog_matrix', type=str, default='absence_cog.txt',
                        help='absence_cog matrix file')
    parser.add_argument('-a', '--absence_presence_matrix', metavar='absence_presence_matrix', dest='absence_presence_matrix', type=str, default='absence_presence.txt',
                        help='absence_presence matrix file')
    
    args = parser.parse_args()
    converter = KO2matrix()
    converter.parse_uproc(args.csv)
    converter.parse_cog(args.cog)
    converter.write_abspres_matrix(args.absence_presence_matrix)
    converter.write_abspres_cog(args.absence_cog_matrix)
    
if __name__ == '__main__':
    main()
