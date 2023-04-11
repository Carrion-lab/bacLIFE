import fnmatch
import os, glob
from Bio import SeqIO
import sys, traceback
import argparse

def rename_fasta(infile, tagfile, outfile, newtagfile):
    counters = {}
    tags = {}
    newseqs=[]
    
    t = open(tagfile,'r')
    tn = open(newtagfile,'w')
    
    for line in t:
        data = line.strip().split('\t')
        tags[data[0]] = data
        
    seqs = SeqIO.parse(infile, 'fasta')
    for s in seqs:
        data = tags[s.id]
        counter = counters.get(data[2], 0)
        counter = counter + 1
        idtag = "%s|%06d"%(data[2], counter)
        counters[data[2]]=counter
        s.id = idtag
        s.description = ""
        tn.write(idtag)
        tn.write('\t')
        tn.write('\t'.join(data))
        tn.write('\n')
        newseqs.append(s)
    tn.close()
    
    try:
        SeqIO.write((newseqs), outfile, 'fasta')
    except:
        print("Exception in user code:")
        print("-"*60)
        traceback.print_exc(file=sys.stdout)
        print("-"*60)


if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', dest = 'fasta', help = 'fasta input file', required = True)
    parser.add_argument('-t', '--tags', dest = 'tags', help = 'tags', required = True)
    parser.add_argument('-o', '--outfile', dest = 'outfile', help = 'outfile', required = True)
    parser.add_argument('-n', '--newtags', dest = 'newtags', help = 'tags', required = True)
    
    args = parser.parse_args()
    rename_fasta(args.fasta, args.tags, args.outfile, args.newtags)
    