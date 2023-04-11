#!/usr/bin/env python
import argparse, os
import Bio
from Bio import SeqIO

def genbank2faa(infile, outfile):
    seq = SeqIO.parse(infile,'genbank')
    fastaseq = []
    if not os.path.exists(outfile):
        print( "Converting :%s"%(infile) )
        #partial CDS_counter
        partial = 0
        
        id_2_tags = "%s.tags"%(outfile)
        
        tagfile = open(id_2_tags,'w')
        
        for s in seq:
            
            features = s.features
            for f in features:
                if f.type == "source":
                    sourcefeature = f
            
            if 'organism' in sourcefeature.qualifiers:
                organism = sourcefeature.qualifiers['organism'][0]
                strain = ''
                
                if 'strain' in sourcefeature.qualifiers:
                    strain = sourcefeature.qualifiers['strain'][0]
                    strain = strain.replace(' ','').replace(':','').replace('|','')
                organism = organism + '_' + strain
                
            #else:
            #    organism = s.annotations['organism'].replace(' ','_').replace(':','_').replace('|','_')
              
            organism = organism.replace(' ','_').replace(':','_').replace('|','_')
            
            ltprefix = estimate_locustag_prefix(features)
            ltprefixnew = "%s_"%(ltprefix)
                   
            
            moltype="C"
            for f in features:
                if f.type == "source":
                    if "plasmid" in f.qualifiers:
                        moltype="P"
            
            cdscount=0
            for f in features:
                if f.type == "CDS":
                    cdscount = cdscount + 1
                    if 'locus_tag' in f.qualifiers:
                        locus_tags = f.qualifiers['locus_tag']
                        locus_tag = locus_tags[-1]
                        if not '_' in locus_tag:
                           locus_tag = locus_tag.replace(ltprefix, ltprefixnew) 
                        
                        if not isinstance(f.location.start, Bio.SeqFeature.ExactPosition):
                            # location is not exact.
                            partial = partial +1
                            locus_tag = "%s_%02d"%(locus_tag, partial)
                                                
                        prot = f.extract(s)
                        prot.seq = prot.seq.translate(table=11, to_stop=True)
                        prot.id = locus_tag
                        
                        tagfile.write(locus_tag)
                        tagfile.write("\t")
                        tagfile.write(moltype)
                        tagfile.write("\t")
                        tagfile.write(organism)
                        tagfile.write("\t")
                        tagfile.write(prot.description)
                        tagfile.write("\n")
                        
                        if 'product' in f.qualifiers:
                            prot.description = f.qualifiers['product'][0]
                        else:
                            prot.description = 'unknown'
                                
                        fastaseq.append(prot)
            print( "Written :%s sequences from %s"%(len(fastaseq), s.description) )
        
        if partial>0:
            print("Warning sequence %s contains %d partial CDS sequences. This might inflate the number of clusters."%(infile, partial))
            
        SeqIO.write(fastaseq, outfile,'fasta')
    else:
        print( "Not Converting :%s"%(infile) )
        print( "%s exists"%(outfile) )
        
        
def genbank2ffn(infile, outfile):
    seq = SeqIO.parse(infile,'genbank')
    fastaseq = []
    print( "Converting :%s"%(infile) )
    #partial CDS_counter
    partial = 0
    
    
    for s in seq:
        features = s.features
        for f in features:
            if f.type == "gene":
                if 'locus_tag' in f.qualifiers:
                    locus_tags = f.qualifiers['locus_tag']
                    locus_tag = locus_tags[-1]
                    if not isinstance(f.location.start, Bio.SeqFeature.ExactPosition):
                        # location is not exact.
                        partial = partial +1
                        locus_tag = "%s_%02d"%(locus_tag, partial)
                    
                    gene = f.extract(s)
                    gene.id = locus_tag
                    if 'gene' in f.qualifiers:
                        gene.description = f.qualifiers['gene'][0]
                    else:
                        gene.description = gene.id
                        
                    fastaseq.append(gene)
        print( "Written :%s sequences from %s"%(len(fastaseq), s.description) )
    
    if partial>0:
        print("Warning sequence %s contains %d partial CDS sequences. This might inflate the number of clusters."%(infile, partial))
        
    SeqIO.write(fastaseq, outfile,'fasta')

def estimate_locustag_prefix(features):
    
    strings = [f.qualifiers['locus_tag'][-1] for f in features if 'locus_tag' in f.qualifiers]
    
    """ Find the longest string that is a prefix of all the strings.
    """
    if not strings:
        return ''
    prefix = strings[0]
    for s in strings:
        if len(s) < len(prefix):
            prefix = prefix[:len(s)]
        if not prefix:
            return ''
        for i in range(len(prefix)):
            if prefix[i] != s[i]:
                prefix = prefix[:i]
                break
    return prefix
            
if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genbank', dest = 'genbank', help = 'genbank input file', required = True)
    parser.add_argument('-f', '--fasta', dest = 'fasta', help = 'fasta CDS output file', required = True)

    args = parser.parse_args()
    genbank2faa(args.genbank, args.fasta)
