GENUS, SPECIES, STR, REPLICON = glob_wildcards("data/{genus}_{species}_{str}_{replicon}.fna")
TEMP_DIR = 'intermediate_files/combined_proteins'
COMBINED_PROTEINS = os.path.join(TEMP_DIR,'combined_proteins.fasta')


COMPARISON_DIR = "intermediate_files/clustering"
CLUSTERING_BINARY_TABLE = 'intermediate_files/clustering/binary_matrix.txt'
CLUSTERING_FASTA = 'intermediate_files/clustering/protein_cluster'

configfile: "config.json"
THREADS = config['threads']
THREADS_prokka_antismash = config['threads_prokka_antismash']
MCL_INFLATION = config['mcl_inflation_value']
LINCLUST_IDENTITY = config['linclust_identity']





GENBANKFILES, = glob_wildcards("intermediate_files/annot/{genome}.gbk")
NEWTAGFILE = "intermediate_files/combined_proteins/id2tags.tsv"

CLUSTERVSGENE = "intermediate_files/clusterVSgene.txt"
PROKKA_ANNOTATION = "intermediate_files/prokka_annotation.txt"
PFAM_ANNOTATION = "intermediate_files/pfam_out/PFAM_annotation.txt"
KEGG_ANNOTATION = "intermediate_files/kegg_annotation/KEGG_annotation_clean.faa.finalkegg"
KEGG_ANNOTATION_CLEAN = "intermediate_files/kegg_annotation/KEGG_annotation_clean2.faa.finalkegg"
KEGG_DESCRIPTIONS = "src/ko_description.txt"
COG_ANNOTATION = "intermediate_files/cog_annotation/COG_annotation_clean.faa.finalcog"
DBCAN_ANNOTATION = "intermediate_files/dbcan_annotation/DBCAN_annotation.txt"
EGGNOG_ANNOTATION = 'intermediate_files/eggnog_annotation/eggnog_annotation.emapper.annotations'
EGGNOG_DATA = 'intermediate_files/mapper_data'
HMM_ANNOTATIONS = "intermediate_files/hmm_annotations.txt"
MEGAMATRIX = "MEGAMATRIX.txt"
MAPPING_FILE = 'mapping_file.txt'
CORRECTED_MAPPING_FILE = 'corrected_mapping_file.txt'
BIGSCAPE = "intermediate_files/BiG-SCAPE/bigscape_output/index.html"

rule final:
        input:
            prokka = expand('intermediate_files/annot/{species}_{str}_{replicon}/{genus}_{species}_{str}_{replicon}.gbk', zip, genus = GENUS, species = SPECIES, str = STR, replicon = REPLICON),
            extprot = expand('intermediate_files/annot/{species}_{str}_{replicon}/{genus}_{species}_{str}_{replicon}.ext_prot.faa', zip, genus = GENUS, species = SPECIES, str = STR, replicon = REPLICON),
            protein_combine = 'intermediate_files/combined_proteins/combined_proteins.fasta.orig',
            protein_rename = COMBINED_PROTEINS,
            clustering_matrix = CLUSTERING_BINARY_TABLE,
            clustering_fasta = CLUSTERING_FASTA,
            extract_prokka = PROKKA_ANNOTATION,
            pfam_annotation = PFAM_ANNOTATION,
            kegg_clean_annotation = KEGG_ANNOTATION,
            cog_clean_annotation = COG_ANNOTATION,
            dbcan_annotation = DBCAN_ANNOTATION,
            hmm_annotations = HMM_ANNOTATIONS,
            megamatrix = MEGAMATRIX,
            antismash = expand('intermediate_files/antismash/{genus}_{species}_{str}_{replicon}/{genus}_{species}_{str}_{replicon}.gbk', zip, genus = GENUS, species = SPECIES, str = STR, replicon = REPLICON),
            bigscape_setup = "intermediate_files/BiG-SCAPE/bigscape.py",
            bigscape = BIGSCAPE,
            binary_table_GCF = 'intermediate_files/BiG-SCAPE/big_scape_binary_table.txt',
            rename_matrix = 'MEGAMATRIX_renamed.txt',
            phylophlan = "intermediate_files/phylophlan/output_phylophlan/RAxML_bestTree.input_refined.tre"
            
rule directories:
        input:
            expand("data/{genus}_{species}_{str}_{replicon}.fna", zip, genus = GENUS, species = SPECIES, str = STR, replicon = REPLICON)
        params:
            annot = "intermediate_files/annot/",
            antismash = "intermediate_files/antismash/",
            eggnog = "intermediate_files/eggnog_annotation/"
        output:
            chk = ".mkdir.chkpnt"
        run:
            shell("mkdir -p {params}")
            shell("touch .mkdir.chkpnt")

rule Prokka_annotation:
        input:
            file = "data/{genus}_{species}_{str}_{replicon}.fna",
            dir = rules.directories.output
        output:
            prokka = 'intermediate_files/annot/{species}_{str}_{replicon}/{genus}_{species}_{str}_{replicon}.gbk',
            prokka2 = 'intermediate_files/annot/{species}_{str}_{replicon}/{genus}_{species}_{str}_{replicon}.faa'
        message: 'executing prokka.'
        params:
            genus = "{genus}",
            species = "{species}",
            str = "{str}",
            outdir = "intermediate_files/annot/{species}_{str}_{replicon}",
            prefix = "{genus}_{species}_{str}_{replicon}"
        threads: config['threads_prokka_antismash']
        priority: 100
        run:
                shell('prokka --force --outdir {params.outdir} --prefix {params.prefix} --locustag {params.str} --addgenes --increment 5 --centre NIOO-KNAW --genus {params.genus} --species {params.species} --str {params.str} --gcode 11 --cpus {THREADS_prokka_antismash} --evalue 1e-03 --rfam {input.file}')

rule extract_proteins:
    input: rules.Prokka_annotation.output.prokka
    output:
        faa = 'intermediate_files/annot/{species}_{str}_{replicon}/{genus}_{species}_{str}_{replicon}.ext_prot.faa',
        tags = 'intermediate_files/annot/{species}_{str}_{replicon}/{genus}_{species}_{str}_{replicon}.ext_prot.faa.tags'
    #log: 'log/{genus}_{species}_{str}_{replicon}_extract_proteins.log'
    #benchmark: 'benchmarks/{genus}_{species}_{str}_{replicon}_extract_proteins.json'
    message: 'Executing genbank_to_protein_fasta on the following files {input}.'
    shell:
        'python ./src/genbank_to_protein_fasta.py --genbank {input} --fasta {output.faa}'

rule protein_combine:
    input:
        #faa = expand('intermediate_files/annot/{genbankfile}.ext_prot.faa',  genbankfile=GENBANKFILES),
        #tags = expand('intermediate_files/annot/{genbankfile}.ext_prot.faa.tags', genbankfile=GENBANKFILES)
        faa = expand(rules.extract_proteins.output.faa, zip, genus = GENUS, species = SPECIES, str = STR, replicon = REPLICON),
        tags = expand(rules.extract_proteins.output.tags, zip, genus = GENUS, species = SPECIES, str = STR, replicon = REPLICON)
    output:
        orig = "%s.orig"%COMBINED_PROTEINS,
        tags = "%s.tags"%COMBINED_PROTEINS
    log: 'log/protein_combine.log'
    run:
        shell('cat {input.faa} >  {output.orig}')
        shell('cat {input.tags} >  {output.tags}')

rule protein_rename:
    input:
        orig = rules.protein_combine.output.orig,
        tags = rules.protein_combine.output.tags
    output: 
        combined_proteins = COMBINED_PROTEINS,
        newtagfile = NEWTAGFILE
    run:
        shell('python ./src/rename_proteins.py -f {input.orig} -t {input.tags} -o {output.combined_proteins} -n {output.newtagfile}')



rule clustering:
    input: COMBINED_PROTEINS
    output: 
        binary_table = CLUSTERING_BINARY_TABLE,
        fasta = CLUSTERING_FASTA
    params:
        mmseq_identity = LINCLUST_IDENTITY,
        mmseq_db = 'intermediate_files/clustering/mmseqDB',
        mmseq_db_clu = 'intermediate_files/clustering/mmseqDB_clu',
        mmseq_temp = 'intermediate_files/clustering/mmseqDB_temp',
        mmseq_tsv = 'intermediate_files/clustering/mmseq_tsv.tsv',
        mmseq_rep = 'intermediate_files/clustering/mmseq_clurep',
        mmseq_fasta = 'intermediate_files/clustering/mmseq_clurep_fasta.fasta',
        diamond_db = 'intermediate_files/clustering/diamond_db',
        diamond_unaligned = 'intermediate_files/clustering/unaligned.fasta',
        diamond_unaligned_headers = 'intermediate_files/clustering/unaligned_headers.txt',
        diamond_results = 'intermediate_files/clustering/diamond_results',
        diamond_results_filtered = 'intermediate_files/clustering/diamond_results_filtered',
        diamond_results_clean = 'intermediate_files/clustering/diamond_results_clean',
        mcl_data = 'intermediate_files/clustering/mcl_data.mci',
        mcl_tab = 'intermediate_files/clustering/mcl_tab.tab',
        mcl_infaltion = MCL_INFLATION,
        mcl_clusters= 'intermediate_files/clustering/out.mcl_data.mci',
        
    run:
        #Extract gene lengths
        shell("bioawk -c fastx '{{ print $name, length($seq) }}' < intermediate_files/combined_proteins/combined_proteins.fasta >intermediate_files/combined_proteins/length_genes.txt")
        
        #Run mmseq2 clustering to 0.95
        shell('mmseqs createdb {input} {params.mmseq_db}')
        shell('mmseqs cluster {params.mmseq_db} {params.mmseq_db_clu} {params.mmseq_temp} --min-seq-id {params.mmseq_identity} --cov-mode 0 -c 0.8')
        shell('mmseqs createtsv {params.mmseq_db} {params.mmseq_db} {params.mmseq_db_clu} {params.mmseq_tsv}')
        
        #Extract mmseq2 representatives
        shell('mmseqs createsubdb {params.mmseq_db_clu} {params.mmseq_db} {params.mmseq_rep}')
        shell('mmseqs convert2fasta {params.mmseq_rep} {params.mmseq_fasta}')
        
        
        #Run diamond
        shell('diamond makedb --in {params.mmseq_fasta} -d {params.diamond_db}')
        shell('diamond blastp -b 10 -p {THREADS} -c 1 -e 0.00001 --un {params.diamond_unaligned} -k 5000 -d {params.diamond_db} -q {params.mmseq_fasta} -o {params.diamond_results} --outfmt 6 qseqid sseqid pident bitscore')
        shell("awk -F'\t' '$3>20' {params.diamond_results} > {params.diamond_results_filtered}")
        shell("cut -f1,2,4 {params.diamond_results_filtered} > {params.diamond_results_clean}")
        
        ##Get unaaligned fasta headers
        shell('grep "^>" {params.diamond_unaligned}  > {params.diamond_unaligned_headers}')
        
        ### Run MCL
        shell('mcxload -abc {params.diamond_results_clean} --stream-mirror -o {params.mcl_data} -write-tab {params.mcl_tab}')
        shell('mcl {params.mcl_data} -I {params.mcl_infaltion} -use-tab {params.mcl_tab} -o {params.mcl_clusters}')
        #shell('mv {params.mcl_clusters} intermediate_files/clustering/')
        
        ##Generate matrix

        shell('Rscript src/MCL_merge.R {params.mmseq_tsv} {params.mcl_clusters} {params.diamond_unaligned_headers} {output.binary_table} {input} {output.fasta}')




rule extract_prokka:
    input:
        id2tag = rules.protein_rename.output.newtagfile,
        faa = rules.protein_combine.output.orig
    params:
        headers = 'intermediate_files/prokka_headers.txt'
    output:
        "intermediate_files/prokka_annotation.txt"
    run:
        shell("grep '^>' {input.faa} >  {params.headers}")
        shell("Rscript ./src/extract_prokka.R {input.id2tag} {params.headers} {output} ")




rule pfam:
    input:
            rules.clustering.output.fasta
    output:
            pfam = PFAM_ANNOTATION
            
    message: 'executing pfam.'
    run:
        shell('hmmsearch --tblout {output.pfam} --cpu {THREADS} -E 1e-5 ./intermediate_files/PFAM/Pfam-A.hmm {input}')

rule EGGNOG:
    input:
        rules.clustering.output.fasta
    params:
        data = EGGNOG_DATA
    output:
        EGGNOG_ANNOTATION
    run:
        shell('emapper.py -i {input} --cpu {THREADS} -o intermediate_files/eggnog_annotation/eggnog_annotation --data_dir {params.data} --pident 30 --query_cover 50 --subject_cover 50 --report_orthologs')

rule KEGG_COG:
    input:
        rules.EGGNOG.output
    params:
        kegg_descriptions = KEGG_DESCRIPTIONS,
        cog_descriptions = "./src/cog_annotation_groups.csv",
    output:
        clean_annotation = 'intermediate_files/eggnog_annotation/eggnog_annotation.emapper.clean.annotations',
        cog = COG_ANNOTATION,
        kegg = KEGG_ANNOTATION,
        kegg_clean = KEGG_ANNOTATION_CLEAN
    run:
        shell("sed '/^#/d' {input} > {output.clean_annotation}")
        shell('Rscript src/COG_KEGG_annotations.R {output.clean_annotation} {params.cog_descriptions} {params.kegg_descriptions} {output.cog} {output.kegg}')
        shell("sed '/^ko/d' {output.kegg}  > {output.kegg_clean}")


rule dbCAN:
        input:
            rules.clustering.output.fasta
        message: 'Retrieving dbCAN annotations.'
        output:
                dbcan = DBCAN_ANNOTATION
        run:
            shell('hmmsearch --tblout {output} -E 1e-5 --cpu {THREADS} ./intermediate_files/DBCAN/dbCAN-HMMdb-V9.txt {input}')

rule process_hmm_annotations:
    input:
        dbcan = rules.dbCAN.output.dbcan,
        pfam = rules.pfam.output.pfam
    params:
        dbcan_family = "./src/CAZyDB.07302020.fam-activities.txt",
        pfam_family = "./src/Pfam-A.clans.tsv"
    output: HMM_ANNOTATIONS
    run:
        shell("Rscript ./src/Process_hmm_annotation.R {input.dbcan} {input.pfam} {params.dbcan_family} {params.pfam_family} {output}")


rule process_annotations:
    input:
        matrix = rules.clustering.output.binary_table,
        cog = rules.KEGG_COG.output.cog,
        kegg = rules.KEGG_COG.output.kegg_clean,
        hmm_annotation = rules.process_hmm_annotations.output,
        prokka = rules.extract_prokka.output
    output:
        MEGAMATRIX
    run:
        shell("set -euo pipefail; Rscript ./src/Process_annotations.R {input.matrix} {input.cog} {input.kegg} {input.hmm_annotation} {input.prokka} {output}")



rule antismash:
        input:
                rules.Prokka_annotation.output.prokka
        output:
                "intermediate_files/antismash/{genus}_{species}_{str}_{replicon}/{genus}_{species}_{str}_{replicon}.gbk"
        params:
            out_dir = 'intermediate_files/antismash/{genus}_{species}_{str}_{replicon}/',
            threads = THREADS
        threads: config['threads_prokka_antismash']
        conda:
            "antismash_bacLIFE"
        shell:
            'antismash --cpus {THREADS_prokka_antismash} --cb-general --cb-knownclusters --cb-subclusters --output-dir {params.out_dir} --asf --pfam2go --genefinding-tool prodigal --smcog-trees {input}'



rule bigscape_exe:
        input: 
            antismash = expand(rules.antismash.output, zip, genus = GENUS, species = SPECIES, str = STR, replicon = REPLICON),
            pfam_hmm = 'intermediate_files/PFAM/Pfam-A.hmm'
        output:
            html = 'intermediate_files/BiG-SCAPE/bigscape_output/index.html',
            clustering = 'intermediate_files/BiG-SCAPE/bigscape_output/network_files/hybrids_glocal/mix/mix_clustering_c0.70.tsv',
            network = 'intermediate_files/BiG-SCAPE/bigscape_output/network_files/hybrids_glocal/mix/mix_c0.70.network',
            annotations = 'intermediate_files/BiG-SCAPE/bigscape_output/network_files/hybrids_glocal/Network_Annotations_Full.tsv'
        threads: THREADS
        params:
            outdir = 'intermediate_files/BiG-SCAPE/bigscape_output/',
            threads = THREADS,
            indir = rules.directories.params.antismash
        conda:
            "bigscape_bacLIFE"
        shell:
            "python ./intermediate_files/BiG-SCAPE/bigscape.py -i {params.indir} -o {params.outdir} --pfam_dir intermediate_files/PFAM/ --mode glocal --mibig --cutoffs 0.3 0.7 --include_singletons --cores {params.threads} --mix; rm -r intermediate_files/BiG-SCAPE/bigscape_output/network_files/hybrids_glocal; mv intermediate_files/BiG-SCAPE/bigscape_output/network_files/*hybrids_glocal intermediate_files/BiG-SCAPE/bigscape_output/network_files/hybrids_glocal"

rule extract_binary_table_GCF:
    input:
        clustering = rules.bigscape_exe.output.clustering,
        network = rules.bigscape_exe.output.network,
        annotations = rules.bigscape_exe.output.annotations
    params:
        output_code_I_network = "intermediate_files/BiG-SCAPE/mix_filtered.network",
        output_code_I_annotations = 'intermediate_files/BiG-SCAPE/GCF_annotation.txt',
        output_code_II = 'intermediate_files/BiG-SCAPE/abs_pres_table.csv',
        names = 'names_equivalence.txt'
    output:
        filtered_network = "intermediate_files/BiG-SCAPE/mix_filtered.network",
        annotations = 'intermediate_files/BiG-SCAPE/annotation.txt',
        merged_annotations = 'intermediate_files/BiG-SCAPE/GCF_annotation.txt',
        abs_presence_list = 'intermediate_files/BiG-SCAPE/abs_pres_table.csv',
        binary_matrix = 'intermediate_files/BiG-SCAPE/big_scape_binary_table.txt'
    run:
        shell("Rscript src/I-BIGSCAPE_revision.R {input.clustering} {input.network} {input.annotations} intermediate_files/BiG-SCAPE/mix_filtered.network intermediate_files/BiG-SCAPE/GCF_annotation.txt {output.annotations} {params.names}")
        shell("python src/II_Absence_Presence.py intermediate_files/BiG-SCAPE/GCF_annotation.txt intermediate_files/BiG-SCAPE/abs_pres_table.csv ")
        shell("Rscript src/III_Absence_Presence_GCF.R intermediate_files/BiG-SCAPE/abs_pres_table.csv {output.binary_matrix}")



rule rename_MEGAMATRIX:
    input:
        genes = rules.process_annotations.output,
        BGCs = rules.extract_binary_table_GCF.output.binary_matrix
    output:
        genes = 'MEGAMATRIX_renamed.txt',
        BGCs = 'intermediate_files/BiG-SCAPE/big_scape_binary_table_renamed.txt',
        mapping_file = 'mapping_file.txt'
    params:
        'names_equivalence.txt'
    run:
        shell('Rscript src/rename_MEGAMATRIX.R {input.genes} {input.BGCs} {params} {output.genes} {output.BGCs} {output.mapping_file}')
        


rule phylophlan:
    input:
        config = "src/supermatrix_aa.cfg",
        to_order= rules.extract_binary_table_GCF.output.binary_matrix
    params:
        database = config['phylo_database'],
        in_file = "intermediate_files/phylophlan/input"
    output:
        out_tree = "intermediate_files/phylophlan/output_phylophlan/RAxML_bestTree.input_refined.tre",
        #out_dir = "intermediate_files/phylophlan/output_phylophlan"
    log: "log/phylophlan.log"
    run:
        shell("mkdir -p intermediate_files/phylophlan/")
        shell("mkdir -p intermediate_files/phylophlan/input/")
        shell("mkdir -p intermediate_files/phylophlan/output_phylophlan/")
        shell("cp -r intermediate_files/annot/*/*O.faa intermediate_files/phylophlan/input/")
        shell("phylophlan -i {params.in_file} -d {params.database} --diversity low -f {input.config} --nproc {THREADS} --output_folder intermediate_files/phylophlan/output_phylophlan/ --databases_folder src/phylophlan_db")
        shell("mv intermediate_files/phylophlan/output_phylophlan/input_{params.database}/* intermediate_files/phylophlan/output_phylophlan/")
        shell("rm -r intermediate_files/phylophlan/output_phylophlan/input_{params.database}")
       


        
        



        
        
