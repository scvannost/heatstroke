# heatstroke
Code in support of my heat stroke project research under Dr. Doug Lauffenburger i.s.m. Dr. Elizabeth Proctor and USARIEM

## Methods
1. Raw RNAseq reads moved to cluster:\~/data/fastq/ as fastq files.
2. Run cluster.sh to map these reads
3. Run import.sh to move results files to local
4. Run pickle_generation.py to prepare data for further analysis
5. Run deseq_prep.py and deseq.R for DESeq analysis
6. Run log2scaled.py for handmade analysis
7. Run zscored.py to zscore the data, then run relevant sections of log2scaled.py for handmade analysis

### Other files
helpers.py -> a file that makes all the imports needed to run log2scaled.py  
mice.json -> a hand-coded JSON object specifying the experimental groups for each mouse

### Files to download
ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz to cluster:\~/data/  
ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz to cluster:\~/data/  
All genes from http://useast.ensembl.org to biomart_genes.tsv  
All transcripts from http://useast.ensembl.org to biomart_isoforms.tsv  
http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2/h.all.v6.2.symbols.gmt

### Other requirements
github.com/scvannost/randompy/blob/master/DataFrameN.py  
github.com/scvannost/randompy/blob/master/printer.py
