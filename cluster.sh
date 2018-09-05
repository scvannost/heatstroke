#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scv@mit.edu
#######

module load r/3.3.1
module load perl
module load rsem
module load star

# gtf: ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz
# fasta: ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
rsem-prepare-reference --gtf ~/data/Mus_musculus.GRCm38.93.gtf \
    --star \
    -p 4 \
    ~/data/Mus_musculus.GRCm38.dna.toplevel.fa \
    ~/data/ref/mouse

for i in ~/data/fastq/*; do
mkdir ~/data/mapped/${i:25:2}
rsem-calculate-expression --star \
    -p 16 \
    --output-genome-bam \
    $i \
    ~/data/ref/mouse \
    ~/data/mapped/${i:25:2}/${i:25:2}
done
