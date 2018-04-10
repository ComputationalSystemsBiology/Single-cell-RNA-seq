#!/bin/sh

# RUN BCL2FASTQ (if appropriate)
# more info: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/bcl2fastq-direct
"
bcl2fastq --use-bases-mask=Y26,I8,Y98 \
  --create-fastq-for-index-reads \
  --minimum-trimmed-read-length=8 \
  --mask-short-adapter-reads=8 \
  --ignore-missing-positions \
  --ignore-missing-controls \
  --ignore-missing-filter \
  --ignore-missing-bcls \
  -r 6 -w 6 \
  -R ${FLOWCELL_DIR} \
  --output-dir=${OUTPUT_DIR} \
  --interop-dir=${INTEROP_DIR} \
  --sample-sheet=${SAMPLE_SHEET_PATH}
"

# EXAMPLE DATA (from https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.3.0/hgmm_100)
wget http://cf.10xgenomics.com/samples/cell-exp/1.3.0/hgmm_100/hgmm_100_fastqs.tar;
tar -xf hgmm_100_fastqs.tar;
cat fastqs/hgmm_100_S1_L00?_R1_001.fastq.gz > hgmm_100_R1.fastq.gz;
cat fastqs/hgmm_100_S1_L00?_R2_001.fastq.gz > hgmm_100_R2.fastq.gz;
mkdir data
mv hgmm_100_R?.fastq.gz data/

# INPUTS 
gtf=/import/kg_csbgn01/genomes/homo_sapiens/hg19/annotation/hg19ens91.gtf.bz2
gff=/import/kg_csbgn01/genomes/homo_sapiens/hg19/annotation/hg19ens91.gff.bz2
genomeRef=/import/kg_csbgn01/genomes/homo_sapiens/hg19/fasta/hg19ens91.fasta.bz2

data=../data/*

configuration=../configuration.txt
workflow=../workflow.xml

# CREATE DESIGN
eoulsan createdesign -p $data $genomeRef $gff $gtf

# RUN (start a 'screen' before)
#screen
eoulsan -conf $configuration exec $workflow design.txt

