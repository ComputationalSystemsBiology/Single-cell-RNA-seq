# Automatic processing of single-cell RNA seq data analysis

This repository contains tools and workflows for the analysis of single cell RNA-Seq data. Currently there are 2 workflows available depending on the sequencing protocol : [smart-seq2](https://www.nature.com/articles/nmeth.2639) and [drop-seq](http://www.cell.com/cell/fulltext/S0092-8674(15)00549-8) (such as [10x genomics](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5241818/)).

In order to deal with such a large amount of data, we will use [Eoulsan](http://www.outils.genomique.biologie.ens.fr/eoulsan/index.html). It is a framework dedicated to high throughput sequencing data analysis on distributed computers. 

## Presentation

## Installation
Please note that core Eoulsan currently runs only on Linux distributions. However, if you use MacOS or Windows, you can install and run Eoulsan through a docker image. For more information about Eoulsan installation, check the dedicated webpage : http://www.outils.genomique.biologie.ens.fr/eoulsan/installing.html 

```{bash, eval=FALSE}
# Installation (Linux only)
wget http://outils.genomique.biologie.ens.fr/eoulsan/eoulsan-2.0.tar.gz
tar xzf eoulsan-2.0.tar.gz
cd eoulsan-2.0

# Installation with docker (Linux / MacOS / Windows)
curl http://outils.genomique.biologie.ens.fr/eoulsan/eoulsan-docker-installer.sh | bash
```

## Quick start

### Create design file

```{bash, eval=FALSE}
genome=/path/to/fasta
gtf=/path/to/gtf
gff=/path/to/gff
data=/path/to/fastq

eoulsan createdesign $data $genome $gtf $gff

# if paired-end reads (e.g. for 10xGenomics)
eoulsan createdesign -p $data $genome $gtf $gff
```

### Customize your analysis with the workflow file

### Run Eoulsan

## Help
