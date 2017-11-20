# Automatic processing of single-cell RNA seq data analysis

This repository contains tools and workflows for the analysis of single cell RNA-Seq data. Currently there are 2 workflows available depending on the sequencing protocol : [smart-seq2](https://www.nature.com/articles/nmeth.2639) and [drop-seq](http://www.cell.com/cell/fulltext/S0092-8674(15)00549-8) (such as [10x genomics](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5241818/)).

In order to deal with such a large amount of data, we will use [Eoulsan](http://www.outils.genomique.biologie.ens.fr/eoulsan/index.html). It is a framework dedicated to high throughput sequencing data analysis on distributed computers. 

## Presentation

## Installation
Please note that core Eoulsan currently runs only on Linux distributions. However, if you use MacOS or Windows, you can install and run Eoulsan through a docker image. For more information about Eoulsan installation, check the dedicated webpage : http://www.outils.genomique.biologie.ens.fr/eoulsan/installing.html 

```{bash, eval=FALSE}
# Installation (Linux only)
wget http://outils.genomique.biologie.ens.fr/eoulsan/eoulsan-2.0-beta5.tar.gz
tar xzf eoulsan-2.0-beta5.tar.gz
cd eoulsan-2.0-beta5

# Installation with docker (Linux / MacOS / Windows)
curl http://outils.genomique.biologie.ens.fr/eoulsan/eoulsan-docker-installer.sh | bash
```

## Quick start

### Create design file

```{bash, eval=FALSE}
genome=/import/kg_csbgn01/genomes/homo_sapiens/hg38/fasta/hg38.fasta
annotation=/import/kg_csbgn01/genomes/homo_sapiens/hg38/annotation/ensembl_Homo_sapiens.GRCh38.84.gtf

eoulsan createdesign /import/kg_csbws03/lehmann/Data/pDC_Soumelis_2017/step1_bcl2fastq/T*R2.fastq $genome $annotation
```

### Customize your analysis with the workflow file

```{bash, eval=FALSE}
workflow=/import/kg_csbws01/lehmann/Eoulsan/pDC-7TP/workflow-10000cells.xml
```

### Run Eoulsan

```{bash, eval=FALSE}
screen
eoulsan.sh -conf eoulsan_conf exec workflow-10000cells.xml design.txt
```

## Help
