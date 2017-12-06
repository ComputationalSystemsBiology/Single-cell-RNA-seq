# README

## Repository 

This repository contains two R functions simulating counts from single cell sequencinge experiment.


## LunSim.R

This function is based on "Pooling across cells to normalize single-cell RNA sequencing data with many zero counts" by A.T Lun et al.
It allows random sampling of counts for a given number of genes and a given number of populations.

## DropoutSim.R

This function extends previous function by adding sampling noise and gene length.
To simulate Smart-seq data you should activate all options.
To simulate Chromium data you shoul desactivate length, and amplification bias options.
The function can use distinct dispersion value to simulate population, you can use dispersion.tsv to have access to dispersion estimated on real dataset.
