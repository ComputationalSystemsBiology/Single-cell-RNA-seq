#! /usr/bin/env Rscript

library('ggplot2')

##### Data pre-processing ===============================================================

##### Reading and ordering non auxillary tables------------------------------------------

## Import Expression Matrix
import_exp_matrix<-function(file){
  ### This function import a matrix, and order it by its column names before returning it.
  ### Name of the file or complete path should be given as input.
  
  count<-read.table(file, header=TRUE, 
                    stringsAsFactors=FALSE, sep='\t', row.names=1) 
  count<-count[, order(colnames(count))]
  return(count)
}

##Import Reduced Design File
import_pheno<-function(file){
  ### This function import a table and named its row by its second column before ordering them and returning it.
  ### Name of the file or complete path should be given as input.
  
  pheno<-read.table(file, header=TRUE, 
                    stringsAsFactors=FALSE, sep='\t') 
  # Design v2 case
  i <- grep("SampleId", colnames(pheno))
  # Design v1 case
  if(length(i) == 0) i <-grep("Name", colnames(pheno), fixed=TRUE)
  # Extract sample names
  nms <- pheno[,i]
  # Check whole numeric case
  if( ! any(is.na(as.numeric(nms)))) nms<-paste0("X", nms)
  # Set names
  rownames(pheno)<-nms
  pheno<-pheno[ order(rownames(pheno)),]
  return(pheno)
}


## Import control genes table
import_controls<-function(file){
  controls<-read.table(file, header=TRUE,
                       stringsAsFactors=FALSE, sep="\t")
  return(controls)
}

#### Control genes extraction functions--------------------------------------------------

## Get position of mitochondrial genes in expression matrix row names
get_MT_positions<-function(controls, count){
  ### This function takes as input a control genes table and a count matrix as produced by Eoulsan.
  ### It then seeks features annotated as mitochondrial ("mt_feature") in the count matrix and 
  ### return a vector containing indices of these features.
  
  positions<-c()
  MT_genes<-rownames(controls[grep("mitochondrial",controls$Type),])
  if (length(MT_genes>0)){
    for (gene in 1:length(MT_genes)){
      name<-MT_genes[gene]
      positions<-c(positions, grep(name, rownames(count)))
    }
  }
  return(as.numeric(positions))
}

## Get position of spike-ins genes in expression matrix row names
get_Spike_positions<-function(controls, count){
  ### This function takes as input a control genes table and a count matrix as produced by Eoulsan.
  ### It then seeks features annotated as Spike features ("spike_feature") in the count matrix and 
  ### return a vector containing indices of these features.
  
  positions<-c()
  Spike_genes<-rownames(controls[grep("spike",controls$Type),])
  if (length(Spike_genes>0)){
    for (gene in 1:length(Spike_genes)){
      name<-Spike_genes[gene]
      positions<-c(positions, grep(name, rownames(count)))
    }
  }
  return(as.numeric(positions))
}

#### Processing Expression Matrix -----------------------------------------------------

## Expression Matrix filtering
Filter_exp_matrix<-function(count, Mt_positions, Spike_positions, option='Nuclear'){
  ### This function takes as input a count matrix and positions to be treated.
  ### It then process the matrix according to three options (see below), removing 
  ### unwanted genes and returning the filtered matrix.
  ### Options : -   Total :   Matrix is not filtered at all
  ###           -Endogenous : Matrix is filtered for exogenous genes (spike-ins)
  ###           - Nuclear :   Matrix is filtered for exogenous genes and mitochondrial genes
  
  if(! grepl('Total', option, fixed=TRUE)){
    if(grepl('Nuclear', option, fixed=TRUE)){
      if(length(Mt_positions)>0 || length(Spike_positions)>0){
        count<-count[-(c(Mt_positions, Spike_positions)),]
      }else{
        warning('Mitochondrial positions vector is of length O. \n Nothing was filtered')
      }
    }else{
      if(grepl('Endogenous', option, fixed=TRUE)){
        if(length(Spike_positions)>0){
          count<-count[-c(Spike_positions),]
        }else{
          warning('Spike positions vector is of length O. \n Nothing was filtered')
        }
      }else{
        stop('Invalid option')
      }
    }
  }
  return(count)
}


#### Calcul and filter functions======================================================================

#### Get metrics-----------------------------------------------------------------------------------

## Calculate proportion of mapped reads
Prop_calculator<-function(count, positions){
  ### This function takes as input an expression matrix (cells as columns, features as rows) 
  ### and a vector of row indices. It then calculate the proportion of reads mapped to the 
  ### given positions.
  
  total<-colSums(count)
  if(length(positions)>0){
    controls<-colSums(count[positions,])
    prop<-controls/total
    return(prop)
  }
  warning(paste0('Positions vector is of length 0.', '\n No feature found.'))
}

## Calculate number of detected genes
get_detected_genes<-function(count, detection){
  nb_genes<-c()
  for (cell in 1:ncol(count)){
    nb_genes<-c(nb_genes, sum(count[,cell]>detection))
  }
  return(nb_genes)
}
#### Filters------------------------------------------------------------------------------

## Filter cell by number of failure
filter_cells<-function(cells, threshold){
  ### This function takes as input a vector with repeated element and an integer value for threshold.
  ### It then detects and returns elements appearing more than the given threshold.
  
  cells<-data.frame(table(cells))
  filt<-which(cells$Freq>=threshold)
  return(as.numeric(levels(cells[,1])[filt]))
}

#### Main function===============================================================================

main<-function(file1, file2, file3,
               detection=10, exp_option='Nuclear',
               exp_threshold=4000, nbr_threshold=500000,
               propmt_threshold=0.2,
               propsp_threshold=0.5,
               nb_filters=1,
               output1, output2, output3){
  ### This function integrates all parameters given by user and realises the quality 
  ### checking of the data.
  ### Set threshold to 0 to disable filtering
  
  # Importing Data-------------------------------------------  
  count<-import_exp_matrix(file1)
  pheno<-import_pheno(file2)
  controls<-import_controls(file3)
  
  # Initializing filtering result vectors--------------------
  exp_cells<-c()
  nbr_cells<-c()
  mt_cells<-c()
  sp_cells<-c()
  
  # Extracting controls positions ---------------------------
  MT_positions<-get_MT_positions(controls, count)
  Spike_positions<-get_Spike_positions(controls, count)
  rm(controls)
  
  # Treating data for basic metrics -------------------------
  fmatrix<-Filter_exp_matrix(count, MT_positions, Spike_positions, option=exp_option)
  exp_features<-get_detected_genes(fmatrix, detection)
  nb_reads<-colSums(fmatrix)
  rm(fmatrix)
  if(exp_threshold){
    exp_cells<-which(exp_features < exp_threshold)
  }
  if(nbr_threshold){
    nbr_cells<-which(nb_reads< nbr_threshold)
  }
  pheno$Nb_features<-exp_features
  rm(exp_features)
  pheno$Nb_reads<-nb_reads
  rm(nb_reads)
  
  # Code from previous version allowing different options for feature and reads counting------  
  # fmatrix<-Filter_exp_matrix(count, MT_positions, Spike_positions, option=exp_option)
  # exp_features<-get_detected_genes(fmatrix, detection)
  # rm(fmatrix)
  # if(exp_threshold){
  #   exp_cells<-which(exp_features < exp_threshold)
  # }
  # pheno$Nb_features<-exp_features
  # rm(exp_features)
  # fmatrix<-Filter_exp_matrix(count, MT_positions, Spike_positions, option=nbr_option)
  # nb_reads<-colSums(fmatrix)
  # rm(fmatrix)
  # if(nb_reads){
  #   nbr_cells<-which(nb_reads< nbr_threshold)
  # }
  # pheno$Nb_reads<-nb_reads
  # rm(nb_reads)
  
  # Treating data for advanced metrics ----------------------
  if(propmt_threshold){
    prop<-Prop_calculator(count, MT_positions)
    pheno$Prop_MT<-prop
    mt_cells<-which(prop > propmt_threshold)
  }
  
  if(propsp_threshold){
    prop<-Prop_calculator(count, Spike_positions)
    pheno$Prop_Sp<-prop
    sp_cells<-which(prop > propsp_threshold)
  }
  
  # Export results ------------------------------------------
  write.table(pheno, output1, row.names=TRUE, sep="\t", col.names=TRUE)
  
  # Extract filtered cells -----------------------------------
  filter<-filter_cells(c(exp_cells, nbr_cells, mt_cells, sp_cells), nb_filters)
  
  
  # Plot Raw data with filtered cells annotated -------------- 
  pdf(paste0('Raw_Cellplot.pdf')) ##### add function to get experiment from input files
  plot(pheno$Nb_reads, pheno$Nb_features, xlab='Number of reads', ylab='Number of detected features', 
       xlim=c(min(min(pheno$Nb_reads), nbr_threshold),max(pheno$Nb_reads)+1000),
       ylim=c(min(min(pheno$Nb_features),exp_threshold), max(pheno$Nb_features)),
       pch=19, col=alpha("gray", 0.5))
  abline(h=exp_threshold, col='red', lty='dashed')
  abline(v=nbr_threshold, col='red', lty='dashed')
  points(x=pheno$Nb_reads[filter], y=pheno$Nb_features[filter], col='red', pch=19)
  dev.off()
  
  # Filtered and export table and matrices--------------------
  if (length(filter) > 0) {
    count<-count[,-filter]
    pheno<-pheno[-filter,]
  }
  write.table(count, output2, row.names=TRUE, sep="\t", col.names=TRUE)
  write.table(pheno, output3, row.names=TRUE, sep="\t", col.names=TRUE)
  
  # Export filtered dot plot ---------------------------------
  pdf(paste0( 'Filtered_cellplot.pdf')) ##### add function to get experiment from input files
  plot(pheno$Nb_reads, pheno$Nb_features, xlab='Number of reads', ylab='Number of detected features',
       pch=19, col=alpha("gray", 0.5))
  dev.off()
}

args<-commandArgs(TRUE)

main(args[1], args[2], args[3],  
     detection=as.numeric(args[4]),
     exp_option = args[5],
     exp_threshold = as.numeric(args[6]),
     nbr_threshold = as.numeric(args[7]),
     propmt_threshold = as.numeric(args[8]), 
     propsp_threshold = as.numeric(args[9]), 
     nb_filters = as.numeric(args[10]),
     args[11], args[12], args[13])
