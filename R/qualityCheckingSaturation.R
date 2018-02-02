#! /usr/bin/env Rscript

library('scales')
library('parallel')

#### Data pre-processing ===============================================================

#### Reading and ordering non auxillary tables------------------------------------------

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
                       stringsAsFactors=FALSE, sep='\t')
  return(controls)
}


#### Control genes extraction functions--------------------------------------------------

## Get position of mitochondrial genes in expression matrix row names
get_MT_positions<-function(controls, count){
  ### This function takes as input a control genes table and a count matrix as produced by Eoulsan.
  ### It then seeks features annotated as mitochondrial ('mt_feature') in the count matrix and 
  ### return a vector containing indices of these features.
  
  positions<-c()
  MT_genes<-rownames(controls[grep('mitochondrial',controls$Type),])
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
  ### It then seeks features annotated as Spike features ('spike_feature') in the count matrix and 
  ### return a vector containing indices of these features.
  
  positions<-c()
  Spike_genes<-rownames(controls[grep('spike',controls$Type),])
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


#### Saturation functions ========================================================================
## Resample function
resample<-function(vector, iterations=10, threshold=10, round=10000, nbPoints=20) {
  ### This function takes as input a vector of count for various element (one element in 
  ### the vector = count for one element)
  ### Options are: 
  ###  iterations:  number of iteration for each step, default value is 10
  ###  threshold:   number of event for an element to be considered as detected, 
  ###               default is 10
  ###  round :      Magnitude for down rounding counts, default is 10 000
  ###  nbPoints :   Number of points to sample
  ### returns : sampling results in a dataframe
  
  # getting down rounded max depth of resampling
  tot_count<-sum(vector)
  tot_count<-floor(tot_count/round)*round
  
  # Creating vector for random sampling
  data<-rep(1:length(vector), vector)
  
  # Set increment
  steps=tot_count/nbPoints
    
  # Create sampling vector
  depths<-seq(0, tot_count, steps)
  
  # Internal function
  ## Draw and count (sampling)
  drawNcount<-function(index, data, depth, threshold){
    samp<-data.frame(table(sample(data, depth, replace=FALSE)))
    count<-sum(samp$Freq>=threshold)
    names(count)<-paste0('draw', index)
    return(count)
  }
  ## Resume sampling
  resume<-function(depth, data, threshold, iterations){
    sampl<-sapply(1:iterations, drawNcount, data=data, 
                  depth=depth, threshold=threshold)
    return(c(depth, mean(sampl), sd(sampl)))
  }
  
  # Resample
  result<-sapply(depths, resume, data=data, threshold=threshold, iterations=10)
  
  #return result
  rownames(result)<-c('depths','means','sd')
  return(data.frame(t(result)))
}


## Fitting function
fitMMModel<-function(table){
  ### This function fits a Michaelis Menten model using a table of resampled features
  ### table :   table of resampled features
  ### return :  a List : max = maximum of the Model,  half = X value to have 1/2 of max
  ###           min = minimum of the fitting function
  
  # function to minimize for non linear model
  michaelis<-function(param, data){
    max<-param[1]
    half<-param[2]
    diff<-sum(abs(data$means-(max*data$depths/(data$depths+half))))
    return(diff)
  }
  # Fit a linear model using Hanes-Woolf equation to obtain near optimal parameters value
  HW<-lm(means~depths, data.frame(depths=table[-which(table[,2]==0),1], 
                                  means=table[-which(table[,2]==0),1]/table[-which(table[,2]==0),2]))
  
  # Extract parameters value
  HW.max<-1/HW$coefficients[2]
  HW.half<-HW$coefficients[1]/HW$coefficients[2]
  
  # Fit non linear model using near optimal starting values
  results<-nlm(michaelis, c(HW.max, HW.half), table)
  return(list(max=results$estimate[1], half=results$estimate[2], min=results$minimum) )
}


## Fucntion for fitting evaluation
testFit<-function(measured, predicted, name){
  ### Fit and plot a linear model such as predicted = f(measured)
  ### measured :  Measured values (numeric vector)
  ### predicted : values predicted by model (numeric vector)
  ### return : R² of linear fit (indicating goodness of fit)
  
  # Fit linear model
  model<-lm(predicted~measured)
  # Get R²
  r<-summary(model)$r.squared
  
  # Plot and save 
  pdf(paste0(name, 'FittingPlot.pdf'))
  plot(measured, predicted, xlab='Measured values', ylab='Predicted values', 
       main=paste0(name, ' goodness of fit'), pch=19)
  abline(model, col='chocolate3')
  dev.off()
  
  #Return R²
  return(r)
}


## Plotting function
plotSaturation<-function(table, predicted, name) {
  ### This function plots and saves saturation resampling given in a dataframe
  ### Parameters : 
  ###     table : table gathering data to plot with sampling on first column means 
  ###             on second column and sd on third column
  ###     name : name of the sample (for plot title)
  
  sampl<-table[,1]
  mean<-table[,2]
  sd<-table[,3]
  
  # data
  plot<-ggplot(data.frame(reads=sampl, genes=mean), 
               aes(x=reads, y=genes))
    # Points and line
    plot<-plot+geom_line(colour='royalblue3')
    plot<-plot+geom_point(colour='royalblue1',size=2)
    # Error bars
    plot<-plot+geom_errorbar(aes(ymin=mean-sd, 
                      ymax=mean+sd), colour='royalblue3', width=(sampl[2]-sampl[1])/5)
    # Adding model curve
    plot<-plot+geom_line(aes(x=sampl, y=predicted), colour='chocolate3')
    # esthetics
    plot<-plot+theme_bw()
    plot<-plot+xlab('Number of reads')+ylab('Number of detected features')
    plot<-plot+ggtitle(paste0(name,' Saturation Curve '))
  ggsave(paste0('Saturation_curve_', name, '.pdf'), plot=plot, device='pdf')
}


## Function to treat columns from modified counting matrix
saturationCol<-function(col, iterations=10, threshold=10, round=10000, nPoints=20) {
  ###This function takes as input an expressionvector and gives as output a
  ###saturation plot, fitting evaluation plot, and a vector containing results 
  ###of the fitting. Requires ggplot2.
  ###Options are: 
  ### iterations: number of iteration for each step, default value is 10
  ### threshold:  number of event for a gene to be considered as detected
  ###             (set as a parameter for the module, default is 10)
  ### round:      Magnitude for down rounding counts data, default is 
  ###             10 000.
  ### nPoints :   number of points to sample for saturation
  ### Return :    a vector containing : model coefficients (maximum and X  
  ###             value for Y=1/2*maximum), fitting value, saturation estimate
  
  # Extract count and name
  l<-length(col)
  name<-col[l]
  exp<-as.numeric(col[-l])
  
  # Resample data
  rs<-resample(exp, threshold=threshold, iterations=iterations, round=round, nbPoints=nPoints)
  # Check sampling, if cell has not enough reads return NA values
  if(nrow(rs)< nPoints) return(rep(NA, 4))
  
  #fit model on resampling
  MM<-fitMMModel(rs)
  # Get predicted values
  predicted<-MM$max*rs$depths/(MM$half+rs$depths)
  
  # Plotdata
  if(! dir.exists('SaturationPlots')) dir.create('SaturationPlots')
  setwd('SaturationPlots')
  plotSaturation(rs, predicted, name)
  setwd('..')
  
  # Test Goodness of fit
  if(! dir.exists('FitEvaluation')) dir.create('FitEvaluation')
  setwd('FitEvaluation')
  r<-testFit(rs$means, predicted, name)
  setwd('..')
  
  # Evaluate saturation
  sat<-sum(exp >= threshold)/MM$max
  
  # Return results
  return(c(MM$max, MM$half, sat, r))
}

## Function for index correction
correctIndex<-function(index, ref){
  ### given an index and removed indices, correct the index
  ### index : a numeric value
  ### ref :   a numeric vector containing removed values (increasing)
  ### return corrected index value
  i<-1
  while(ref[i] <= index & i <=length(ref)){
    i<-i+1
    index<-index+1
  }
  return(index)
}

## Function identifying outliers (won't work if more than 5000 values are passed)
getOutliers<-function(values, normThreshold=0.001){
  ### This function identify long left tail outliers
  ### values : numerical vector
  ### normThreshold : probability under gaussian expectation to reach
  ### return : vector containing indices of removed values
  
  # Initialize variables
  p<-shapiro.test(values)$p.value
  proba<-p
  rmv<-NULL
  
  # Remove outliers
  while(p<normThreshold) {
    if( ! is.null(rmv)) {
      mini<-which.min(values[-rmv])
      mini<-correctIndex(mini, rmv)
      rmv<-c(rmv, mini)
      rmv<-rmv[order(rmv, decreasing=FALSE)]
    }
    else {
      mini<-which.min(values)
      rmv<-mini
    }
    if( max(values[-rmv]) != min(values[-rmv]) & length(rmv)<length(values)-3) {
      p<-shapiro.test(values[-rmv])$p.value
      proba<-c(proba, p)
    } else {
      warning("Outliers removal failed")
      break
    }
  }
  
  
  # Plot P values
  pdf('gaussianProbability.pdf')
  plot(1:length(proba)-1, proba, type='b', pch=19, 
       xlab='Number of removal', ylab='P value (Shapiro)')
  dev.off()
  
  # return indices
  return(rmv)
}


## Saturation global function
saturation <- function (count, satThreshold=0.7, fitThreshold=0.97, iterations = 10, 
                        threshold=10, round=10000, nPoints=20, nCores=2) {
  ###This function takes as input an expression matrixand gives as output a
  ###saturation plot and a fitting plot for heach column of the matrix (Requires ggplot2)
  ###
  ###Options are: 
  ### iterations: number of iteration for each step, default value is 10
  ### threshold:  number of event for a gene to be considered as detected
  ###             (set as a parameter for the module, default is 10)
  ### round:      Magnitude for down rounding counts data, default is 
  ###             10 000.
  ### nPoints :   number of points to sample for saturation
  ### nCores :    number of threads to run
  ### Return :    a vector containing : model coefficients (maximum and X  
  ###             value for Y=1/2*maximum), fitting value, saturation estimate
  
  # Add a line containing Sample name to input matrix
  mat<-rbind(count, colnames(count))
  
  # Prepare parallelized computing
  clus<-makeCluster(nCores)
  clusterEvalQ(clus, library('ggplot2'))
  clusterExport(clus, c('iterations', 'threshold', 'round', 'nPoints'), envir=environment())
  clusterExport(clus, c( 'resample', 'fitMMModel', 'plotSaturation', 'testFit', 'saturationCol'))
  
  # Following operations are applied to each cell
  metrics<-parCapply(clus, mat, saturationCol, iterations=iterations, threshold=threshold,
                     round=round, nPoints=nPoints )
  stopCluster(clus)
  metrics<-matrix(metrics, nrow=4, byrow=FALSE)
  colnames(metrics)<-colnames(count)
  rownames(metrics)<-c('max', 'half', 'saturation', 'r.squared')
  
  # Change Working Directory to save plots
  setwd('FitEvaluation')
  
  # Get cells to suppress indices
  i<-which(apply(metrics[1:2,], 2, function(x){any(x <0, na.rm=TRUE)}) | metrics[3,] < satThreshold | apply(metrics, 2, anyNA))
  outliers<-which(metrics[4,] < fitThreshold) 
  i<-c(i, outliers)
  
  # Plot data before and after removal
  rmv<-rep('grey', ncol(metrics))
  rmv[i]<-'red'
  
    ## Full set
  pdf('qualityAll.pdf')
  plot(t(metrics[c(1,4),]), xlab='Maximum', ylab='R squared', pch=19, 
       col=rmv, main='Saturation quality \n (all cells)')
  dev.off()
  
  pdf('densityAll.pdf')
  plot(density(metrics[4,], na.rm=TRUE), main="R squared density \n (all cells)")
  dev.off()
  
    ## Filtered cell
  if(length(i) > 0) {
    pdf('qualityFiltered.pdf')
    plot(t(metrics[c(1,4),-i]), xlab='Maximum', ylab='R squared', pch=19, main='Saturation quality \n (good cells)')
    dev.off()
  
    pdf('densityFiltered.pdf')
    plot(density(metrics[4, -i], na.rm=TRUE), main="R squared density \n (good cells)")
    dev.off()
  }
  
  # Back to main directory
  setwd('..')
  
  # return indices  
  return(list(indices=i, metrics=t(metrics)))
}


#### Main function===============================================================================
main<-function(file1, file2, file3, detection=10,  exp_option='Nuclear',
               satThreshold=0.7, fitThreshold=0.01, round=10000, nCores=2,
               propmt_threshold=0.2, propsp_threshold=0.5, nb_filters=1,
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
 
  # Filtering Metrix ----------------------------------------
  fmatrix<-Filter_exp_matrix(count, MT_positions, Spike_positions, option=exp_option)
  
  # Finding unsaturated cells -------------------------------
  unsat<-saturation(fmatrix, threshold = detection, satThreshold = satThreshold, 
                      fitThreshold = fitThreshold, round = round, nCores= nCores)
  
  # Extracting basic metrics --------------------------------
  exp_features<-get_detected_genes(fmatrix, detection)
  nb_reads<-colSums(fmatrix)
  
  # Saving in phenotype table -------------------------------
  pheno$Nb_features<-exp_features
  rm(exp_features)
  pheno$Nb_reads<-nb_reads
  rm(nb_reads)

  
  # Treating data for advanced metrics ----------------------
  if(propmt_threshold & length(MT_positions) > 0){
    prop<-Prop_calculator(count, MT_positions)
    pheno$Prop_MT<-prop
    mt_cells<-which(prop > propmt_threshold)
  }
  
  if(propsp_threshold & length(Spike_positions) > 0){
    prop<-Prop_calculator(count, Spike_positions)
    pheno$Prop_Sp<-prop
    sp_cells<-which(prop > propsp_threshold)
  }
  
  
  pheno<-cbind(pheno, unsat$metrics)
  
  # Export results ------------------------------------------
  write.table(pheno, output1, row.names=TRUE, sep='\t', col.names=TRUE)
  
  # Extract filtered cells -----------------------------------
  filter<-filter_cells(c(unsat$indices, mt_cells, sp_cells), nb_filters)
 
  
  # Plot Raw data with filtered cells annotated --------------
  colour<-rep(alpha('gray', 0.5), nrow(pheno))
  colour[filter]<-'red'
  pdf(paste0('Raw_Cellplot.pdf'))
  plot(pheno$Nb_reads, pheno$Nb_features, xlab='Number of reads', ylab='Number of detected features', 
         xlim=c(min(pheno$Nb_reads)-1000,max(pheno$Nb_reads)+1000),
         ylim=c(min(pheno$Nb_features)-50, max(pheno$Nb_features)+50),
         pch=19, col=colour)
  dev.off()
  
  # Filtered and export table and matrices--------------------
  if (length(filter) > 0) {
    count<-count[,-filter]
    pheno<-pheno[-filter,]
  }
  write.table(count, output2, row.names=TRUE, sep='\t', col.names=TRUE)
  write.table(pheno, output3, row.names=TRUE, sep='\t', col.names=TRUE)
  
  # Export filtered dot plot ---------------------------------
  pdf(paste0( 'Filtered_cellplot.pdf')) 
  plot(pheno$Nb_reads, pheno$Nb_features, xlab='Number of reads', ylab='Number of detected features',
         pch=19, col=alpha('gray', 0.5))
  dev.off()
}

#### Running script==============================================================================
args<-commandArgs(TRUE)

main(args[1], args[2], args[3], 
     detection=as.numeric(args[4]), exp_option = args[5],
     satThreshold = as.numeric(args[6]), fitThreshold = as.numeric(args[7]),
     round = as.numeric(args[8]), nCores = as.numeric(args[9]),
     propmt_threshold = as.numeric(args[10]), 
     propsp_threshold = as.numeric(args[11]), 
     nb_filters = as.numeric(args[12]),
     args[13], args[14], args[15])
