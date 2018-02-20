library(Seurat)

##### Data pre-processing ===============================================================

#### Testing parameters

as.valid.logical<-function(test) {
  res<-as.logical(test)
  if(is.na(res)) stop("Logical parameters must be : one of true, TRUE, True, False, false, or FALSE")
  return(res)
}


##### Reading and ordering non auxillary tables------------------------------------------

## Import Expression Matrix
import_exp_matrix<-function(file){
  ### This function import a matrix, and order it by its column names before returning it.
  ### Name of the file or complete path should be given as input.
  
  count<-read.table(file, header=TRUE, 
                    stringsAsFactors=FALSE, sep="\t", row.names=1) 
  count<-count[, order(colnames(count))]
  return(count)
}

##Import Reduced Design File
import_pheno<-function(file){
  ### This function import a table and named its row by its second column before ordering them and returning it.
  ### Name of the file or complete path should be given as input.
  
  pheno<-read.table(file, header=TRUE, 
                    stringsAsFactors=FALSE, sep="\t") 
  rownames(pheno)<-pheno[,2]
  pheno<-pheno[ order(rownames(pheno)),]
  return(pheno)
}

## Import control genes table
import_controls<-function(file){
  controls<-read.table(file, header=TRUE,
                       stringsAsFactors=FALSE, sep="\t")
  return(controls)
}

#### Apply corrections =================================================

normalize.seurat<-function(seurat.object, file) {
  pheno<-import_pheno(file)
  seurat.object@data<-log(t(t(seura.object@raw.data)/pheno$Size_factors) + 1)
}

#### HVG detection =====================================================

### Return highly variable genes identify by scran
hvg.with.scran<-function(seurat.object, genes, spike = TRUE, 
                         low.mean = 0.01, high.mean = 5, bio = 1 ) {
  library(scran)
  # construct object
  SCE <- newSCESet(exprsData = seurat.object@data, logged = TRUE)
  
  # find spikes if needed (can be long)
  if(spike && any(grepl("spike", genes$Type))) {
    spikes <- rownames(genes)[grep("spike", genes$Type)]
    isSpike(SCE) <- sapply(spike, grep, x = rownames(SCE))
  }
  
  # find HVGs
  fit <- trendVar(SCE, trend = "loess",  use.spikes = (! is.null(isSpike(SCE))) )
  dc <- decomposeVar(SCE, fit)
  
  # Plot results
  pdf("scranFittedVariability.pdf")
  plot(dc$mean, dc$total, xlab="Mean log-expression", 
         ylab="Variance", pch ="19", col = "gray")
  o <- order(dc$mean)
  lines(dc$mean[o], dc$tech[o], col="blue", lwd=2)
  if (! is.null(isSpike(SCE))) points(fit$mean, fit$var, col="blue", pch=19)
  dev.off()
  
  # get and return hvg
  hvg<-rownames(dc)[which( dc$mean > low.mean && dc$mean <= high.mean && dc$bio >= bio)]
  return(hvg)
}


#### Cluster quality metrics ===========================================

#### Silhouette plot ---------------------------------------------------

### Calculate and plot silhouette from seurat object with clusters
### parameters :
###   seurat.object : seurat object to use
###   pc.used : principal component used for clustering
###   distance : distance to use for silhouette calculation
###   plot : weither to plot silhouette or not
### returns : silhouette object from cluster package

silhouette.seurat<-function(seurat.object, pc.used, distance="euclidean", plot=TRUE) {
  library(cluster)
  coord <- seurat.object@pca.rot
  coord <- coord[,pc.used]
  clusters <- seurat.object@ident
  d <- dist(coord, method=distance)
  sil<-silhouette(as.numeric(clusters), dist=d)
  if(plot){
    plot(sil, col=as.factor(clusters[order(clusters, decreasing=FALSE)]), main="Silhouette plot of Seurat clustering", lty=2)
    abline(v=mean(sil[,3]), col="red4", lty=2)
  }
  return(sil)
}


#### Main function 
main<-function(matrix.file, cell.file, genes.file, # input files
               normalize = TRUE, scale = TRUE, # normalization and scaling
               n.cells = 0, detection = 10, # filtering genes on detection
               hvg.detection = "none", low.mean, high.mean, var.threshold, spike = FALSE, # filtering on variability
               jackstraw.replicate = 100, jackstraw.prop = 0.1, sig.threshold = 0.05, score.threshold = 1e-5, # selecting significative dimensions
               resolution = 0.6, k.param = 30, k.scale = 25, algorithm="Louvain", sparse = FALSE, # clustering options
               genesOutput, cellsOutput) {
  
  library(methods)
  
  # Checking parameters
  if( ! hvg.detection %in% c("scran", "seurat", "none")) stop("hvg.detection must be one of : scran, seurat, or none")
  if( ! algorithm %in% c("Louvain", "Louvain.multilevel", "SLM")) stop("algorithm must be one of : Louvain, Louvain.multilevel, or SLM")
  
  ## import data and construct Seurat object
  count <- import_exp_matrix(matrix.file)
  object <- new("seurat", raw.data=count)
  
  ## Normalize data
  if(normalize) object <- normalize.seurat(object, cell.file)
  else object@data <- log( object@raw.data + 1 )
 
  ## Filter genes on expression
  if(n.cells) {
    kept <- which(rowSums( object@raw.data >= detection) >=  n.cells )
    object@data <- object@data[kept,]
  } 
  
  ## Scale data
  if(scale) object@scale.data <- t(scale(t(object@data)))
  else object@scale.data <- object@data
 
  ## Select HVG
  genes<-import_controls(genes.file)
  if( hvg.detection == "none") hvg <- rownames(object@scale.data)
  if( hvg.detection == "scran") hvg <- hvg.with.scran(object, genes, spike, low.mean, 
                                                      high.mean, bio = var.threshold ) 
  if( hvg.detection == "seurat") {
    pdf("seuratFittedVariation.pdf")
    object <- MeanVarPlot(object, fxn.x = expMean, fxn.y = logVarDivMean, 
                          x.low.cutoff = low.mean, x.high.cutoff = high.mean,
                          y.cutoff = var.threshold, do.contour = FALSE)
    dev.off()
    hvg <- object@var.genes
  } 
    ## If no genes are kept keep all genes
  if(length(hvg) == 0) {
    hvg <- rownames(object@scale.data)
    warning("No Highly variable gene found. Continue with all genes.")
  }
  
  ## Save genes used for dimensionnality reduction
  genes$Used.for.analysis <- rep(FALSE, nrow(genes))
  genes$Used.for.analysis[sapply(hvg, grep, x = rownames(genes))] <- TRUE
  
  write.table(genes, genesOutput, col.names = TRUE, row.names= TRUE, sep="\t")
  
  ## Perform PCA and JackStraw
  object@cell.names <- colnames(object@scale.data)
  object <- PCA(object, pc.genes = hvg, do.print = FALSE, 
                pcs.store = min(ncol(object@scale.data), length(hvg)))
  object <- JackStraw(object, num.pc = dim(object@pca.rot)[2], prop.freq = jackstraw.prop,
                      num.replicate = jackstraw.replicate, do.print = FALSE)
  
  ## Select significative PCs
  
  pc.scores <- object@jackStraw.empP
  p.value <- 0
  i <- 1
  pcs<-NULL
  while( p.value <= sig.threshold && i < ncol(pc.scores)) {
    n.obs <- sum(pc.scores[,i] <= score.threshold)
    if(n.obs == 0 ) {
      p.value <- 1
      break
    }
    tot <- nrow(pc.scores)
    n.theo <- tot * score.threshold
    p.value <- prop.test(c(n.obs, n.theo), c(tot, tot))$p.val
    if(p.value <= sig.threshold) pcs <- c(pcs, i)
    i <- i + 1
  }
  rm(pc.scores)
  
  # Treat few pcs 
  if(length(pcs) < 2) pcs <- 1:2
  
  ## plot and save results
  p<-JackStrawPlot(object, PCs = c(pcs, length(pcs)+1), score.thresh = score.threshold)
  ggplot2::ggsave("PCSignificativityPlots.pdf", p, device = "pdf")
  
  ## Produce PCs heatmaps
  ncells<-min(200, ncol(object@scale.data))
  pdf("PCsHeatmap.pdf")
  PCHeatmap(object, pc.use = pcs, cells.use = ncells, 
            do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
  dev.off()
  
  ## Cluster cells
  if(algorithm == "Louvain") algorithm <- 1
  if(algorithm == "Louvain.multilevel") algorithm <- 2
  if(algorithm == "SLM") algorithm <- 3
  
  object@data.info <- data.frame(cells.names = object@cell.names)
  object <- FindClusters(object, pc.use = pcs, resolution = resolution,
                         print.output = FALSE, save.SNN = TRUE,  k.param = k.param, 
                         k.scale = k.scale, algorithm = algorithm, do.sparse = sparse)
  
  ## Visualize using TSNE
  object <- RunTSNE(object, dims.use = pcs, do.fast=FALSE)
  p <- TSNEPlot(object, do.return = TRUE)
  ggplot2::ggsave("clusterPlot.pdf", plot = p, device = "pdf")
  
  ## Evaluate with silhouette
  pdf("silhouettePlot.pdf")
  sil<-silhouette.seurat(object, pc.used = pcs, distance="euclidean", plot=TRUE)
  dev.off()
  
  ## Save clustering
  pheno <- import_pheno(cell.file)
  pheno$Clusters <- object@ident
  write.table(pheno, cellsOutput, row.names = TRUE, col.names = TRUE, sep = "\t")
  saveRDS(object, file = "seurat.RDS")
}

## Treat Parameters
args <- commandArgs(TRUE)
arg4 <- as.valid.logical(args[4])
arg5 <- as.valid.logical(args[5])
arg12 <- as.valid.logical(args[12])
arg21 <- as.valid.logical(args[21])

## Launch main function
main(args[1], args[2], args[3], # input files
     arg4, arg5, # correction options
     as.numeric(args[6]), as.numeric(args[7]), # filtering genes on detection
     args[8], as.numeric(args[9]), as.numeric(args[10]), as.numeric(args[11]), arg12, #filtering on variability
     as.numeric(args[13]), as.numeric(args[14]), as.numeric(args[15]), as.numeric(args[16]), # jackstraw options
     as.numeric(args[17]), as.numeric(args[18]), as.numeric(args[19]), args[20], arg21, # clustering options
     args[22], args[23])
