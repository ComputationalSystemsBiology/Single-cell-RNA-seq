### simulateCountLun : this function reproduces count simulation as in "Pooling across cells 
###    to normalize single-cell RNA sequencing data with many zero counts" by A.T Lun et al.
### parameters :
###    - nbgenes : total number of simulated genes
###    - nbDEgenes : number of differencially expressed genes
###    - nbpop : number of population defined by DE genes
###    - nbcells : number of cells by population
###    - pup : probability for a gene to be upregulated
###    - fcup : fold change for upregulated genes
###    - fcdown : fold change for downregulated genes
###    - cellBiases : weither to add or not cell specific bias


simulateCountLun<-function(nbgenes=10000, nbDEgenes = 3000 ,nbpop = 3, 
                           nbcells = 100, pup = 0.4, fcup = 5, fcdown = 0,
                           cellBiases=TRUE) {
  
  # Check parameters
  if(fcdown > 1) stop("Fold-change for down regulated genes must be less than 1.")
  if(nbDEgenes > nbgenes/nbpop) stop("Too many differencially expressed genes in each population")
  if(fcup < 1) stop("Fold change for upregulated genes must be more than 1.")
  if(pup > 1) stop("More than 100% of DE genes upregulated.")
  
  # Initialize return variables
  expMatrix<-c()
  biases<-c()
  
  # Creating gene caracteristics
  expVector<-rgamma(nbgenes, shape=2, rate=2)
 
  
  # Select differentially expressed genes
  DEvector<-rep(0, nbgenes)
  genesInd<-c(1:nbgenes)
  for (p in 1: nbpop) {
    DEgenesInd<-sample(genesInd, nbDEgenes, replace=FALSE)
    genesInd<-genesInd[-DEgenesInd]
    DEvector[DEgenesInd]<-p
  }
  
  
  cellsName<-c()
  expMatrix<-seq(from=1, to=nbgenes)
  # Adding population caracteristics
  for(p in 1: nbpop ) {
    popExpCount<-expVector
    for(i in 1:nbgenes) {
      foldChange<-1
      if(DEvector[i] == p) {
        if(rbinom(1, 1, pup) == 1) {
          foldChange<-fcup
          DEvector[i]<-paste0("Up in Pop", p)
        } else {
          foldChange<-fcdown
          DEvector[i]<-paste0("Down in Pop", p)
        }
        popExpCount[i]<-expVector[i]*foldChange
      }
    }
    
    # Adding cells caracteristics
    for(c in 1:nbcells) {
      cellBias<-1
      cellsName<-c(cellsName , paste0("Pop", p, "-cell", c))
      if(cellBiases) {
        cellBias<-2^(rnorm(1, 0, 0.25))
      }
      biases<-c(biases, cellBias)
      cellCount<-sapply(popExpCount, function(x) {rnbinom(1, mu=cellBiases*x, size=0.1) })
      expMatrix<-cbind(expMatrix, cellCount) 
    }
  }

  # Process Matrix
  rownames(expMatrix)<-expMatrix[,1]
  expMatrix<-expMatrix[,-1]
  colnames(expMatrix)<-cellsName
  
  if(length){
    return(list(expression= expMatrix, biases = biases, 
                DE.Genes = DEvector, ExpectedCount = expVector, Lengths=sizeVector))
  }
  return(list(expression= expMatrix, biases = biases, 
              DE.Genes = DEvector, ExpectedCount = expVector))
}
