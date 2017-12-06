### simulateCountDropout : this function extend data simulation by Lun et al. adding
###     amplification biases, gene length, variations in nb size parameter.
### parameters :
###    - nbgenes : total number of simulated genes
###    - nbDEgenes : number of differentially expressed genes
###    - nbpop : number of population defined by DE genes
###    - nbcells : number of cells by population
###    - pup : probablity for a gene to be upregulated
###    - fcup : fold change for upregulated genes
###    - fcdown : fold change for downregulated genes
###    - cellBiases : weither to add or not cell specific bias
###    - ampliBiases : weither to add or not amplification bias
###    - length : weither to add or not genes length
###    - sizesDisp : dispertion of gaussian distribution from which cell biases are drawn
###    - popsizes : weither to use or not different size parameters for each population negative binomial
###    - nbsize : size parameter of the negative binomial (can be a vector)
###    - g.shape : shape parameter for the gamma distribution of mean expression
###    - g.rate : rate parameter for the gamma distribution of mean expression
###    - dropout : weither to use or not dropout simulation
###    - dropout.rate : probability of dropout (default value will mimic observations by Teichlab)


simulateCountDropout<-function(nbgenes=10000, nbDEgenes = 3000 ,nbpop = 3, 
                           nbcells = 100, pup = 0.4, fcup = 5, fcdown = 0,
                           cellBiases=TRUE, ampliBiases=TRUE, length=TRUE,
                           sizesDisp = 0.25, popsizes=TRUE, nbsize=0.1, g.shape=2, g.rate=2,
                           dropout=TRUE, dropout.rate=1/c(rep(1:3,1800), rep(100:300,20), 2400:4000)) {
  
  # Check parameters
  if(fcdown > 1) stop("Fold-change for down regulated genes must be less than 1.")
  if(nbDEgenes > nbgenes/nbpop) stop("Too many differencially expressed genes in each population")
  if(fcup < 1) stop("Fold change for upregulated genes must be more than 1.")
  if(pup > 1) stop("More than 100% of DE genes upregulated.")
  
  # Initialize return variables
  expMatrix<-c()
  bioMatrix<-c()
  biases<-c()
  drop.rate<-1:nbgenes
  
  # Creating gene caracteristics
  expVector<-rgamma(nbgenes, shape=g.shape, rate=g.rate)
    #sizes come from Mus Musculus genomes, and account for almost 95% of transcripts
  sizeVector<-sample(150:20000, nbgenes, replace=TRUE)
  nbSizes<-sample(nbsize, nbgenes, replace=TRUE)
  
 
  
  # Select differentially expressed genes
  DEvector<-rep(0, nbgenes)
  genesInd<-c(1:nbgenes)
  for (p in 1: nbpop) {
    DEgenesInd<-sample(genesInd, nbDEgenes, replace=FALSE)
    genesInd<-genesInd[-DEgenesInd]
    DEvector[DEgenesInd]<-p
    if(popsizes) nbSizes<-cbind(nbSizes, sample(nbsize, nbgenes, replace=TRUE))
  }
  if(popsizes) nbSizes<-nbSizes[,-1]
  else nb.sizes=nbSizes
  
  cellsName<-c()
  expMatrix<-seq(from=1, to=nbgenes)
  bioMatrix<-1:nbgenes
  dropedMatrix<-1:nbgenes
  lengthMatrix<-1:nbgenes
  
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
    if(popsizes) nb.sizes=nbSizes[,p]
    # Adding cells caracteristics
    for(c in 1:nbcells) {
      cellBias<-1
      cellsName<-c(cellsName , paste0("Pop", p, "-cell", c))
      if(cellBiases) {
        cellBias<-2^(rnorm(1, 0, sizesDisp))
      }
      biases<-c(biases, cellBias)
      cellCount<-sapply(1:nbgenes, function(i) {rnbinom(1, mu=cellBias*popExpCount[i], size=nb.sizes[i]) })
      
      # save "biological" expression 
      bioMatrix<-cbind(bioMatrix, cellCount)
      
      # Dropout : sample for each gene in each cell
      if(dropout){
        rate<-sample(dropout.rate, nbgenes, replace=TRUE)
        drop.rate<-cbind(drop.rate, rate)
        cellCount<-sapply(1:(length(cellCount)), function(i){
          rbinom(1, round(cellCount[i]), rate[i])
        })
      }
      
      # Save detection affected by dropout
      dropedMatrix<-cbind(dropedMatrix, cellCount)
      
      # Length
      if(length) {
        cellCount<-cellCount*sizeVector
      }
      
      # Save reads converted expression
      lengthMatrix<-cbind(lengthMatrix, cellCount)
      
      # Amplification biases
      if(ampliBiases) {
        ampliBias<-rgamma(n=nbgenes, shape=1, rate=2) # initial shape=2.5, lead to many over-amplification
        cellCount<-round(cellCount*ampliBias)
      }
      
      expMatrix<-cbind(expMatrix, cellCount) 
    }
  }

  # Process Matrices
  rownames(expMatrix)<-rownames(bioMatrix)<-rownames(dropedMatrix)<-rownames(lengthMatrix)<-expMatrix[,1]
  expMatrix<-expMatrix[,-1]
  bioMatrix<-bioMatrix[,-1]
  dropedMatrix<-dropedMatrix[,-1]
  lengthMatrix<-lengthMatrix[,-1]
  colnames(expMatrix)<-colnames(bioMatrix)<-colnames(dropedMatrix)<-colnames(lengthMatrix)<-cellsName
  
  # Return simulation results
  l<-list(expression= expMatrix, biases = biases, DE.Genes = DEvector, 
          expected.count = expVector, nb.sizes=nbSizes, bio=bioMatrix, 
          droped=dropedMatrix, reads=lengthMatrix)
  
  if(length) l<-c(l, list(lengths=sizeVector))
  if(dropout) l<-c(l, list(dropout.rate=drop.rate))
  return(l)
}
