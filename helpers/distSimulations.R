### Calculate distribution ranges for the provided allele frequencies and calculate posterior probabilities for each class
### Author - Daniel Fernandes (with posterior probability code by John Finarelli from 'tkrelated' package of Fernandes et al. 2017 - doi:10.1038/srep41529)
### Instructions:
### On line 494, you need to provide input for at least the compulsory argument 'sampleVec' when calling the function. Run the whole code afterwards.

distSimulations = function(sampleVec, numSimPairs=2000, freqFileHeader = FALSE) {
  # Description of arguments #
  # sampleVec = Vector of input frequency files to work on ("commSAMP1____SAMP2.frq")
  # numSimPairs = Default 2000. Number of simulated pairs of individuals per class to be generated
  # freqFileHeader = Default FALSE. When false a header will be added to the loaded .frq file
  
  for(samInVec in sampleVec) {
    unrelated_vec = c()
    fullSib_vec = c()
    halfSib_vec = c()
    counter = 1
    
    file = samInVec
    dyadi = paste0(strsplit(strsplit(file,"comm")[[1]][2],".frq")[[1]][1],"__",numSimPairs)
    samp1 = strsplit(strsplit(strsplit(file,"comm")[[1]][2],".frq")[[1]][1], "____")[[1]][1]
    samp2 = strsplit(strsplit(strsplit(file,"comm")[[1]][2],".frq")[[1]][1], "____")[[1]][2]
    
    ## Read allele frequencies
    alFreq = read.csv(file,stringsAsFactors=FALSE,sep="",colClasses = c(character(),character(),character(),character(),numeric(),numeric()),header=freqFileHeader)
    if(freqFileHeader == FALSE){      colnames(alFreq) = c("CHR", "SNP", "A1",  "A2",  "MAF", "NCHROBS") }
    if(length(which(0 == alFreq$A2)) > 0){      alFreq = alFreq[-which(0 == alFreq$A2),] }
    alFreq$NCHROBS=NULL
    ## Remove SNPs with fixed alleles or missing data
    alFreq = alFreq[alFreq$MAF != 0,]
    alFreq$AFa2=format((1-as.numeric(alFreq$MAF)),digits=3)
    alFreq[5]=lapply(alFreq[5],round,3)
    ## Work with 3 decimal places for frequencies and make sure MAF+AFa2=1
    maf=format(as.numeric(alFreq$MAF),signif=3)
    alFreq$MAF=maf
    af2=format((1-as.numeric(maf)),signif=3)
    alFreq$AFa2=af2
    alFreq$CHR=NULL
    
    alFrequ=alFreq
    alFrequ$MAF = as.double(alFrequ$MAF)
    alFrequ$AFa2 = as.double(alFrequ$AFa2)
    
    print(paste0("Generating simulated individuals based on allele frequencies from ",samInVec))
    
    library(data.table, quietly = T)
    
    ## Unrelated
    unrelated0 = c()
    unrelated4 = c()
    unrelatedratio1 = c()
    unrelatedratio2 = c()
    ## Generate individual+allele columns PLUS choose allele according to RANDOM value
    numIndividuals = numSimPairs*2
    seqIndividuals = seq(1:numIndividuals)
    ## Generate individuals
    nm1 <- paste0("IND", rep(seqIndividuals, each = 2), rep(letters[1:2], length(seqIndividuals)))
    system.time({
      setDT(alFrequ)
      for(j in seq_along(nm1)){
        alFrequ[, nm1[j] := A2][runif(.N, 0, 1) < MAF , nm1[j] := A1][]
      }
    })
    setDF(alFrequ)
    
    ## Create new data.frame with forced homozygous
    alFrequFH = alFrequ
    remove(alFrequ)
    library(data.table, quietly = T)
    nm2 = c(colnames(alFrequFH)[1:5],paste0("IND", rep(seqIndividuals)))
    nm1 = colnames(alFrequFH)[6:length(colnames(alFrequFH))]
    system.time({
      setDT(alFrequFH)
      for(j in seq(from=1,to=numSimPairs*4,by=2)){
        alFrequFH[, "a1" := alFrequFH[[nm1[j]]]]
        alFrequFH[, "b1" := alFrequFH[[nm1[j+1]]]]
        alFrequFH[, nm1[j] := a1][runif(.N, 0, 1) < 0.5 , nm1[j] := b1][]
        alFrequFH[, nm1[j+1] := NULL]
      }
    })
    setDF(alFrequFH)
    alFrequFH = alFrequFH[1:(length(colnames(alFrequFH))-2)]
    colnames(alFrequFH) = nm2
    
    ## Estimates relatedness using Queller and Goonight's (1989) Rxy estimator
    dyads = seq(1,length(alFrequFH)-5,by=2)
    iter = dyads[1]
    for(dyad in dyads) {
      ind1=paste0("IND",dyad)
      ind2=paste0("IND",dyad+1)
      commonDataFrameRxy = alFrequFH[1:5]
      
      ## Add one combination to dataframe
      commonDataFrameRxy$S1x = eval(alFrequFH[[ind1]])
      commonDataFrameRxy$S2x = eval(alFrequFH[[ind2]])
      
      ## Remove SNPs with fixed alleles
      commonDataFrameRxy = commonDataFrameRxy[commonDataFrameRxy$MAF != 0,]
      
      ## Make sure all variants are SNPs
      ## Remove them from alFreq and commonDataFrame
      rem1 = which(nchar(commonDataFrameRxy$A1) >1)
      rem2 = which(nchar(commonDataFrameRxy$A2) >1)
      rem3 = c(rem1,rem2)
      remSnps = commonDataFrameRxy[rem3,2]
      if(length(rem3) >0) { commonDataFrameRxy = commonDataFrameRxy[-rem3,] }
      
      ## Calculate components
      commonDataFrameRxy["PaXy"] = ifelse(as.character(commonDataFrameRxy$S1x) == commonDataFrameRxy$A1, commonDataFrameRxy$MAF, commonDataFrameRxy$AFa2)
      commonDataFrameRxy["PcXy"] = ifelse(as.character(commonDataFrameRxy$S2x) == commonDataFrameRxy$A1, commonDataFrameRxy$MAF, commonDataFrameRxy$AFa2)
      commonDataFrameRxy["IacADbcBD"] = ifelse(as.character(commonDataFrameRxy$S1x) == as.character(commonDataFrameRxy$S2x), 1, 0)*4
      
      commonDataFrameRxy$rxyL = ((0.5*(as.numeric(commonDataFrameRxy$IacADbcBD)))-(as.numeric(commonDataFrameRxy$PaXy)*2)) / (1+1-(as.numeric(commonDataFrameRxy$PaXy)*2))
      commonDataFrameRxy$ryxL = ((0.5*(as.numeric(commonDataFrameRxy$IacADbcBD)))-(as.numeric(commonDataFrameRxy$PcXy)*2)) / (1+1-(as.numeric(commonDataFrameRxy$PcXy)*2))
      rxyL = sum(commonDataFrameRxy$rxyL)
      ryxL = sum(commonDataFrameRxy$ryxL)
      commonDataFrameRxy$Av = (commonDataFrameRxy$rxyL + commonDataFrameRxy$ryxL) /2
      commonDataFrameRxy$SrQ1 = ((0.5*(as.numeric(commonDataFrameRxy$IacADbcBD)))-(as.numeric(commonDataFrameRxy$PaXy)*2))
      commonDataFrameRxy$SrQ2 = ((0.5*(as.numeric(commonDataFrameRxy$IacADbcBD)))-(as.numeric(commonDataFrameRxy$PcXy)*2))
      commonDataFrameRxy$SrQ1MinusSum = commonDataFrameRxy$SrQ1 - sum(commonDataFrameRxy$SrQ1)
      commonDataFrameRxy$SrQ2MinusSum = commonDataFrameRxy$SrQ2 - sum(commonDataFrameRxy$SrQ2)
      commonDataFrameRxy$WrQ1 = (1+1-(as.numeric(commonDataFrameRxy$PaXy)*2))
      commonDataFrameRxy$WrQ2 = (1+1-(as.numeric(commonDataFrameRxy$PcXy)*2))
      commonDataFrameRxy$WrQ1MinusSum = commonDataFrameRxy$WrQ1 - sum(commonDataFrameRxy$WrQ1)
      commonDataFrameRxy$WrQ2MinusSum = commonDataFrameRxy$WrQ2 - sum(commonDataFrameRxy$WrQ2)
      a=mean(commonDataFrameRxy$SrQ1MinusSum/commonDataFrameRxy$WrQ1MinusSum)
      b=mean(commonDataFrameRxy$SrQ2MinusSum/commonDataFrameRxy$WrQ2MinusSum)
      relatednessR=(a+b)/2
      unrelated_vec = c(unrelated_vec,relatednessR)
      unrelated0 = c(unrelated0, length(which(0 == commonDataFrameRxy$IacADbcBD)))
      unrelated4 = c(unrelated4, length(which(4 == commonDataFrameRxy$IacADbcBD)))
      
    }
    unrelatedratio1 = c(unrelated0/unrelated4)
    unrelatedratio2 = c(unrelated4/unrelated0)
    
    ## Full Siblings - 1st degree
    fullSibs0 = c()
    fullSibs4 = c()
    fullSibsratio1 = c()
    fullSibsratio2 = c()
    
    alFrequ=alFreq
    alFrequ$MAF = as.double(alFrequ$MAF)
    alFrequ$AFa2 = as.double(alFrequ$AFa2)
    ## Generate individual+allele columns PLUS choose allele according to RANDOM value
    numIndividuals = numSimPairs*2 
    seqIndividuals = seq(1:numIndividuals)
    
    library(data.table, quietly = T)
    nm1 <- paste0("IND", rep(seqIndividuals, each = 2), rep(letters[1:2], length(seqIndividuals)))
    system.time({
      setDT(alFrequ)
      for(j in seq_along(nm1)){
        alFrequ[, nm1[j] := A2][runif(.N, 0, 1) < MAF , nm1[j] := A1][]
      }
    })
    setDF(alFrequ)
    
    ## Choosing one allele at random from parent
    namesFS = paste0("FSIB", rep(seq(from=1,to=length(seqIndividuals),by=2),each=4),"-",rep(seq(from=2,to=length(seqIndividuals),by=2),each=4),rep(c("I","II"),each=2,times=numSimPairs),"_",rep(letters[1:2], length(seqIndividuals)-2))
    fullSib = data.frame(matrix(data="N",ncol = numSimPairs*4,nrow = length(alFrequ$SNP)))
    colnames(fullSib) = namesFS
    system.time({
      setDT(alFrequ)
      setDT(fullSib)
      for(j in seq(from=1,to=numSimPairs*4,by=4)){
        fullSib[, "a1" := alFrequ[[nm1[j]]]]
        fullSib[, "b1" := alFrequ[[nm1[j+1]]]]
        fullSib[, "a2" := alFrequ[[nm1[j+2]]]]
        fullSib[, "b2" := alFrequ[[nm1[j+3]]]]
        fullSib[, namesFS[j] := a1][runif(.N, 0, 1) < 0.5 , namesFS[j] := b1 ][]
        fullSib[, namesFS[j+1] := a2][runif(.N, 0, 1) < 0.5 , namesFS[j+1] := b2][]
        fullSib[, namesFS[j+2] := a1][runif(.N, 0, 1) < 0.5 , namesFS[j+2] := b1 ][]
        fullSib[, namesFS[j+3] := a2][runif(.N, 0, 1) < 0.5 , namesFS[j+3] := b2][]
      }
    })
    setDF(alFrequ)
    setDF(fullSib)
    fullSib = fullSib[1:(numSimPairs*4)]
    
    ## Create a list of the pairs of siblings
    listFullSibs = list()
    numSiblingPairs = numSimPairs
    seqSiblings = seq(1:numSiblingPairs)
    iter=1
    for(i in seqSiblings) {
      nextInd=iter+1
      ind = paste0("FSIB",iter,"-",nextInd,"I") 
      ind2 = paste0("FSIB",iter,"-",nextInd,"II") 
      listFullSibs = c(listFullSibs, ind,ind2)
      iter=iter+2
    }
    
    ## Force homozygous
    fullSibFH = cbind(alFrequ[,1:5],fullSib)
    remove(alFrequ)                         
    library(data.table, quietly = T)
    nm2 = c(colnames(fullSibFH)[1:5],paste0("FSIB", rep(seq(from=1,to=length(seqIndividuals),by=2),each=2),"-",rep(seq(from=2,to=length(seqIndividuals),by=2),each=2),rep(c("I","II"),each=1,times=numSimPairs)))
    nm1 = colnames(fullSibFH)[6:length(colnames(fullSibFH))]
    system.time({
      setDT(fullSibFH)
      for(j in seq(from=1,to=numSimPairs*4,by=2)){
        fullSibFH[, "a1" := fullSibFH[[nm1[j]]]]
        fullSibFH[, "b1" := fullSibFH[[nm1[j+1]]]]
        fullSibFH[, nm1[j] := a1][runif(.N, 0, 1) < 0.5 , nm1[j] := b1][]
        fullSibFH[, nm1[j+1] := NULL]
      }
    })
    setDF(fullSibFH)
    fullSibFH = fullSibFH[1:(length(colnames(fullSibFH))-2)]
    colnames(fullSibFH) = nm2
    
    ## Estimates relatedness using Queller and Goonight's (1989) Rxy estimator
    dyadsOther = as.character(listFullSibs)[c(T,F)]
    iter = 1
    for(dyad in dyadsOther) {
      ind1=dyad
      ind2=as.character(listFullSibs[iter+1])
      iter=iter+2
      commonDataFrameRxy = fullSibFH[1:5]
      
      ## Add one combination to dataframe
      commonDataFrameRxy$S1x = eval(fullSibFH[[ind1]])
      commonDataFrameRxy$S2x = eval(fullSibFH[[ind2]])
      
      ## Calculate components
      commonDataFrameRxy["PaXy"] = ifelse(as.character(commonDataFrameRxy$S1x) == commonDataFrameRxy$A1, commonDataFrameRxy$MAF, commonDataFrameRxy$AFa2)
      commonDataFrameRxy["PcXy"] = ifelse(as.character(commonDataFrameRxy$S2x) == commonDataFrameRxy$A1, commonDataFrameRxy$MAF, commonDataFrameRxy$AFa2)
      commonDataFrameRxy["IacADbcBD"] = ifelse(as.character(commonDataFrameRxy$S1x) == as.character(commonDataFrameRxy$S2x), 1, 0)*4
      
      commonDataFrameRxy$rxyL = ((0.5*(as.numeric(commonDataFrameRxy$IacADbcBD)))-(as.numeric(commonDataFrameRxy$PaXy)*2)) / (1+1-(as.numeric(commonDataFrameRxy$PaXy)*2))
      commonDataFrameRxy$ryxL = ((0.5*(as.numeric(commonDataFrameRxy$IacADbcBD)))-(as.numeric(commonDataFrameRxy$PcXy)*2)) / (1+1-(as.numeric(commonDataFrameRxy$PcXy)*2))
      rxyL = sum(commonDataFrameRxy$rxyL)
      ryxL = sum(commonDataFrameRxy$ryxL)
      commonDataFrameRxy$Av = (commonDataFrameRxy$rxyL + commonDataFrameRxy$ryxL) /2
      commonDataFrameRxy$SrQ1 = ((0.5*(as.numeric(commonDataFrameRxy$IacADbcBD)))-(as.numeric(commonDataFrameRxy$PaXy)*2))
      commonDataFrameRxy$SrQ2 = ((0.5*(as.numeric(commonDataFrameRxy$IacADbcBD)))-(as.numeric(commonDataFrameRxy$PcXy)*2))
      commonDataFrameRxy$SrQ1MinusSum = commonDataFrameRxy$SrQ1 - sum(commonDataFrameRxy$SrQ1)
      commonDataFrameRxy$SrQ2MinusSum = commonDataFrameRxy$SrQ2 - sum(commonDataFrameRxy$SrQ2)
      commonDataFrameRxy$WrQ1 = (1+1-(as.numeric(commonDataFrameRxy$PaXy)*2))
      commonDataFrameRxy$WrQ2 = (1+1-(as.numeric(commonDataFrameRxy$PcXy)*2))
      commonDataFrameRxy$WrQ1MinusSum = commonDataFrameRxy$WrQ1 - sum(commonDataFrameRxy$WrQ1)
      commonDataFrameRxy$WrQ2MinusSum = commonDataFrameRxy$WrQ2 - sum(commonDataFrameRxy$WrQ2)
      a=mean(commonDataFrameRxy$SrQ1MinusSum/commonDataFrameRxy$WrQ1MinusSum)
      b=mean(commonDataFrameRxy$SrQ2MinusSum/commonDataFrameRxy$WrQ2MinusSum)
      relatednessR=(a+b)/2
      fullSib_vec = c(fullSib_vec,relatednessR)
      fullSibs0 = c(fullSibs0, length(which(0 == commonDataFrameRxy$IacADbcBD)))
      fullSibs4 = c(fullSibs4, length(which(4 == commonDataFrameRxy$IacADbcBD)))
    }
    fullSibsratio1 = c(fullSibs0/fullSibs4)
    fullSibsratio2 = c(fullSibs4/fullSibs0)
        
    ## Half siblings - 2nd degree
    halfSibs0 = c()
    halfSibs4 = c()
    halfSibsratio1 = c()
    halfSibsratio2 = c()
    
    alFrequ=alFreq
    alFrequ$MAF = as.double(alFrequ$MAF)
    alFrequ$AFa2 = as.double(alFrequ$AFa2)
    ## Generate individual+allele columns PLUS choose allele according to RANDOM value
    numIndividuals = numSimPairs*2
    numIndividuals = numIndividuals +1 ##Correcting for permutations
    seqIndividuals = seq(1:numIndividuals)
    
    library(data.table, quietly = T)
    nm1 <- paste0("IND", rep(seqIndividuals, each = 2), rep(letters[1:2], length(seqIndividuals)))
    system.time({
      setDT(alFrequ)
      for(j in seq_along(nm1)){
        alFrequ[, nm1[j] := A2][runif(.N, 0, 1) < MAF , nm1[j] := A1][]
      }
    })
    setDF(alFrequ)
    
    ## Choosing one allele at random from parent
    namesHS = paste0("HSIB", rep(seq(from=1,to=length(seqIndividuals)-2,by=2),each=4),"-",rep(seq(from=2,to=length(seqIndividuals),by=1),each=2),rep(c("I","II"),each=2,times=numSimPairs),"_",rep(letters[1:2], length(seqIndividuals)-2))
    halfSib = data.frame(matrix(data="N",ncol = numSimPairs*4,nrow = length(alFrequ$SNP)))
    colnames(halfSib) = namesHS
    system.time({
      setDT(alFrequ)
      setDT(halfSib)
      for(j in seq(from=1,to=numSimPairs*4,by=4)){
        halfSib[, "a1" := alFrequ[[nm1[j]]]]
        halfSib[, "b1" := alFrequ[[nm1[j+1]]]]
        halfSib[, "a2" := alFrequ[[nm1[j+2]]]]
        halfSib[, "b2" := alFrequ[[nm1[j+3]]]]
        halfSib[, "a3" := alFrequ[[nm1[j+4]]]]
        halfSib[, "b3" := alFrequ[[nm1[j+5]]]]
        halfSib[, namesHS[j] := a1][runif(.N, 0, 1) < 0.5 , namesHS[j] := b1 ][]
        halfSib[, namesHS[j+1] := a2][runif(.N, 0, 1) < 0.5 , namesHS[j+1] := b2][]
        halfSib[, namesHS[j+2] := a1][runif(.N, 0, 1) < 0.5 , namesHS[j+2] := b1 ][]
        halfSib[, namesHS[j+3] := a3][runif(.N, 0, 1) < 0.5 , namesHS[j+3] := b3][]
      }
    })
    setDF(alFrequ)
    setDF(halfSib)
    halfSib = halfSib[1:(numSimPairs*4)]
    
    ## Create a list of the pairs of siblings
    listHalfSibs = list()
    numSiblingPairs = numSimPairs
    seqSiblings = seq(1:numSiblingPairs)
    iter=1
    for(i in seqSiblings) {
      nextInd=iter+1
      nextNextInd=nextInd+1
      ind = paste0("HSIB",iter,"-",nextInd,"I") 
      ind2 = paste0("HSIB",iter,"-",nextNextInd,"II") 
      listHalfSibs = c(listHalfSibs, ind,ind2)
      iter=iter+2
    }
    
    ## Force homozygous for Rxy/2
    halfSibFH = cbind(alFrequ[,1:5],halfSib)
    nm2 = c(colnames(halfSibFH)[1:5],paste0("HSIB", rep(seq(from=1,to=length(seqIndividuals)-2,by=2),each=2),"-",rep(seq(from=2,to=length(seqIndividuals),by=1),each=1),rep(c("I","II"),each=1,times=numSimPairs)))
    nm1 = colnames(halfSibFH)[6:length(colnames(halfSibFH))]
    system.time({
      setDT(halfSibFH)
      for(j in seq(from=1,to=numSimPairs*4,by=2)){
        halfSibFH[, "a1" := halfSibFH[[nm1[j]]]]
        halfSibFH[, "b1" := halfSibFH[[nm1[j+1]]]]
        halfSibFH[, nm1[j] := a1][runif(.N, 0, 1) < 0.5 , nm1[j] := b1][]
        halfSibFH[, nm1[j+1] := NULL]
      }
    })
    setDF(halfSibFH)
    halfSibFH = halfSibFH[1:(length(colnames(halfSibFH))-2)]
    colnames(halfSibFH) = nm2
    
    ## Estimates relatedness using Queller and Goonight's (1989) Rxy estimator
    dyadsOther = as.character(listHalfSibs)[c(T,F)]
    iter = 1
    for(dyad in dyadsOther) {
      ind1=dyad
      ind2=as.character(listHalfSibs[iter+1])
      iter=iter+2
      commonDataFrameRxy = halfSibFH[1:5]
      
      ## Add one combination to dataframe
      commonDataFrameRxy$S1x = eval(halfSibFH[[ind1]])
      commonDataFrameRxy$S2x = eval(halfSibFH[[ind2]])
      
      ## Calculate components
      commonDataFrameRxy["PaXy"] = ifelse(as.character(commonDataFrameRxy$S1x) == commonDataFrameRxy$A1, commonDataFrameRxy$MAF, commonDataFrameRxy$AFa2)
      commonDataFrameRxy["PcXy"] = ifelse(as.character(commonDataFrameRxy$S2x) == commonDataFrameRxy$A1, commonDataFrameRxy$MAF, commonDataFrameRxy$AFa2)
      commonDataFrameRxy["IacADbcBD"] = ifelse(as.character(commonDataFrameRxy$S1x) == as.character(commonDataFrameRxy$S2x), 1, 0)*4
      
      commonDataFrameRxy$rxyL = ((0.5*(as.numeric(commonDataFrameRxy$IacADbcBD)))-(as.numeric(commonDataFrameRxy$PaXy)*2)) / (1+1-(as.numeric(commonDataFrameRxy$PaXy)*2))
      commonDataFrameRxy$ryxL = ((0.5*(as.numeric(commonDataFrameRxy$IacADbcBD)))-(as.numeric(commonDataFrameRxy$PcXy)*2)) / (1+1-(as.numeric(commonDataFrameRxy$PcXy)*2))
      rxyL = sum(commonDataFrameRxy$rxyL)
      ryxL = sum(commonDataFrameRxy$ryxL)
      commonDataFrameRxy$Av = (commonDataFrameRxy$rxyL + commonDataFrameRxy$ryxL) /2
      commonDataFrameRxy$SrQ1 = ((0.5*(as.numeric(commonDataFrameRxy$IacADbcBD)))-(as.numeric(commonDataFrameRxy$PaXy)*2))
      commonDataFrameRxy$SrQ2 = ((0.5*(as.numeric(commonDataFrameRxy$IacADbcBD)))-(as.numeric(commonDataFrameRxy$PcXy)*2))
      commonDataFrameRxy$SrQ1MinusSum = commonDataFrameRxy$SrQ1 - sum(commonDataFrameRxy$SrQ1)
      commonDataFrameRxy$SrQ2MinusSum = commonDataFrameRxy$SrQ2 - sum(commonDataFrameRxy$SrQ2)
      commonDataFrameRxy$WrQ1 = (1+1-(as.numeric(commonDataFrameRxy$PaXy)*2))
      commonDataFrameRxy$WrQ2 = (1+1-(as.numeric(commonDataFrameRxy$PcXy)*2))
      commonDataFrameRxy$WrQ1MinusSum = commonDataFrameRxy$WrQ1 - sum(commonDataFrameRxy$WrQ1)
      commonDataFrameRxy$WrQ2MinusSum = commonDataFrameRxy$WrQ2 - sum(commonDataFrameRxy$WrQ2)
      a=mean(commonDataFrameRxy$SrQ1MinusSum/commonDataFrameRxy$WrQ1MinusSum)
      b=mean(commonDataFrameRxy$SrQ2MinusSum/commonDataFrameRxy$WrQ2MinusSum)
      relatednessR=(a+b)/2
      halfSib_vec = c(halfSib_vec,relatednessR)
      halfSibs0 = c(halfSibs0, length(which(0 == commonDataFrameRxy$IacADbcBD)))
      halfSibs4 = c(halfSibs4, length(which(4 == commonDataFrameRxy$IacADbcBD)))
    }
    halfSibsratio1 = c(halfSibs0/halfSibs4)
    halfSibsratio2 = c(halfSibs4/halfSibs0)
    
    histoUN = unrelated_vec
    histoHS = halfSib_vec
    histoFS = fullSib_vec
    
    if(max(histoUN) < min(histoHS) || max(histoHS) < min(histoFS))  {print("Warning: Some of the curves do not seem to be overlapping. It is advised to re-run the script with a larger number of simulated individuals for more accurate probabilities.")}
    
    n_size = length(histoUN)+length(histoHS)+length(histoFS)
    
    ## Read in the relatedness coefficient for this pair of samples and number of SNPs used from "TKGWV2_Results.txt"
    fileIn = read.csv("TKGWV2_Results.txt",stringsAsFactors = F,sep="\t")
    fileIn$Merge = paste0(fileIn$Sample1,"____",fileIn$Sample2)
    coefficientRelatedness = fileIn[which(strsplit(dyadi,paste0("__",numSimPairs))[[1]][1] == fileIn$Merge),paste("HRC")]
    commsnps = fileIn[which(strsplit(dyadi,paste0("__",numSimPairs))[[1]][1] == fileIn$Merge),paste("Used_SNPs")]
    
    roundUp <- function(x)  {round(x+5,-1)}
    pdf(file=paste0(getwd(),"/",dyadi,".pdf"),width=10,height=8, pointsize=3,onefile = FALSE) 
    par(mar=c(5.5,5.5,3,2),mgp=c(3.2,2,1))
    
    allTemp = data.frame(matrix(ncol=1,nrow = length(histoUN)*3))
    allTemp[,1] = c(histoFS,histoHS,histoUN)
    breaksAll = seq(min(allTemp[,1]),max(allTemp[,1]),0.01)
    hist(allTemp[,1],xlim=c(-0.15,0.40),col=rgb(42,77,110,190, maxColorValue = 255), breaks=length(breaksAll), xlab="", ylab="", main=paste0("Distribution ranges for ",samp1,"____",samp2," using ",numSiblingPairs," simulated relationships per class"),xaxt="n")
    if(length(histoUN) > 0) {
      breaksUN = seq(min(histoUN),max(histoUN),0.01)
      hist(histoUN,xlim=c(-0.15,0.40),col=rgb(42,77,110,190, maxColorValue = 255), breaks=length(breaksUN), xlab="", ylab="", main=paste0("Distribution ranges for ",samp1,"____",samp2," using ",numSiblingPairs," simulated relationships per class"),xaxt="n", yaxt="n", ylim=c(0,roundUp(par("usr")[4])),cex.main=2)
    }
    if(length(histoFS) > 0) {
      breaksFS = seq(min(histoFS),max(histoFS),0.01)
      hist(histoFS,add=TRUE,col=rgb(168,56,59,190, maxColorValue = 255),breaks=length(breaksFS))
    }
    if(length(histoHS) > 0) {
      breaksHS = seq(min(histoHS),max(histoHS),0.01)
      hist(histoHS,add=TRUE,col=rgb(154,166,55,190,maxColorValue = 255),breaks=length(breaksHS))
    }
    axis(side=1, at=c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4),cex.axis =1.8, cex.lab = 2)
    mtext("Halved Relatedness Coefficient", side=1, line=4, cex=1.8)
    axis(side=2,cex.axis =1.8, cex.lab = 2)
    mtext(paste0("N=",as.numeric(n_size)), side=2, line=4, cex=1.8)
    
    ## Define small database with data to use only for the orders being used
    dbTable = data.frame(matrix(nrow = 3,ncol= 7))
    dbTable$X1 = c("Unrelated","Second Order","First Order")
    dbTable$X2 = c(rgb(42,77,110,190, maxColorValue = 255),rgb(154,166,55,190,maxColorValue = 255),rgb(168,56,59,190, maxColorValue = 255))
    dbTable$X3 = c("UN","HS","FS")
    dbTable$X4 = c("UNnl","HSnl","FSnl")
    dbTable$X5 = c(1,2,3)
    exi = c()
    for(i in grep("histo",ls())) {
      if(length(get(ls()[i])) > 0) {
        exi = c(exi,strsplit(ls()[i],"histo")[[1]][2])
      }
    }
    exi = exi[order(match(exi,dbTable$X3))]
    lascolores = dbTable[which(dbTable$X3 %in% exi),][2]
    osnomes = dbTable[which(dbTable$X3 %in% exi),][1]
    
    legend("topright",
           bty="n",
           legend=eval(osnomes)[[1]],
           pch=21,
           pt.bg=eval(lascolores)[[1]],
           col="black",
           y.intersp=1,
           cex=2,
           pt.cex=4)
    
    ## Draw line and text for relatedness coefficient
    segments(x0=coefficientRelatedness, y0=-0.5, y1= par("usr")[4]/20, col="black",lwd=12)
    segments(x0=coefficientRelatedness, y0=-0.5, y1= par("usr")[4]/20, col="gold",lwd=10)
    text(coefficientRelatedness, -par("usr")[4]/50, paste0("HRC=",format(coefficientRelatedness,digits = 3)),cex=2)
    
    ## Calculate posterior probabilities
    normal.likelihood <- function(params,x){
      mu <- params[1]
      var <- params[2]
      logl <- -0.5*log(var) - (1/(2*var)) * ((x-mu)**2)
      return(logl)
    }
    if(length(histoUN) > 0) {
      paramsUN <- c((sum(as.numeric(histoUN))/length(as.numeric(histoUN))),((length(as.numeric(histoUN))-1)/length(as.numeric(histoUN)))*(sd(as.numeric(histoUN))**2))
      UNnl = normal.likelihood(paramsUN,coefficientRelatedness)
    }
    if(length(histoHS) > 0) {
      paramsHS <- c((sum(as.numeric(histoHS))/length(as.numeric(histoHS))),((length(as.numeric(histoHS))-1)/length(as.numeric(histoHS)))*(sd(as.numeric(histoHS))**2))
      HSnl = normal.likelihood(paramsHS,coefficientRelatedness)
    }
    if(length(histoFS) > 0) {
      paramsFS <- c((sum(as.numeric(histoFS))/length(as.numeric(histoFS))),((length(as.numeric(histoFS))-1)/length(as.numeric(histoFS)))*(sd(as.numeric(histoFS))**2))
      FSnl = normal.likelihood(paramsFS,coefficientRelatedness)
    }
    
    log.like = c()
    for(i in dbTable[which(dbTable$X3 %in% exi),][4][,1]) {
      dbTable$X6[dbTable[,5][which(dbTable$X4 %in% i)]] = get(i)
      log.like = c(log.like,get(i))
    }
    
    likelihoods = exp(log.like)
    post.probs = sprintf("%.3f",(as.numeric(likelihoods/(sum(likelihoods)))),3)
    
    counter = 1
    for(i in exi) {
      dbTable$X7[which(dbTable$X3 %in% i)] = post.probs[counter]
      counter = counter + 1
    }
    
    ## Add posterior probabilities to plot
    if(length(histoUN) > 0) {
      text(median(histoUN),(round(par("usr")[4])/4),labels = dbTable$X7[1],col="white", cex=2)
    }
    if(length(histoHS) > 0) {
      text(median(histoHS),(round(par("usr")[4])/4),labels = dbTable$X7[2],col="white", cex=2)
    }
    if(length(histoFS) > 0) {
      text(median(histoFS),(round(par("usr")[4])/4),labels = dbTable$X7[3],col="white", cex=2)
    }
    dev.off()
  }
}  

distSimulations(sampleVec = c("commMySamp1____MySamp2.frq", "commMySamp1____MySamp3.frq"))
