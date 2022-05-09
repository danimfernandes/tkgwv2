### Calculate error rates for specific subsets of SNPs by generating simulated individuals based on provided allele frequencies
### Author - Daniel Fernandes
### Instructions:
### Please provide input for the following 3 variables, and then run the whole code.
thins = c(100,500,1000,2000,3000,4000)    # SNP subset sizes
simSeqs = c(100,250,500,1000,2000,4000)   # Number of simulated relationships
freqsFi = "1000GP3_1240K.frq"   # Frequencies file in plink format with header (frq)

simsPair = function(freqFile, numSimPairs=2500, unrelated=TRUE, fullsibs=TRUE, halfsibs=TRUE, freqFileHeader = FALSE, suffix = FALSE) {
  unrelated_vec = c()
  fullSib_vec = c()
  halfSib_vec = c()
  counter = 1
  
  file = freqFile
  dyadi = ""
  
  ################# Read allele frequencies ######
  print(paste0("Loading: ",file))
  alFreq = read.csv(file,stringsAsFactors=FALSE,sep="",colClasses = c(character(),character(),character(),character(),numeric(),numeric()),header=freqFileHeader)
  if(freqFileHeader == FALSE){      colnames(alFreq) = c("CHR", "SNP", "A1",  "A2",  "MAF", "NCHROBS") }
  if(length(which(0 == alFreq$A2)) > 0){      alFreq = alFreq[-which(0 == alFreq$A2),] }
  alFreq$NCHROBS=NULL
  ## Remove SNPs with fixed alleles or missing data
  alFreq = alFreq[alFreq$MAF != 0,]
  alFreq$AFa2=format((1-as.numeric(alFreq$MAF)),digits=3)
  alFreq[5]=lapply(alFreq[5],round,3)
  ##Work with 3 decimal places for frequencies and make sure MAF+AFa2=1
  maf=format(as.numeric(alFreq$MAF),signif=3)
  alFreq$MAF=maf
  af2=format((1-as.numeric(maf)),signif=3)
  alFreq$AFa2=af2
  alFreq$CHR=NULL
  
  alFrequ=alFreq
  alFrequ$MAF = as.double(alFrequ$MAF)
  alFrequ$AFa2 = as.double(alFrequ$AFa2)
  
  ######### UNRELATED ##########
  if(unrelated == TRUE) {
    unrelated0 = c()
    unrelated4 = c()
    unrelatedratio1 = c()
    unrelatedratio2 = c()
    ## Generate individual+allele columns PLUS choose allele according to RANDOM value
    numIndividuals = numSimPairs*2
    seqIndividuals = seq(1:numIndividuals)
    ## Generate individuals
    library(data.table, quietly = T)
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
  }
  
  ######### FULL SIBLINGS ##########
  if(fullsibs == TRUE) {
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
  }
  
  ######### HALF SIBLINGS ##########
  if(halfsibs == TRUE) {
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
    library(data.table, quietly = T)
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
  }
  if( suffix == FALSE) {
    suffix = ""
  }
  return(list(c(fullSib_vec,halfSib_vec,unrelated_vec),c(round(length(which(fullSib_vec < 0.1875)) / sim,digits = 5),round((length(which(halfSib_vec > 0.1875)) + length(which(halfSib_vec < 0.0625))) / sim,digits = 5),round(length(which(unrelated_vec > 0.0625)) / sim,digits = 5))))
}

df1 = data.frame(matrix(nrow = length(thins),ncol = length(simSeqs)+1))
colnames(df1) = c("SNPs/Sims",simSeqs)
df2 = data.frame(matrix(nrow = length(thins),ncol = length(simSeqs)+1))
colnames(df2) = c("SNPs/Sims",simSeqs)
dfU = data.frame(matrix(nrow = length(thins),ncol = length(simSeqs)+1))
colnames(dfU) = c("SNPs/Sims",simSeqs)
counter = 1
for(i in thins) {
  comm1 = paste0("shuf -n ",i," ",freqsFi," > Freqs_",i,"sub.frq")
  system(comm1)
  freqsNew = paste0("Freqs_",i,"sub.frq")
  print(paste0("SNPS: ",i))
  df1[counter,"SNPs/Sims"] = i
  df2[counter,"SNPs/Sims"] = i
  dfU[counter,"SNPs/Sims"] = i
 
  for(sim in simSeqs) {
    print(paste0("Simulated relationships: ",sim))
    errorRates = simsPair(freqFile = freqsNew, numSimPairs = sim, suffix = paste0(i,"SNP_",sim))
    
    df1[counter,paste(sim)] = errorRates[[2]][1]
    df2[counter,paste(sim)] = errorRates[[2]][2]
    dfU[counter,paste(sim)] = errorRates[[2]][3]
  }
  counter = counter + 1
}

write.table(rbind("1stD",df1,"2ndD",df2,"Unr",dfU),"SimulationErrorRates.txt", sep="\t", col.names = T, row.names = F, quote = F)
pdf(file = "SimulationErrorRates.pdf",width = 10, height = 8, pointsize = 1)
par(mar = c(5,5,1,1))
plot(-1,type = "b",xlim = c(0,max(thins)), ylim = c(0,0.5), xlab="",ylab="",axes = F)
axis(2,c(0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.4,0.45,0.5),lwd=0,lwd.ticks = 1, cex.axis =1.8, cex.lab = 2)
mtext("Fraction of erroneously assigned relationships", side=2, line=3.3, cex=1.8)
axis(1,thins,lwd=0,lwd.ticks=1, cex.axis = 1.8)
mtext("SNPs used", side=1, line=3.3, cex=1.8)
box()

#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#colUse = sample(color, length(simSeqs))
library(RColorBrewer)
colUse = brewer.pal(length(simSeqs), "Accent")

cit=1
for(sim in simSeqs) {
  lines(x = df1$SNPs,y = df1[,paste(sim)], lty = 3, lwd = 2, col = colUse[cit])
  lines(x = df2$SNPs,y = df2[,paste(sim)], lty = 2, lwd = 2, col = colUse[cit])
  lines(x = dfU$SNPs,y = dfU[,paste(sim)], lty = 1, lwd = 2, col = colUse[cit])
  cit = cit+1
}

abline(h=0.01,lty=3)
if(length(as.vector(which(rowMeans(df2[,2:length(colnames(df2))]) < 0.01))) != 0) {
  abline(v=df2[as.vector(which(rowMeans(df2[,2:length(colnames(df2))]) < 0.01))[1],1],col="grey65",lwd=2)
}

legend("topright",legend = c("# Simulations",c(paste(simSeqs)),"","Unrelated","2nd-degree","1st-degree"), lty = c(NA,rep(NA,length(simSeqs)),NA,1,2,3), pch = c(NA,rep(22,length(simSeqs)),NA,NA,NA,NA), pt.cex = 5.5, cex = 1.8, pt.bg=c(NA,colUse))
dev.off()
