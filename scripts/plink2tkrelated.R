### Estimating relatedness from individual text-PLINKs
### Author - Daniel Fernandes
# Programs with system-wide installation required # 
# 1. plink
# R package required # 
# 1. data.table  
#
# Description of arguments #
# dyads = Default TRUE will analyse all possible dyads between all individuals in current folder.
#         Otherwise, user can input 'dyads' as a tab-spaced text file with each line containing a pair to be analysed.
#         Sample1 Sample2
#         Sample1 Sample3
#         (...)
# freqFile = Allele frequencies file in PLINK format (FRQ) for the same SNP set used in bam2plink() (or for example for the 1240K dataset)
# ignoreThresh = Default 1. Threshold for the minimum number of SNPs allowed to estimate relatedness.
# verbose = Default FALSE. Use verbose mode and print extended information throughout file processing.
#
# Input:
# - individual text-PLINK
# - pairwise FRQ
# Output:
# - pairwise relatedness coefficient result text file  

args = commandArgs(trailingOnly=TRUE)
pywd = args[length(args)]
args = args[-length(args)]

if(length(grep("bam2plink",args)) != 0) { 
  args = args[-grep("bam2plink",args)]
} else { 
  cat("\n ################################################################################")
  cat("\n ### TKGWV2 - An ancient DNA relatedness pipeline for ultra-low coverage data ###")
  cat("\n ## Version 1.0b - Released 07/2022")
  cat("\n #")
}

cat(paste0("\n # [",Sys.time(),"] ","Running 'plink2tkrelated' on folder ", getwd()))
cat("\n\t # Text-PLINK >> Pairwise transposed text-PLINK >> Relatedness estimates ")

if(length(grep("freqFile",args)) != 0) {
  freqFile = strsplit(args[grep("freqFile",args)],"=")[[1]][2]
} else {
  cat("\n\t # ERROR: The required argument 'freqFile' could not be found\n")
  quit()
}
dyads = TRUE
if(length(grep("dyads",args)) != 0) { dyads = strsplit(args[grep("dyads",args)],"=")[[1]][2] }

ignoreThresh = 1
if(length(grep("ignoreThresh",args)) != 0) { ignoreThresh = strsplit(args[grep("ignoreThresh",args)],"=")[[1]][2] }

verbose = FALSE
if(length(grep("verbose",args)) != 0) { verbose = TRUE }

if(dyads == TRUE) {
  list_samples = c()
  for(i in list.files(pattern = "\\.ped$")) {
    samp = as.character(strsplit(i,"\\.ped"))
    list_samples = c(list_samples,samp)
  }
  if(length(list_samples) < 2) {
    cat("\n\t # ERROR: Less than 2 '.ped' files found in current folder. A minimum of 2 samples is required\n"); quit()
  } else {combos = data.frame(t(combn(list_samples,2)),stringsAsFactors = F)}
} else {
  combos = read.table(dyads,header = F, stringsAsFactors = F)
}
sampleList = unique(c(combos[,1],combos[,2]))
if(verbose == TRUE) {  cat("######\nMatrix of dyads:\n");  print(combos);  cat("\nList of unique samples:\n");  cat(sampleList);  cat("\n######\n") }

if(length(grep("bam2plink",commandArgs(trailingOnly=TRUE))) != 0) { 
  } else { 
    if(length(list.files(pattern = "\\.ped$")) == 0) {cat("\n\t # ERROR: No '.ped' files found in current folder. A minimum of 2 samples is required\n"); quit()}
    cat("\n\t # Files to be processed:\n")
    if(length(list.files(pattern = "\\.ped$")) == 1) {cat(paste0("\t\t",list.files(pattern = "\\.ped$")[1],"\n\t\tERROR: Only one '.ped' file found in current folder. A minimum of 2 samples is required\n")); quit()}
    if(length(list.files(pattern = "\\.ped$")) > 1) {
      for(i in sampleList) {
        cat(paste0("\t\t ",i))
        if(length(list.files(pattern = paste0("^",i,".ped$"))) == 0) {
          cat(paste0("\tERROR: No '.ped' file ",i," found in current folder\n")); quit()
        }
        if(length(list.files(pattern = paste0("^",i,".map$"))) == 0) {
          cat(paste0("\tERROR: No '.map' file for ",i,".ped found in current folder\n")); quit()
        }
        if(length(list.files(pattern = paste0("^",i,".map$"))) > 0 & length(list.files(pattern = paste0("^",i,".ped$"))) > 0) {
          cat(paste0(paste0(".ped\t",i,".map"),"\n"))
        }
        }
    }
}

cat("\t # Arguments used:\n")
for(i in gsub("=","\t",args)) { cat(paste0("\t\t ",i,"\n"))}
cat("\n")

if(file.exists(freqFile) == FALSE) {
  cat(paste0("\t # ERROR: The frequencies file ",freqFile," could not be found\n")); quit()
}

toExport = data.frame(matrix(nrow = length(combos$X1),ncol=7),stringsAsFactors = F)
colnames(toExport) = c("Sample1","Sample2","Used_SNPs","HRC","counts0","counts4","Relationship")

for(line in seq(1:length(combos[,1]))) {
  samp1 = combos[line,1]
  samp2 = combos[line,2]
  
  cat(paste0("\t # Estimating coefficient of relatedness Rxy for   ",samp1,"   ",samp2,"\t"))
  
  #Find common SNPs
  comm4 = paste0("grep -F -x -f ",samp1,".map ",samp2,".map | awk '{print $2}' > comm",samp1,"_",samp2,"_SNPs")
  system(comm4)
  if(verbose == TRUE) {  cat(paste0("\n\t###### Number of common SNPs for ",samp1,"_",samp2," before QC: ", as.integer(strsplit(system(paste0("wc -l comm",samp1,"_",samp2,"_SNPs"), intern=T)," ")[[1]][1])),"\n") }
  
  if(file.size(paste0("comm",samp1,"_",samp2,"_SNPs")) > 0) {
    #Use PLINK to reduce samples to common SNPs and merge
    comm5 = paste0("plink --file ",samp1," --extract comm",samp1,"_",samp2,"_SNPs --make-bed --out ",samp1,"short --allow-no-sex")
    system(comm5, ignore.stdout = TRUE, ignore.stderr = TRUE)
    comm6 = paste0("plink --file ",samp2," --extract comm",samp1,"_",samp2,"_SNPs --make-bed --out ",samp2,"short --allow-no-sex")
    system(comm6, ignore.stdout = TRUE, ignore.stderr = TRUE)
    comm7 = paste0("plink --bfile ",samp1,"short --bmerge ",samp2,"short --geno 0.1 --recode transpose --out ",samp1,"____",samp2," --allow-no-sex")
    system(comm7, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    #Generate shorter FRQ file
    comm8 = paste0("awk '{print $2}' ",samp1,"____",samp2,".tped > comm",samp1,"_",samp2,"_SNPs2")
    system(comm8)
    if(verbose == TRUE) {  cat(paste0("\t###### Number of common SNPs after QC1: ", as.integer(strsplit(system(paste0("wc -l comm",samp1,"_",samp2,"_SNPs2"), intern=T)," ")[[1]][1])),"\n") }
    #Get frequencies from .frq file
    comm9 = paste0("grep -w -F -f comm",samp1,"_",samp2,"_SNPs2 ",freqFile," > comm",samp1,"____",samp2,".frq")
    system(comm9)
    ## The next commands can be useful when the freqs file is generated by merging individuals in PLINK while excluding variants detected by -missnp
    comm10a = paste0("wc -l comm",samp1,"____",samp2,".frq")
    comm10b = paste0("wc -l ",samp1,"____",samp2,".tped")
    if(as.integer(strsplit(system(comm10a, intern=T)," ")[[1]][1]) != as.integer(strsplit(system(comm10b, intern=T)," ")[[1]][1])) {
      if(verbose == TRUE) {  cat("\t###### Detected mismatch in SNP numbers in common .tped and .frq files. Reducing to common SNPs\n") }
      comm10c = paste0("awk '{print $2",'" "',"}' comm",samp1,"____",samp2,".frq > comm",samp1,"_",samp2,"_SNPs3")
      system(comm10c)
      comm10d = paste0("plink --tfile ",samp1,"____",samp2," --extract comm",samp1,"_",samp2,"_SNPs3 --recode transpose --out ",samp1,"____",samp2," --allow-no-sex")
      system(comm10d, ignore.stdout = TRUE, ignore.stderr = TRUE)
      if(verbose == TRUE) {  cat(paste0("\t###### Number of common SNPs after reduction: ", as.integer(strsplit(system(paste0("wc -l comm",samp1,"_",samp2,"_SNPs3"), intern=T)," ")[[1]][1])),"\n") }
    }
    unlink(c(list.files(getwd(),paste0(samp1,"short\\.(bed|bim|fam|log|nosex)")),list.files(getwd(),paste0(samp2,"short\\.(bed|bim|fam|log|nosex)"))))
    unlink(list.files(getwd(),paste0(samp1,"____",samp2,"\\.(bed|bim|fam|log|nosex|tfam)")))
    unlink(list.files(getwd(),paste0("comm",samp1,"_",samp2,"_(SNPs|SNPs2|SNPs3)")))
    
    ## Read pairwise .TPED
    if(verbose == TRUE) {  cat(paste0("\t###### Trying to load ",samp1,"____",samp2,".tped\n")) }
    pairLoadNew = read.table(paste0(samp1,"____",samp2,".tped"),header=FALSE,stringsAsFactors = FALSE, sep=" ")
    if(ncol(pairLoadNew) != 8) {  cat(paste0("\t ERROR: Expected number of columns 8. Detected ",ncol(pairLoadNew),". Possible cause is duplicated IDs in .ped files\n"));  quit() }
    pairLoadNew = pairLoadNew[,-3]
    
    ## Read allele frequencies
    alFreqNew = read.csv(paste0("comm",samp1,"____",samp2,".frq"),stringsAsFactors=FALSE,sep="",colClasses = c("character","character","character","character","numeric","numeric"),header = F,col.names = c("CHR","SNP","A1","A2","MAF","NCHROBS"))
    if(verbose == TRUE) {  cat(paste0("\t###### Trying to load ",samp1,"____",samp2,".frq\n")) }
    alFreqNew$NCHROBS=NULL
    alFreqNew[5]=lapply(alFreqNew[5],round,5)
    alFreqNew$AFa2 = as.numeric(0)
    alFreqNew$AFa2=(1-alFreqNew$MAF)
    alFreqNew[6]=lapply(alFreqNew[6],round,5)
    alFrequNewSAFE = alFreqNew
    
    ## Add allele frequency data to main dataframe
    if((length(pairLoadNew$V1) == length(alFreqNew$CHR)) == TRUE) {
      pairLoadNew$A1Mi = alFreqNew$A1
      pairLoadNew$A2Ma = alFreqNew$A2
      pairLoadNew$al1freq = alFreqNew$MAF
      pairLoadNew$al2freq = alFreqNew$AFa2
    } 

    if (!all(pairLoadNew$V1 == alFreqNew$CHR)) {
      stop("Invalid chromosome sort order between allele frequencies and output tped file. Ensure your input dataset is numerically sorted across chromosomes and try again.")
    }

    ## Remove SNPs with fixed alleles
    pairLoadNew = pairLoadNew[pairLoadNew$al1freq != 0,]
    ## Remove SNPs with no data on .FRQ file
    if(length(which(is.na(pairLoadNew$al1freq)) != 0)) {  pairLoadNew = pairLoadNew[-(which(is.na(pairLoadNew$al1freq))),]  }
    if(verbose == TRUE) {  cat(paste0("\t###### Number of common SNPs after QC2: ",length(pairLoadNew$al1freq),"\n")) }
    
    ## Make sure all variants are SNPs
    ## Remove them from alFreq and commonDataFrame
    rem1 = which(nchar(alFreqNew$A1) >1)
    rem2 = which(nchar(alFreqNew$A2) >1)
    rem3 = c(rem1,rem2)
    remSnps = alFreqNew[rem3,2]
    if(length(rem1) >0) { alFreqNew = alFreqNew[-rem1,] }
    if(length(rem2) >0) { alFreqNew = alFreqNew[-rem2,] }
    rem10 = which(pairLoadNew$V2 %in% remSnps)
    if(length(rem10) >0){  pairLoadNew = pairLoadNew[-rem10,] }
    rem4 = which(nchar(pairLoadNew$V5) >1)
    rem5 = which(nchar(pairLoadNew$V6) >1)
    rem6 = which(nchar(pairLoadNew$V7) >1)
    rem7 = which(nchar(pairLoadNew$V8) >1)
    rem8 = c(rem4,rem5,rem6,rem7)
    if(length(rem8) >0){  pairLoadNew = pairLoadNew[-rem8,] }
    if(verbose == TRUE) {  cat(paste0("\t###### Number of common SNPs after QC3: ",length(pairLoadNew$al1freq),"\n")) }
    
    ## Test if both samples have all homozygous SNPs, and force so if FALSE
    library(data.table, quietly = T)
    if(all(pairLoadNew$V5 == pairLoadNew$V6) == FALSE) {
      setDT(pairLoadNew)
      pairLoadNew[, V5 := V5][runif(.N, 0, 1) < 0.5 , V5 := V6][]
      pairLoadNew[, V6 := V5][runif(.N, 0, 1) < 0.5 , V6 := V6][]
      setDF(pairLoadNew)
    }
    if(all(pairLoadNew$V7 == pairLoadNew$V8) == FALSE) {
      setDT(pairLoadNew)
      pairLoadNew[, V7 := V7][runif(.N, 0, 1) < 0.5 , V7 := V8][]
      pairLoadNew[, V8 := V7][runif(.N, 0, 1) < 0.5 , V8 := V8][]
      setDF(pairLoadNew)
    }
    
    ## Work with 5 decimal places for frequencies and make sure MAF+AFa2=1
    maf=format(as.numeric(pairLoadNew$al1freq),signif=5)
    pairLoadNew$al1freq=maf
    af2=format((1-as.numeric(maf)),signif=5)
    pairLoadNew$al2freq=af2
    
    if(ncol(pairLoadNew) != 11) {  cat(paste0("\t ERROR: Expected the following 11:\n"))
      cat("\t V1\tV2\tV4\tV5\tV6\tV7\tV8\tA1Mi\tA2Ma\tal1freq\tal2freq\n")
      cat(paste0("\t Detected the following ",ncol(pairLoadNew),":\n"))
      cat(paste0("\t ",colnames(pairLoadNew)),"\n")
      cat(paste0("\t Missing V5/V6 and V7/V8 can be due to problem loading genotypes from the 1st/2nd individual, respectively\n"))
      cat(paste0("\t Please check .ped and .tped files or enable verbose mode to try to identify possible errors in the data\n"))
      quit() }
    colnames(pairLoadNew) = c("chr","snp","pos","S1x","S1y","S2x","S2y","A1Mi","A2Ma","al1freq","al2freq")
    
    ## Remove triallelic SNPs
    toRemove=c()
    ## Check in individual 1
    it=1
    for(i in pairLoadNew$S1x) {
      if((i!=pairLoadNew$A1Mi[it])==TRUE && (i!=pairLoadNew$A2Ma[it])==TRUE) {
        toRemove=c(toRemove,pairLoadNew$snp[it])
      }
      it=it+1}
    ## Check in individual 2
    it=1
    for(i in pairLoadNew$S2x) {
      if((i!=pairLoadNew$A1Mi[it])==TRUE && (i!=pairLoadNew$A2Ma[it])==TRUE) {
        toRemove=c(toRemove,pairLoadNew$snp[it])
      }
      it=it+1}
    ## Remove those SNPs
    toRemove = unique(toRemove)
    poses = which(pairLoadNew$snp %in% toRemove)
    if(is.null(toRemove) == FALSE) {
      pairLoadNew=pairLoadNew[-poses,]
    }
    if(verbose == TRUE) {  cat("\t###### Final number of used SNPs after QC4: ",paste0(length(pairLoadNew$chr),"\t\n")) }
    if(verbose == FALSE) {  cat(paste0(length(pairLoadNew$chr)," SNPs","\t")) }
    
    toExport[line,] = c(samp1,samp2,length(pairLoadNew$snp),NA,NA,NA,NA)
    if(length(pairLoadNew$chr) >= ignoreThresh) {
      pairLoadNewRxy = pairLoadNew
      ## Remove duplicated allele
      pairLoadNewRxy = pairLoadNewRxy[,-5]
      pairLoadNewRxy = pairLoadNewRxy[,-6]
      ## Calculate components
      pairLoadNewRxy["PaXy"] = ifelse(as.character(pairLoadNewRxy$S1x) == pairLoadNewRxy$A1Mi, pairLoadNewRxy$al1freq, pairLoadNewRxy$al2freq)
      pairLoadNewRxy["PcXy"] = ifelse(as.character(pairLoadNewRxy$S2x) == pairLoadNewRxy$A1Mi, pairLoadNewRxy$al1freq, pairLoadNewRxy$al2freq)
      pairLoadNewRxy["IacADbcBD"] = ifelse(as.character(pairLoadNewRxy$S1x) == as.character(pairLoadNewRxy$S2x), 1, 0)*4
      pairLoadNewRxy$rxyL = ((0.5*(as.numeric(pairLoadNewRxy$IacADbcBD)))-(as.numeric(pairLoadNewRxy$PaXy)*2)) / (1+1-(as.numeric(pairLoadNewRxy$PaXy)*2))
      pairLoadNewRxy$ryxL = ((0.5*(as.numeric(pairLoadNewRxy$IacADbcBD)))-(as.numeric(pairLoadNewRxy$PcXy)*2)) / (1+1-(as.numeric(pairLoadNewRxy$PcXy)*2))
      rxyL = sum(pairLoadNewRxy$rxyL)
      ryxL = sum(pairLoadNewRxy$ryxL)
      pairLoadNewRxy$Av = (pairLoadNewRxy$rxyL + pairLoadNewRxy$ryxL) /2
      pairLoadNewRxy$SrQ1 = ((0.5*(as.numeric(pairLoadNewRxy$IacADbcBD)))-(as.numeric(pairLoadNewRxy$PaXy)*2))
      pairLoadNewRxy$SrQ2 = ((0.5*(as.numeric(pairLoadNewRxy$IacADbcBD)))-(as.numeric(pairLoadNewRxy$PcXy)*2))
      pairLoadNewRxy$SrQ1MinusSum = pairLoadNewRxy$SrQ1 - sum(pairLoadNewRxy$SrQ1)
      pairLoadNewRxy$SrQ2MinusSum = pairLoadNewRxy$SrQ2 - sum(pairLoadNewRxy$SrQ2)
      pairLoadNewRxy$WrQ1 = (1+1-(as.numeric(pairLoadNewRxy$PaXy)*2))
      pairLoadNewRxy$WrQ2 = (1+1-(as.numeric(pairLoadNewRxy$PcXy)*2))
      pairLoadNewRxy$WrQ1MinusSum = pairLoadNewRxy$WrQ1 - sum(pairLoadNewRxy$WrQ1)
      pairLoadNewRxy$WrQ2MinusSum = pairLoadNewRxy$WrQ2 - sum(pairLoadNewRxy$WrQ2)
      relatednessR = ((mean(pairLoadNewRxy$SrQ1MinusSum/pairLoadNewRxy$WrQ1MinusSum))+(mean(pairLoadNewRxy$SrQ2MinusSum/pairLoadNewRxy$WrQ2MinusSum)))/2
      
      if(verbose == TRUE) { cat(paste0("\t###### Estimated HRC and dyad number: ",round(relatednessR,3)," (",line,"/",length(combos[,1]),")\n"))}
      if(verbose == FALSE) { cat(paste0(round(relatednessR,3)," HRC\t(",line,"/",length(combos[,1]),")\n")) }
      
      if(is.na(relatednessR) == FALSE) {
        if(relatednessR < 0.0625) { degree = "Unrelated"}
        if(relatednessR >= 0.0625 && relatednessR < 0.1875) { degree = "2nd degree"}
        if(relatednessR >= 0.1875 && relatednessR < 0.3126) { degree = "1st degree"}
        if(relatednessR >= 0.3126) { degree = "Same individual/Twins"}
        ## Append to toExport table
        if(relatednessR == 1) {toExport[line,] = c(samp1,samp2,length(pairLoadNewRxy$snp),sprintf("%.4f", round(relatednessR,digits=4)),length(which(pairLoadNewRxy$IacADbcBD %in% 0)),"0",degree)
        } else {toExport[line,] = c(samp1,samp2,length(pairLoadNewRxy$snp),sprintf("%.4f", round(relatednessR,digits=4)),length(which(pairLoadNewRxy$IacADbcBD %in% 0)),length(which(pairLoadNewRxy$IacADbcBD %in% 4)),degree)}
      } else if(is.na(relatednessR) == TRUE) { toExport[line,] = c(samp1,samp2,length(pairLoadNewRxy$snp),"NA","0","0","NA") } 
    }
  } else if(file.size(paste0("comm",samp1,"_",samp2,"_SNPs")) == 0) {
    cat(paste0("SKIPPING due to no overlapping SNPs\t(",line,"/",length(combos[,1]),")\n"))
    toExport[line,] = c(samp1,samp2,"0","NA","0","0","NA")
  }  
}
## Export table with compilation of results
write.table(toExport,"TKGWV2_Results.txt",quote = F,row.names = F,col.names = T,sep="\t")

cat(paste0(" # [",Sys.time(),"] All dyads processed\n"))
cat(paste0(" # [",Sys.time(),"] Results exported to TKGWV2_Results.txt\n"))
