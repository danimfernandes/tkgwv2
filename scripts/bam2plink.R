### Converting BAMs to individual text-PLINKs
### Author - Daniel Fernandes
# Programs with system-wide installation required # 
# 1. samtools
# 2. plink
# 3. Python3
#
# Description of arguments #
# referenceGenome = Path to index fasta file used to generate BAM files
# gwvList = List of positions to generate pileup file from in BED format
# gwvPlink = PLINK dataset with same positions as gwvList, to generate MAP files. A dummy dataset with a single individual is provided
# bamExtension = Default ".bam". Specific extension/suffix of the BAM files to be used. Anything before this will be considered as the sample ID
# minMQ = Default 30. Minimum mapping quality
# minBQ = Default 30. Minimum base quality
# excludeTerminalReadBases = Default FALSE. Exclude terminal read positions from pileup file
#
# Input:
# - individual BAM
# Output:
# - individual text-PLINK

args = commandArgs(trailingOnly=TRUE)
pywd = args[length(args)]
args = args[-length(args)]

cat("\n ################################################################################")
cat("\n ### TKGWV2 - An ancient DNA relatedness pipeline for ultra-low coverage data ###")
cat("\n ## Version 1.0 - Released 05/2021")
cat("\n #")
cat(paste0("\n # [",Sys.time(),"] ","Running 'bam2plink' on folder ", getwd()))
cat("\n\t # BAM >> Pileup >> Text-PLINK ")

error0 = c()
if(length(grep("referenceGenome",args)) != 0) {
  referenceGenome = strsplit(args[grep("referenceGenome",args)],"=")[[1]][2]
} else {
  error0 = c(error0,"referenceGenome")
}
if(length(grep("gwvList",args)) != 0) {
  gwvList = strsplit(args[grep("gwvList",args)],"=")[[1]][2]
} else {
  error0 = c(error0,"gwvList")
}
if(length(grep("gwvPlink",args)) != 0) {
  gwvPlink = strsplit(args[grep("gwvPlink",args)],"=")[[1]][2]
} else {
  error0 = c(error0,"gwvPlink")
}
if(length(error0)>0) {
  cat("\n\t # ERROR: The following required arguments could not be found:\n")
  cat(paste0("\t\t"),paste0(error0,collapse = ", "))
  cat("\n")
  quit()
}

bamExtension = ".bam"
if(length(grep("bamExtension",args)) != 0) { bamExtension = strsplit(args[grep("bamExtension",args)],"=")[[1]][2] }
minMQ = 30
if(length(grep("minMQ",args)) != 0) { minMQ = strsplit(args[grep("minMQ",args)],"=")[[1]][2] }
minBQ = 30
if(length(grep("minBQ",args)) != 0) { minBQ = strsplit(args[grep("minBQ",args)],"=")[[1]][2] }
excludeTerminalReadBases = FALSE
if(length(grep("excludeTerminalReadBases",args)) != 0) { excludeTerminalReadBases = TRUE }

cat("\n\t # Files to be processed:\n")
for(i in list.files(pattern = paste0(bamExtension,"$"))) { cat(paste0("\t\t ",i,"\n"))}
cat("\t # Arguments used:\n")
for(i in gsub("="," = ",args)) { cat(paste0("\t\t ",i,"\n"))}

error1 = c()
if(file.exists(referenceGenome) == F) {error1 = c(error1,referenceGenome)}
if(file.exists(gwvList) == F) {error1 = c(error1,gwvList)}
if(file.exists(paste0(gwvPlink,".bed")) == F) {error1 = c(error1,gwvPlink) }

if(length(error1)>0) {
  cat("\t # ERROR: The following required file(s) could not be found:\n")
  for(e in error1) { cat("\t\t ",e,"\n") }
} else {
  if(length(error1)==0) {
    for(i in list.files(pattern = paste0(bamExtension,"$"))) {
      sid = strsplit(i,bamExtension)[[1]]
      
      ## Generate pileup using samtools
      commA = paste0("samtools mpileup -Q ",minBQ," -q ",minMQ," -B -f ",referenceGenome," -l ",gwvList," ",i," > ",sid,".pileupsamtools.gwv.txt")
      system(commA, ignore.stderr=T)
      
      ## Convert sample positions to PLINK's range format ### INSTANTANEOUS
      comm2 = paste0("sed 's/chr//' ",sid,".pileupsamtools.gwv.txt | awk '{print $1,$2,$2,NR}' > ",sid,".plinkRange")
      system(comm2)
      
      ## Get map file using PLINK, and delete ped file: ### SLOW but still faster than awk
      comm3 = paste0("plink --bfile ",gwvPlink," --extract range ",sid,".plinkRange --recode --out ",sid," --allow-no-sex")
      system(comm3,ignore.stdout = TRUE)
      
      unlink(paste0(sid,".ped"))
      unlink(paste0(sid,".plinkRange"))
      unlink(paste0(sid,".log"))
      unlink(paste0(sid,".nosex"))
    }
    ## Convert pileup to text-PLINK
    if(excludeTerminalReadBases == FALSE) { excludeTerminalReadBases = "False"}
    if(excludeTerminalReadBases == TRUE) { excludeTerminalReadBases = "True"}
    if(file.exists(paste0(pywd,"/scripts/pileup2ped_updatedMay2021.py")) == F) {cat("\t # ERROR: File 'pileup2ped.py' could not be found in 'scripts' folder\n")} else {
      comm3b = paste0(pywd,"/scripts/pileup2ped_updatedMay2021.py ",excludeTerminalReadBases)
      system(comm3b)
    }
  }
  cat(paste0(" # [",Sys.time(),"] All BAM files processed\n"))
}

