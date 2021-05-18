############### Set of helper commands
### Downsize BAM files for faster analysis
subsampleN = 1500000      # Number of reads to downsample to
for(i in list.files(getwd(),pattern = "\\.bam$")) {
  sid = strsplit(i,"\\.bam")[[1]]
  comm0 = paste0("samtools idxstats ",i," | awk {'print $3'}")
  mapLen = sum(as.integer(system(comm0, intern=T)))
  rat = subsampleN/mapLen
  if(mapLen <= subsampleN) {
    rat2 = 0.999999
  } else if(mapLen > subsampleN) {
    rat2 = strsplit(as.character(rat),"\\.")[[1]][2]
  }
  comm1 = paste0("samtools view -s ",round(runif(1,1,10000),0),".",rat2," -b ",i," > ",sid,"_subsampled.bam")
  system(comm1)
}

### Extract and convert all samples from dataset into individual text-PLINK
dataset="v435_forKinship"
datasetLoad=read.table(paste0(dataset,".fam"), header=FALSE, sep="",as.is=TRUE, stringsAsFactors = FALSE)
it=1
for(i in datasetLoad$V2) {
  to_write=paste(datasetLoad[it,1],i)
  write.table(to_write,"to_keep",quote=FALSE,col.names = FALSE,row.names = FALSE)
  cmd1=paste0("plink --bfile ",dataset," --keep to_keep --recode --out ",i)
  system(cmd1)
  it=it+1
}
