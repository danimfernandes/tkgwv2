
### Downsample BAM files for faster analysis
downsampleN = 1500000    # Maximum number of reads

for(i in list.files(getwd(),pattern = "\\.bam$")) {   # Edit this line and next if you want to target only BAMS with a more specific extension/suffix than ".bam"
  sid = strsplit(i,"\\.bam")[[1]]
  comm0 = paste0("samtools idxstats ",i," | awk {'print $3'}")
  mapLen = sum(as.integer(system(comm0, intern=T)))
  rat = downsampleN/mapLen
  if(mapLen <= downsampleN) {
    rat2 = 0.999999
  } else if(mapLen > downsampleN) {
    rat2 = strsplit(as.character(rat),"\\.")[[1]][2]
  }
  comm1 = paste0("samtools view -s ",round(runif(1,1,10000),0),".",rat2," -b ",i," > ",sid,"_downsampled.bam")   # Edit this line if you want your downsampled BAMs with a different suffix other than "_downsampled.bam"
  system(comm1)
}

### Extract and convert all samples from binary PLINK dataset into individual text-PLINK
dataset="1240K_forKinship_TKGWV2"   # Dataset name without extensions
dataLoad=read.table(paste0(dataset,".fam"), header=FALSE, sep="",as.is=TRUE, stringsAsFactors = FALSE)

it=1
for(i in dataLoad$V2) {
  to_write=paste(dataLoad[it,1],i)
  write.table(to_write,"to_keep",quote=FALSE,col.names = FALSE,row.names = FALSE)
  cmd1=paste0("plink --bfile ",dataset," --keep to_keep --keep-allele-order --recode --out ",i)
  system(cmd1)
  it=it+1
}
unlink("to_keep")
