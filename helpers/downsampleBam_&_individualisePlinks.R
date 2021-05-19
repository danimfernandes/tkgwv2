### Downsize BAM files for faster analysis
downsampleBAM = function(downsampleN = 1500000, extensionBam = "\\.bam$", extensionDownBam = "_subsampled.bam") {
  # Description of arguments #
  # downsampleN = Default 1500000. Maximum number of reads to downsample BAM files to
  # extensionBam = Default "\\.bam". Work on BAM files with this extension/suffix
  # extensionDownBam = Default "_subsampled.bam". Extenstion/suffix for new downsampled BAM files
  #
  # Make sure your BAMs are indexed
  
  for(i in list.files(getwd(),pattern = extensionBam)) {
    sid = strsplit(i,paste0(extensionBam,"$"))[[1]]
    comm0 = paste0("samtools idxstats ",i," | awk {'print $3'}")
    mapLen = sum(as.integer(system(comm0, intern=T)))
    rat = downsampleN/mapLen
    if(mapLen <= downsampleN) {
      rat2 = 0.999999
    } else if(mapLen > downsampleN) {
      rat2 = strsplit(as.character(rat),"\\.")[[1]][2]
    }
    comm1 = paste0("samtools view -s ",round(runif(1,1,10000),0),".",rat2," -b ",i," > ",sid,extensionDownBam)
    system(comm1)
  }
}


### Extract and convert all samples from binary PLINK dataset into individual text-PLINK
individualisePlinks = function(dataset) {
  # Description of arguments #
  # dataset = Target dataset name without extensions
  
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
}

# Edit these lines and (un)comment to run the wanted functions
downsampleBAM(downsampleN = 1500000, extensionBam = "\\.bam$", extensionDownBam = "_subsampled.bam")
individualisePlinks(dataset = "v443_1240K_forKinship")
