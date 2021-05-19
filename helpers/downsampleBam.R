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

downsampleBAM(downsampleN = 1500000, extensionBam = "\\.bam$", extensionDownBam = "_subsampled.bam")
