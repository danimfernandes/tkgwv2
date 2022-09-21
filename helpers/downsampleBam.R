### Downsize BAM files for faster analysis
### Author - Daniel Fernandes
### Instructions:
### There are no mandatory arguments that need input information, but if you want to edit any, do so on line 31, and run the whole code afterwards.

downsampleBam = function(downsampleN = 1500000, extensionBam = "\\.bam$", suffixDownBam = "_subsampled") {
  # Description of arguments #
  # downsampleN = Default 1500000. Maximum number of reads to downsample BAM files to
  # extensionBam = Default "\\.bam$". Work on BAM files with this extension/suffix
  # suffixDownBam = Default "_subsampled". Suffix for new downsampled BAM files
  #
  # Make sure your BAMs are indexed
  
  for(i in list.files(getwd(),pattern = extensionBam)) {
    sid = strsplit(i,paste0(extensionBam,"$"))[[1]]
    comm0 = paste0("samtools idxstats ",i," | awk {'print $3'}")
    mapLen = sum(as.integer(system(comm0, intern=T)))
    rat = downsampleN/mapLen
    if(mapLen <= downsampleN) {
      rat2 = 999999
    } else if(mapLen > downsampleN) {
      rat2 = strsplit(as.character(rat),"\\.")[[1]][2]
    }
    comm1 = paste0("samtools view -s ",round(runif(1,1,10000),0),".",rat2," -b ",i," > ",sid,suffixDownBam,".bam")
    system(comm1)
  }
}

# Call the function if this script is at the root of the stack (i.e. not called as a source file)
if (sys.nframe() == 0) {
  downsampleBam(downsampleN = 1500000, extensionBam = "\\.bam$", suffixDownBam = "_subsampled")
}
