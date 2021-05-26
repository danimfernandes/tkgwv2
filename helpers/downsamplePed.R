### Downsize text-PLINK (ped/map) for faster analysis
### Author - Daniel Fernandes

downsamplePED = function(downsampleN = 80000, extensionPed = "\\.ped$", suffixDownPed = "_subsampled") {
  # Description of arguments #
  # downsampleN = Default 80000. Maximum number of SNPs to downsample text-PLINK dataset to
  # extensionPed = Default "\\.ped$". Work on PED/MAP files with this extension/suffix
  # suffixDownPed = Default "_subsampled". Extenstion/suffix for new downsampled files
  
  for(i in list.files(getwd(),pattern = extensionPed)) {
    sid = strsplit(i,paste0(extensionPed,"$"))[[1]]
    comm0 = paste0("wc -l ",sid,".map")
    mapLen = as.integer(strsplit(system(comm0, intern=T), " ")[[1]][1])
    #rat = downsampleN/mapLen
    if(mapLen <= downsampleN) {
      rat2 = mapLen
    } else if(mapLen > downsampleN) {
      rat2 = downsampleN
    }
    comm1 = paste0("plink --file ",sid," --thin-count ",rat2," --recode --out ",sid,suffixDownPed," --allow-no-sex --keep-allele-order")
    system(comm1)
  }
}

downsamplePED(downsampleN = 80000, extensionPed = "\\.ped$", suffixDownPed = "_subsampled")
