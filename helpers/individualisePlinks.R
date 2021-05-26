### Extract and convert all samples from binary PLINK dataset into individual text-PLINK
### Author - Daniel Fernandes
### Instructions:
### On line 22, you need to provide input for the compulsory argument 'dataset' when calling the function. Run the whole code afterwards.

individualisePlinks = function(dataset) {
  # Description of arguments #
  # dataset = Target dataset name without extensions
  
  dataLoad=read.table(paste0(dataset,".fam"), header=FALSE, sep="",as.is=TRUE, stringsAsFactors = FALSE)
  it=1
  for(i in dataLoad$V2) {
    to_write=paste(dataLoad[it,1],i)
    write.table(to_write,"to_keep",quote=FALSE,col.names = FALSE,row.names = FALSE)
    cmd1=paste0("plink --bfile ",dataset," --keep to_keep --geno 0.01 --keep-allele-order --allow-no-sex --recode --out ",i)
    system(cmd1)
    it=it+1
  }
  unlink("to_keep")
}

individualisePlinks(dataset = "v443_1240K_forKinship")
