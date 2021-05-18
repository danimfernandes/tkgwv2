# TKGWV2
# An ancient DNA relatedness pipeline for ultra-low coverage whole genome shotgun data

## Citation
Fernandes, DM, Cheronet, O, Gelabert, P, Pinhasi, R. TKGWV2: An ancient DNA relatedness pipeline for ultra-low coverage whole genome shotgun data (2021). Preprint.







- PLINK frequency file was generated from 1000 Genomes Phase 3 CEU population. The only QC step required is to exclude fixed SNPs, reducing from 77 to 22 million variants
(maybe have default EUR files, where both the BED file and FRQ file have no fixed SNPs, will make things much faster)




You can access support files to run TKGWV2 from genome-wide BAM files and from 1240K data from here:
https://drive.google.com/drive/folders/1Aw-0v_7CUorHJOLpCJ0QVCwEdH43HlO4?usp=sharing

### Converting BAMs to individual text-PLINKs
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
#
# Input:
# - individual text-PLINK
# - pairwise FRQ
# Output:
# - pairwise relatedness coefficient result text file  
