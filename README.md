# TKGWV2
# An ancient DNA relatedness pipeline for ultra-low coverage whole genome shotgun data

# Description
TKGWV2 is a pipeline to estimate biological relatedness (1st, 2nd, and unrelated degrees) between individuals specifically aimed at ultra-low coverage ancient DNA data obtained from whole genome sequencing.
It is a massive update to the original TKrelated method published in 2017 (Fernandes et al.), with ease of use and efficiency as major concerns.
From the initial suggested minimum coverage threshold of 0.1X, TKGWV2 can be applied to BAM files of individuals with 0.025X, or, in some cases, as little as 0.005X when the other individual has above 0.03X. These characteristics have the potential to offer relatedness estimation during screening sequencing steps at early stages of an ancient DNA project, and can therefore be very useful for project planning.

Although mainly designed for shotgun data, TKGWV2 also works with, for example, the widely-used 1240K capture set, as long as allele frequencies are provided.

# Requirements
TKGWV2 was developed for Linux. The user is required to have the following software available as a system-wide installation:
- Python 3
- PLINK 1.9
- samtools (tested on version 1.7)
- R (tested on version 3+)
- R package: data.table (install.packages("data.table"))

As input, TKGWV2 can take either individual BAM or text-PLINK (ped/map) files .

# Installation and pipeline
Download and unzip this package. Keep the folder structure and make sure that all Python and R scripts are executable.
The TKGWV2 pipeline is divided in two main parts:
- Generating text Pileup information from aligned BAM files and converting it to text-PLINK (ped/map) > **'bam2plink.R'**
- Generating pairwise transposed text-PLINK (tped) and allele frequency (frq) files > **'plink2tkrelated.R'**

This allows TKGWV2 to take as input either individual BAM or text-PLINK files.

![alt text](https://github.com/danimfernandes/tkgwv2/blob/master/tkgwv2_pipeline.png?raw=true)

Depending on what input data you use, TKGWV2 will require different support files.
### If starting from BAM files (orange point in diagram above):
1. The index genome (fasta) to which the files where aligned to;
2. A list of biallelic and non-fixed SNPs (bed) in 'bed' format, to generate Pileup information;
3. A binary PLINK dataset (bed/bim/fam) covering the exact same positions as the 'bed' file, to generate the text-PLINK 'map' files.

With these files, TKGWV2 will generate a set of text-PLINK files (ped/map) per individual. If you already have your data in PLINK format (e.g. from a 1240K dataset), you can convert it to individual text-PLINK sets and avoid running 'bam2plink.R'.

### If starting from individual text-PLINKs (ped/map):
1. A PLINK frequencies file (frq) containing the same SNPs (or a subset) of the ones in the ped/map files, obtained from a relevant large dataset.

We provide a sample set of these support files that can be downloaded from here:
https://drive.google.com/drive/folders/1Aw-0v_7CUorHJOLpCJ0QVCwEdH43HlO4?usp=sharing

Crucially, the allele frequencies must be representative of the ancient individuals being tested. Although the files we provide were generated using modern European CEU data from the 1000 Genomes Project Phase 3, we showed in the original publication of TKGWV2 that they work well for ancient populations from Europe spanning from the Neolithic to Medieval times.

For other untested regions and chronologies, or if the user wishes to use a different set of SNPs, it is very important to run a confirmation analysis on pairs of individuals with previously known relationships, to make sure that the used allele frequencies are not introducing substantial bias. Further below, in 'Generation of own support files', we provide a walkthrough for these situations.

# Usage, arguments, and examples
TKGWV2 includes R, Python, and Bash code. The two main utilities - 'bam2plink' and 'plink2tkrelated' can be easily used and combined through a Python wrapper script - TKGWV2.py.

    $ ./TKGWV2.py 
    
    Main analyses: ./ngsrelate  [options] 
    
    Utilities and options:
      bam2plink  
        [mandatory arguments]
        - r, --referenceGenome            <path> Path to index fasta file used to generate BAM files
        - L, --gwvList                    <path> List of positions to generate pileup file from in BED format
        - p, --gwvPlink                   <path> Path to PLINK dataset with same positions as --gwvList, to generate MAP files. A dummy dataset with a single individual is provided
        [optional arguments]
        - e, --bamExtension               <str> Default ".bam". Specific extension/suffix of the BAM files to be used. Anything before this will be considered as the sample ID
        - m, --minMQ                      <int> Default 30. Minimum mapping quality
        - b, --minBQ                      <int> Default 30. Minimum base quality
        - t, --excludeTerminalReadBases   Exclude terminal read positions from pileup file, otherwise include by default
        - h, --help                       Displays this help

      plink2tkrelated
        [mandatory arguments]
        - f, --freqFile                   <path> Allele frequencies file in PLINK format (FRQ) containing the same SNPs (or a subset) of the ones in the input PED/MAP files
        [optional arguments]
        - d, --dyads                      <path> Tab-spaced text file with each pair to be analysed per line, otherwise will run on every possible pair
        - i, --ignoreThresh               <int> Default 1. Threshold for the minimum number of SNPs allowed to estimate relatedness.

### Example 1 - Starting from BAM files and running 'bam2plink' and then 'plink2tkrelated':

    $ ~/Software/tkgwv2-master/TKGWV2.py bam2plink --referenceGenome ~/Data/hg19/full_karyo.fa --gwvList ~/Software/tkgwv2-master/support/1000GP3_22M_noFixed_noChr.bed --bamExtension final.bam --gwvPlink ~/Software/tkgwv2-master/support/DummyDataset_EUR_22M_noFixed plink2tkrelated --freqFile ~/Software/tkgwv2-master/support/1000GP3_EUR_22M_noFixed.frq
    
### Example 2 - Starting from text-PLINKs (e.g. for 1240K data) and running 'plink2tkrelated':

    $ ~/Software/tkgwv2-master/TKGWV2.py plink2tkrelated --freqFile ~/Software/tkgwv2-master/support/1000GP3_EUR_1240K.frq
    
# Results output
All relationships will be written to a 7-column file 'TKGWV2_Results.txt':

    Sample1	Sample2	Used_SNPs	HRC	counts0	counts4	Relationship
    Ind1	Ind2	17260	0.2248	1256	16004	1st degree
    Ind1	Ind3	18039	-0.0043	1701	16338	Unrelated
    Ind2	Ind3	16987	0.0569	1509	15478	Unrelated
   
Sample1 - First individual of tested pair
Sample2 - Second individual of tested pair
Used_SNPs - Used SNPs to calculate relatedness
HRC - Halved Relatedness Coefficient. Unrelated < 0.0625. 2nd Degree between 0.0625 and 0.1875. 1st Degree > 0.1875 
counts0 - Number of non-shared SNPs
counts4 - Number of shared SNPs
Relationship - Descriptive relationship based on HRC value
   
   
   
   
   
   
# Generation of own support files






- PLINK frequency file was generated from 1000 Genomes Phase 3 CEU population. The only QC step required is to exclude fixed SNPs, reducing from 77 to 22 million variants
(maybe have default EUR files, where both the BED file and FRQ file have no fixed SNPs, will make things much faster)






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




## Citation
Fernandes, DM, Cheronet, O, Gelabert, P, Pinhasi, R. TKGWV2: An ancient DNA relatedness pipeline for ultra-low coverage whole genome shotgun data (2021). Preprint.
