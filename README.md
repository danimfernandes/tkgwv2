# TKGWV2
# An ancient DNA relatedness pipeline for ultra-low coverage whole genome shotgun data

# Description
TKGWV2 is a pipeline to estimate biological relatedness (1st, 2nd, and unrelated degrees) between individuals specifically aimed at ultra-low coverage ancient DNA data obtained from whole genome sequencing.
It is a massive update to the original TKrelated method published in 2017 (Fernandes et al.), with ease of use and efficiency as major concerns.
From the initial suggested minimum coverage threshold of 0.1X, TKGWV2 can be applied to BAM files of individuals with 0.025X, or, in some cases, as little as 0.005X when the other individual has above 0.03X. 

These characteristics have the potential to offer relatedness estimation during screening sequencing steps at early stages of an ancient DNA project, and can therefore be very useful for project planning.

Although mainly designed for shotgun data, TKGWV2 also works with, for example, the widely-used 1240K capture set, as long as adequate allele frequencies are provided.

# Requirements
TKGWV2 was developed for Linux. The following software needs to be available as system-wide installations:
- Python 3
- PLINK 1.9
- samtools (tested on version 1.7)
- R (tested on version 3+)
- R package: data.table (install.packages("data.table"))

As input, TKGWV2 can take either individual BAM or text-PLINK (ped/map) files .

# Installation and pipeline
Download and unzip this package. Keep the folder structure and make sure that all Python and R scripts are executable.

    wget https://github.com/danimfernandes/tkgwv2/archive/refs/heads/master.zip 
    unzip master.zip

The TKGWV2 pipeline is divided in two main parts:
- Generating text Pileup information from aligned BAM files and converting it to text-PLINK (ped/map) > **'bam2plink.R'**
- Generating pairwise transposed text-PLINK (tped) and allele frequency (frq) files > **'plink2tkrelated.R'**

This allows TKGWV2 to take as input either individual BAM or text-PLINK files.

![alt text](https://user-images.githubusercontent.com/22391172/118780982-e31a0e80-b88c-11eb-906a-655c04742859.png?raw=true)

Depending on what input data you use, TKGWV2 will require different support files (circled in diagram above).
### If starting from BAM files (orange point in diagram):
1. The reference genome (fasta) to which the files were aligned to;
2. A list of biallelic and non-fixed SNPs (bed) in 'bed' format, to generate Pileup information;
3. A binary PLINK dataset (bed/bim/fam) covering the exact same positions as the 'bed' file, to generate the text-PLINK 'map' files.

With these files, TKGWV2 will generate a set of text-PLINK files (ped/map) per individual. If you already have your data in PLINK format (e.g. from a 1240K dataset), you can convert it to individual text-PLINK sets and avoid running 'bam2plink.R'.

### If starting from individual text-PLINKs (ped/map) (green point in diagram):
1. A PLINK frequencies file (frq) containing the same SNPs (or a subset) of the ones in the ped/map files, obtained from a relevant large dataset.

We provide a sample set of these support files that can be downloaded from here:
https://drive.google.com/drive/folders/1Aw-0v_7CUorHJOLpCJ0QVCwEdH43HlO4?usp=sharing

Crucially, the allele frequencies provided must be representative of the populations the ancient individuals being tested belong to. Although the files we provide were generated using modern European CEU data from the 1000 Genomes Project Phase 3, we showed in the TKGWV2 publication that they work well for ancient populations from Europe spanning from the Neolithic to Medieval times. Due to this extensive geographical and temporal working window, you can, for example, use TKGWV2 on a pair of individuals from a Neolithic site in Poland at the same time that you use it for a Bronze Age pair from Spain.

For other untested regions and chronologies, or if you wish to use a different set of SNPs, it is very important to run a confirmation analysis on pairs of individuals with previously known relationships, to make sure that the used allele frequencies are not introducing any biases. Further below we provide a walkthrough for generating your own support files, in 'Generating own support files'.

# Usage, arguments, and examples
TKGWV2 includes R, Python, and Bash code. The two main utilities - 'bam2plink' and 'plink2tkrelated' - can be easily used and combined through the Python wrapper script TKGWV2.py.

    $ ./TKGWV2.py 
        
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
        - d, --dyads                      <path> Tab-spaced text file with each specified pair to be analysed per line, otherwise will run on every possible pair. Useful when analysing temporal or geographically distant pairs (as long as both are within the variation captured by the --freqFile)
        - i, --ignoreThresh               <int> Default 1. Threshold for the minimum number of SNPs allowed to estimate relatedness

### Example 1 - Starting from BAM files and running 'bam2plink' and then 'plink2tkrelated':

    $ ~/Software/tkgwv2-master/TKGWV2.py bam2plink --referenceGenome ~/Data/hg19/full_karyo.fa --gwvList ~/Software/tkgwv2-master/support/genomeWideVariants_hg19/1000GP3_22M_noFixed_noChr.bed --bamExtension final.bam --gwvPlink ~/Software/tkgwv2-master/support/genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed plink2tkrelated --freqFile ~/Software/tkgwv2-master/support/1000GP3_EUR_22M_noFixed.frq
    
### Example 2 - Starting from text-PLINKs (e.g. for 1240K data) and running 'plink2tkrelated':

    $ ~/Software/tkgwv2-master/TKGWV2.py plink2tkrelated --freqFile ~/Software/tkgwv2-master/support/1240K/1000GP3_EUR_1240K.frq
    
# Results and interpretation
All relationships will be written to a 7-column file named 'TKGWV2_Results.txt':

    Sample1	Sample2	Used_SNPs	HRC	counts0	counts4	Relationship
    Ind1	Ind2	17260	0.2248	1256	16004	1st degree
    Ind1	Ind3	18039	-0.0043	1701	16338	Unrelated
    Ind2	Ind3	16987	0.0569	1509	15478	Unrelated
   
Sample1 - First individual of tested pair<br/>
Sample2 - Second individual of tested pair<br/>
Used_SNPs - Used SNPs to calculate relatedness<br/>
HRC - Halved Relatedness Coefficient. Unrelated < 0.0625. 2nd Degree between 0.0625 and 0.1875. 1st Degree > 0.1875<br/>
counts0 - Number of non-shared alleles<br/>
counts4 - Number of shared alleles<br/>
Relationship - Descriptive relationship based on HRC value

TKGWV2 is only able to detect 1st and 2nd degree relationships due to their typical distribution ranges. The 3rd degree class, for example, by definition would have a mean HRC value of 0.0625, which would be located at the exact transition point between Unrelated and 2nd degree. With a theoretical distribution 0f 0.0625+-0.0625, the overlap with the 2 neighbouring classes would be substantial, and therefore would lead to uncertain classifications.<br/><br/>

![alt text](https://user-images.githubusercontent.com/22391172/118780776-ae0dbc00-b88c-11eb-9a7a-9552b909aa26.png?raw=true)

For these reasons, TKGWV2 does not include 3rd degree relationship estimation, and all real 3rd degree relatives will be assigned as 2nd degree or Unrelated, depending on whether their HRC value is slightly above or below 0.0625, respectively.

However, as part of a provided helper script to produce simulated ranges like the ones above ('distSimulations.R'), you have the option to include 3rd degree curves and calculate the probability of an ancient pair's HRC belonging to that class and the overlapping ones. *Do be aware that this must only be used as suggestive evidence for investigation purposes, and never as definite proof of 3rd degree relationships.*

# Helper scripts
We provide 4 R functions to help automate some situations you might come across while preparing your data for TKGWV2 or analysing the results.
- 'downsampleBAM.R' can be used to downsample BAM files for faster processing. Check the 'Tips and suggestions' section below for more details on this.
- 'individualisePlink.R' can be used to extract and convert all samples from a binary PLINK dataset into individual text-PLINK sets (ped/map). This is useful, for example, if you want to use TKGWV2 on 1240K data. After running this function you can run TKGWV2 directly from 'plink2tkrelated'.
- 'distSimulations.R' can be used to generate simulated distribution curves and posterior probabilities for estimated HRC values based on an input plink frequencies file (frq).
- 'simsForRanges.R' can be used to generate simulated distribution curves over a set of increasing SNP numbers in order to assess error rates and SNP thresholds for an input plink frequencies file (frq).





# Generation of own support files
- PLINK frequency file was generated from 1000 Genomes Phase 3 CEU population. The only QC step required is to exclude fixed SNPs, reducing from 77 to 22 million variants
(maybe have default EUR files, where both the BED file and FRQ file have no fixed SNPs, will make things much faster)
- ................................................






# Tips and suggestions
- *Downsample your data*<br/>
In the TKGWV2 publication, we showed that from 15000 used SNPs the error rates were under 0.50%, therefore one of our main suggestions is to downsample your data for a more time-efficient analyses. When starting from BAM files, we showed that downsampling them to a maximum of 1.3-1.6 million reads was sufficient to obtain an average of 23000 used SNPs per pair, and consequently run TKGWV2 with very low error rates at a much higher speed than if larger BAM files were to be used.

- *Use simulations to get posterior probabilities for <15000 or <10000 SNPs*<br/>
Although at 10000 SNPs the estimated error rates when using the provided support files are at around 3%, we provide a helper script to run simulations on the specific set of SNPs used to estimate that pair's relatedness. This information is contained in the files named as 'commInd1_Ind2.frq'.<br/><br/>
The simulations will generate distribution ranges for the 3 relatedness classes, and use them them to calculate the probability of the pair's estimate to represent each class.<br/><br/>
The R script for this, 'distSimulations.R', can be found in the folder 'helpers'.<br/><br/>
# EDIT THIIIIIIS WITH NEW 1240K RESULTS
The 10-15000 SNPs thresold was identified for the genome wide data set based on shotgun data, which likely includes less-informative variants than, for example, the curated 1240K SNP set. Consequently, for the latter, the minimum number of SNPs that would correspond to similar error rates should be substantially lower. On an updated version of TKGWV2 we will include information on error rates and thresholds for the 1240K dataset.

- *Always run confirmation analysis for previously untested datasets and frequencies*<br/>
As mentioned above, it is important that the allele frequencies provided properly reflect the genetic composition of the ancient population the individuals being tested belonged to. Whether you use closely related modern populations or a large set of ancient individuals to obtain allele frequencies from, make sure to confirm that those frequencies are appropriate by testing the pipeline on previously published individuals from the same (or minimally closely related) population with known relationships.<br/><br/>
This is an important caveat that needs to be considered, and might prove challenging for regions or periods for which there is a lack of published ancient relatives. <br/><br/>
However, with the ongoing exponential increase in availability of both modern and aDNA data around the world, TKGWV2 can potentially be applied to the great majority of situations.<br/>

   
   
   
   
   




# Citation and contact
Fernandes, DM, Cheronet, O, Gelabert, P, Pinhasi, R. TKGWV2: An ancient DNA relatedness pipeline for ultra-low coverage whole genome shotgun data (2021). bioRxiv.

Feel free to contact Daniel Fernandes for questions and/or suggestions - dani.mag.fernandes(at)gmail.com
