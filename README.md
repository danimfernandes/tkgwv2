# TKGWV2
# An ancient DNA relatedness pipeline for ultra-low coverage whole genome shotgun data
# Description
TKGWV2 is a pipeline to estimate biological relatedness (1st, 2nd, and unrelated degrees) between individuals specifically aimed at ultra-low coverage ancient DNA data obtained from whole genome sequencing.
It is a massive update to the original TKrelated method published in 2017 (Fernandes et al.), with ease of use and efficiency as major concerns.
From the initial suggested minimum coverage threshold of 0.1X, TKGWV2 can be applied to BAM files of individuals with 0.025X, or, in some cases, as little as 0.005X when the other individual has above 0.03X. 

These characteristics have the potential to offer relatedness estimation during screening sequencing steps at early stages of an ancient DNA project, and can therefore be very useful for project planning.

Although mainly designed and thoroughly tested for shotgun data, TKGWV2 can also work with, for example, the widely-used 1240K capture set, as long as adequate allele frequencies are provided and some precautions are taken (see notes at the end of this document for more details).

# Requirements
### Version 1.0a - Released 06/2021
TKGWV2 was developed for Linux. The following software needs to be available as system-wide installations:
- Python 3
- PLINK 1.9
- samtools (tested on version 1.7)
- R (tested on version 3+)
- R package: data.table (install.packages("data.table"))

As input, TKGWV2 can take either individual BAM or text-PLINK (ped/map) files .

# Installation and pipeline
Download and unzip this package. Keep the folder structure and make sure that all Python and R scripts have permissions and are executable.

    wget https://github.com/danimfernandes/tkgwv2/archive/refs/heads/master.zip 
    unzip master.zip

The TKGWV2 pipeline is divided in two main parts:
- Generating text Pileup information from aligned BAM files and converting it to text-PLINK (ped/map) > **'bam2plink.R'**
- Generating pairwise transposed text-PLINK (tped) and allele frequency (frq) files > **'plink2tkrelated.R'**

This allows TKGWV2 to take as input either individual BAM or text-PLINK files.

![alt text](https://user-images.githubusercontent.com/22391172/118780982-e31a0e80-b88c-11eb-906a-655c04742859.png?raw=true)

Depending on what input data you use, TKGWV2 will require different support files (circled in diagram above).
***
### If starting from BAM files (orange point in diagram):
1. The reference genome (fasta) to which the files were aligned to;
2. A list of biallelic and non-fixed SNPs (bed) in 'bed' format (chr pos pos), to generate Pileup information. Not to be confused with the PLINK dataset .bed file;
3. A binary PLINK dataset (bed/bim/fam) covering the exact same positions as the 'bed' file from point 2, to generate the text-PLINK 'map' files.

With these files, TKGWV2 will generate a set of text-PLINK files (ped/map) per individual. If you already have your data in PLINK format (e.g. from a 1240K dataset), you can convert it to individual text-PLINK sets and avoid running 'bam2plink.R'.
***
### If starting from individual text-PLINKs (ped/map) (green point in diagram):
1. A PLINK frequencies file (frq) containing the same SNPs (or a subset) of the ones in the ped/map files, obtained from a relevant large dataset.
***
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

<br/>

***
### Example 1 - Starting from BAM files and running 'bam2plink' and then 'plink2tkrelated':

    $ ~/Software/tkgwv2-master/TKGWV2.py bam2plink --referenceGenome ~/Data/hg19/full_karyo.fa --gwvList ~/Software/tkgwv2-master/support/genomeWideVariants_hg19/1000GP3_22M_noFixed_noChr.bed --bamExtension final.bam --gwvPlink ~/Software/tkgwv2-master/support/genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed plink2tkrelated --freqFile ~/Software/tkgwv2-master/support/1000GP3_EUR_22M_noFixed.frq

### Example 1.1 - Starting from BAM files, downsampling, and running 'bam2plink' and then 'plink2tkrelated':

Use the 'downsampleBam.R' helper script:

    ### Downsize BAM files for faster analysis
    downsampleBAM = function(downsampleN = 1500000, extensionBam = "\\.bam$", suffixDownBam = "_subsampled") {}
    downsampleBAM(downsampleN = 1500000, extensionBam = "\\.bam$", suffixDownBam = "_subsampled")
 
 Now you can use the downsampled BAM files with suffix "\_subsampled":

    $ ~/Software/tkgwv2-master/TKGWV2.py bam2plink --referenceGenome ~/Data/hg19/full_karyo.fa --gwvList ~/Software/tkgwv2-master/support/genomeWideVariants_hg19/1000GP3_22M_noFixed_noChr.bed --bamExtension _subsampled.bam --gwvPlink ~/Software/tkgwv2-master/support/genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed plink2tkrelated --freqFile ~/Software/tkgwv2-master/support/1000GP3_EUR_22M_noFixed.frq

<br/>

***
### Example 2 - Starting from text-PLINKs (e.g. for 1240K data) and running 'plink2tkrelated':

    $ ~/Software/tkgwv2-master/TKGWV2.py plink2tkrelated --freqFile ~/Software/tkgwv2-master/support/1240K/1000GP3_EUR_1240K.frq

### Example 2.1 - Starting from 1240K dataset, individualising text-PLINKs, downsampling, and running 'plink2tkrelated':
Use the 'individualisePlinks.R' helper script:
    
    ### Extract and convert all samples from binary PLINK dataset into individual text-PLINK
    individualisePlinks = function(dataset) {}
    individualisePlinks(dataset = "v443_1240K_forKinship")
    
Use the 'downsamplePed.R' helper script:

    ### Downsize text-PLINK (ped/map) for faster analysis
    downsamplePED = function(downsampleN = 80000, extensionPed = "\\.ped$", suffixDownPed = "_subsampled") {}
    downsamplePED(downsampleN = 80000, extensionPed = "\\.ped$", suffixDownPed = "_subsampled")

Now either generate a text file listing the downsampled file pairs and use that as the argument for the 'dyads' option of 'plink2tkrelated', or move/delete the original PED/MAP files, as by default 'plink2tkrelated' will work on all PED/MAP file pairs in the working folder:

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

***Disclaimer:*** For these reasons, TKGWV2 does not include 3rd degree relationship estimation, and all real 3rd degree relatives will be assigned as 2nd degree or Unrelated, depending on whether their HRC value is above or below 0.0625, respectively. In the TKGWV2 manuscript we showed that 3rd degree relatives tend to have an HRC between 0.0625 and ~0.090, so you might want to consider referring to any result within this interval as “2nd- or 3rd-degree” to cover this possibility. This interval was obtained based on only 11 3rd-degree relatives, therefore is only informative. Furthermore, it relies on the their correct assessment per the original publication.

# Helper scripts
We provide 4 R functions to help automate some situations you might come across while preparing your data for TKGWV2 or analysing the results.
- 'downsampleBam.R' can be used to downsample BAM files for faster processing. Check the 'Tips and suggestions' section below for more details on this.
- 'downsamplePed.R' can be used to downsample text-PLINK files for faster processing. Check the 'Tips and suggestions' section below for more details on this.
- 'individualisePlink.R' can be used to extract and convert all samples from a binary PLINK dataset into individual text-PLINK sets (ped/map). This is useful, for example, if you want to use TKGWV2 on 1240K data. After running this function you can run TKGWV2 directly from 'plink2tkrelated'.
- 'distSimulations.R' can be used to generate simulated distribution curves and posterior probabilities for estimated HRC values based on an input plink frequencies file (frq).
- 'simsForErrorRates.R' can be used to generate simulated distribution curves over a set of increasing SNP numbers in order to assess error rates and SNP thresholds for an input plink frequencies file (frq). The output is an error rates vs SNPs plot and a text file with the data used to generate it. A gray vertical line is drawn at the lowest number of SNPs tested with an average error <=1%. The text file is a concatenation of the tables of values for 1st degree, 2nd degree, and Unrelated, with the SNP number as the first column, and the general header as the number of simulations. 


# Description and generation of support files
For the genome-wide application of TKGWV2 a few support files are required, as described above. For the 1240K application, only a PLINK frequency file (FRQ) from the reference population(s) of interest is required.

- *Genome-wide support files*<br/>
To generate Pileups with _samtools mpileup_, a list of <ins>biallelic</ins>, <ins>non-fixed</ins>, and <ins>unique SNPs</ins>, in 'bed' format (chr pos-1 pos) is required:
```
1 10538 10539
1 11007 11008
1 11011 11012
1 13109 13110
1 13115 13116
(...)
```
These positions can initially be obtained from the 1000 Genomes VCF files (e.g. http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/), using _bcftools_, and QCed down the line to exclude unwanted variants. Be sure to avoid fixed SNPs, SNPs with different names but same genomic positions, insertions, deletions, etc.:
```
#Note: Dummy coordinates and SNPids in PLINK's 'bim' file format
1       rs000000000     0       11111   A       C
1       ss111111111     0       11111   G       C

3       rs222222222     0       22222   G       GAAA
8       rs333333333     0       33333   AATTAG       AATT
14       rs444444444;rs55555555     0       44444   A       G
```

Then, to help generate the individual 'map' files, a binary PLINK dataset covering the <ins>exact</ins> same positions as the 'bed' file above is needed. PLINK's _--vcf_ can be used to convert from VCF.GZ to binary PLINK. For increased efficiency and speed, we suggest to generate this dataset with a single individual.

'bim' file:
```
1       rs537182016     0       10539   A       C
1       rs575272151     0       11008   G       C
1       rs544419019     0       11012   G       C
1       rs540538026     0       13110   A       G
1       rs62635286      0       13116   G       T
(...)
```

'fam' file:
```
HG00096 HG00096 0 0 0 -9
```
Lastly, after isolating the individuals from the relevant population to source the allele frequencies from - again for the <ins>exact</ins> same positions as in the other files - the FRQ file can be obtained using PLINK's _--frq_ argument:
```
1   rs537182016    A    C     0.000994     1006
1   rs575272151    G    C      0.08847     1006
1   rs544419019    G    C      0.08847     1006
1   rs540538026    A    G      0.05666     1006
1    rs62635286    G    T       0.1869     1006
```

- *1240K support files*<br/>
The required FRQ file is likely to be easily retrievable with PLINK's _--frq_ argument from available datasets such as the Allen Ancient DNA Resource (AADR), curated by the Reich Lab in Harvard:<br/>
https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data.

 
# Tips, suggestions, and notes
- *Downsample your data*<br/>
In the TKGWV2 publication we showed that, with the genome-wide SNP set, from 15000 used SNPs the error rates were under 0.50%, therefore one of our main suggestions is to downsample your data for a more time-efficient analyses. When starting from BAM files, we showed that downsampling them to a maximum of 1.3-1.6 million reads per file was sufficient to obtain an average of 23000 used SNPs per pair, and consequently run TKGWV2 with very low error rates at a much higher speed than if larger BAM files were to be used.<br/><br/>
Similarly, for 1240K data, 60-100K SNPs per individual produced pairwise estimates with average used SNPs between 2400-6700, and error rates between ~3 and <0.5%, respectively (based on simulated data). Doing this on individuals with ~800K SNPs reduced the running time of 'plink2tkrelated' by 36 times, from 330 to 9 seconds, for 10 relationships.<br/><br/>
We provide 2 helper scripts to downsample BAM files and text-PLINK sets, so you can downsample whole-genome ('downsampleBam.R') and 1240K data ('downsamplePed.R'). Their default subsampling numbers are 1.5 million reads for BAM files, and 80K SNPs for 1240K text-PLINKs.

- *Use simulations to get posterior probabilities for <15000 or <10000 SNPs*<br/>
The estimated error rates when using the provided support files with genome-wide SNP set are at 3% for 10000 SNPs and 0.5% at 15000 SNPs, so we provide a helper script in order to see if these potentially less accurate estimates overlap with more than one relatedness class. This helper script runs simulations on the specific set of SNPs used to estimate that pair's relatedness. This information is contained in the files named as 'commInd1_Ind2.frq'.<br/><br/>
The simulations will generate distribution ranges for the 3 relatedness classes, and use them them to calculate the probability of the pair's estimate to represent each class.<br/><br/>
The R script for this, 'distSimulations.R', can be found in the folder 'helpers'.<br/><br/>
*Note:* The 10-15000 SNPs thresold was identified for the provided genome-wide SNP set based on shotgun data, which likely includes less-informative variants than, for example, the curated 1240K SNP set. Consequently, for the latter, the minimum number of SNPs that would correspond to similar error rates of <1% are substantially lower at around 3000 SNPs instead of between 10-15000. However, the use of a curated SNP set might introduce small HRC biases that the user needs to take into consideration, especially for HRC estimates that are too close to the outer edges of a class's distribution range.

- *Always run confirmation analysis and calculate error rates for previously untested datasets and frequencies*<br/>
As mentioned above, it is important that the allele frequencies provided properly reflect the genetic composition of the ancient population the individuals being tested belonged to. Whether you use closely related modern populations or a large set of ancient individuals to obtain allele frequencies from, make sure to confirm that those frequencies are appropriate by testing the pipeline on previously published individuals from the same (or minimally closely related) population with known relationships.<br/><br/>
This is an important caveat that needs to be considered, and might prove challenging for regions or periods for which there is a lack of published ancient relatives. <br/><br/>
However, with the ongoing exponential increase in availability of both modern and aDNA data around the world, TKGWV2 can potentially be applied to the great majority of situations.<br/>

- *For when imperfect allele frequencies are used, try a normalization step*<br/>
If the alelle frequencies used are not perfect or seem biased (based on comparisons with previously known relatives or based on a non-zero average of unrelated HRCs), it might be possible to correct the results by applying a normalization step like in READ (Monroy-Kuhn et al. 2018). As the theoretical average HRC for a set of unrelated individuals is supposed to be 0, the user can find the average HRC of all unrelated pairs and use that value (i.e. its deviation from 0) to correct the HRC estimates. 
In most situations the unrelated average HRC is usually above 0, and suggest incorrect or biased allele frequencies, due to a lack of very low frequency (but not fixed) SNPs, or simply substantially different frequencies. Examples of these situations are the use of modern European allele frequencies on ancient hunter-gatherer individuals from Europe, as well as when using curated SNP sets such as the 1240K dataset (see next note).

- *If you are using data restricted to a curated dataset, such as the 1240K SNP set*<br/>
A curated SNP set such as the 1240K includes ancestry SNPs selected based on the existing variability of modern populations. Therefore, issues and biases such as the one mentioned in the previous note are more likely to happen in these capture sets, than when using whole-genome shotgun data. The introduced bias is likely small (HRC inflations of 0.02-0.03 have been observed), but can nevertheless produce false positives for HRC estimates that are too close to the other classes' distribution ranges. Therefore, providing allele frequencies from genetically close ancient population for these situations might be preferable, when available.




# Citation and contact
Fernandes, DM, Cheronet, O, Gelabert, P, Pinhasi, R. TKGWV2: an ancient DNA relatedness pipeline for ultra-low coverage whole genome shotgun data. Sci Rep 11, 21262 (2021). https://doi.org/10.1038/s41598-021-00581-3

Feel free to contact Daniel Fernandes for questions and/or suggestions - dani.mag.fernandes(at)gmail.com

# Changelog

04-02-2022
- Fixed error in "plink2tkrelated.py" for dyads with 0 common SNPs
- Small QOL changes to "plink2tkrelated.py" stdout, which is now also printing HRC values and skipped dyads

27-06-2021
- Corrected missing " and wrong link to "pileup2ped.py" in "bam2plink.R"
- Removed duplicated SNP from support WGV files
