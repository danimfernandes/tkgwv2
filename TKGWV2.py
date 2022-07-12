#!/usr/bin/python3

__author__ = 'Daniel Fernandes'

import os
import sys
import re
import getopt

cwd = os.getcwd()
pywd = os.path.dirname(os.path.realpath(__file__))

argList = sys.argv[1:]
if "bam2plink" in argList and "plink2tkrelated" in argList:
	#input()
	if argList.index("bam2plink") < argList.index("plink2tkrelated"):
		pass
	elif argList.index("bam2plink") > argList.index("plink2tkrelated"):
		print("\n Error: When using both 'bam2plink' and 'plink2tkrelated' please always input 'plink2tkrelated' and its options after 'bam2plink'")
		quit()
if "bam2plink" in argList:
	argList.remove("bam2plink")
if "plink2tkrelated" in argList:
	argList.remove("plink2tkrelated")

short_options = "hr:L:e:p:m:b:tf:d:i:v"
long_options = ["help", "referenceGenome=", "gwvList=", "bamExtension=", "gwvPlink=", "minMQ=", "minBQ=", "excludeTerminalReadBases", "freqFile=", "dyads=", "ignoreThresh=", "verbose"]

try:
    options, args = getopt.getopt(argList, short_options, long_options)
except:
    print("\n Error: sys.argv error. Please make sure that you include the mandatory arguments and that they include the respective value, where necessary.")
    quit()

bam2plinkArgs = []
plink2tkrelatedArgs = []
for opt, value in options:
	if opt in ['-h', '--help']:
		print("""TKGWV2.py 
Version 1.0b - Released 07/2022

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
    - v, --verbose                    Use verbose mode and print extended information throughout file processing

Example run command:
  $ ./TKGWV2.py bam2plink --referenceGenome full_karyo.fa --gwvList 1000GP3_22M_noFixed_noChr.bed --bamExtension final.bam --gwvPlink DummyDataset_EUR_22M_noFixed plink2tkrelated --freqFile 1000GP3_EUR_22M_noFixed.frq""")
		
	elif opt in ['-r', '--referenceGenome']:
		bam2plinkArgs.append("=".join(['--referenceGenome',value]))
	elif opt in ['-L', '--gwvList']:
		bam2plinkArgs.append("=".join(['--gwvList',value]))
	elif opt in ['-e', '--bamExtension']:
		bam2plinkArgs.append("=".join(['--bamExtension',value]))
	elif opt in ['-p', '--gwvPlink']:
		bam2plinkArgs.append("=".join(['--gwvPlink',value]))
	elif opt in ['-t', '--excludeTerminalReadBases']:
		bam2plinkArgs.append('--excludeTerminalReadBases')
	elif opt in ['-m', '--minMQ']:
		bam2plinkArgs.append("=".join(["--minMQ",value]))
	elif opt in ['-b', '--minBQ']:
		bam2plinkArgs.append("=".join(["--minBQ",value]))
	elif opt in ['-f', '--freqFile']:
		plink2tkrelatedArgs.append("=".join(['--freqFile',value]))
	elif opt in ['-d', '--dyads']:
		plink2tkrelatedArgs.append("=".join(["--dyads",value]))
	elif opt in ['-i', '--ignoreThresh']:
		plink2tkrelatedArgs.append("=".join(["--ignoreThresh",value]))
	elif opt in ['-v', '--verbose']:
		plink2tkrelatedArgs.append("--verbose")

if len(bam2plinkArgs) > 0 and "bam2plink" not in sys.argv:
	print("\n Warning: Arguments for unused 'bam2plink' found. Will be ignored.")
if len(plink2tkrelatedArgs) > 0 and "plink2tkrelated" not in sys.argv:
	print("\n Warning: Arguments for unused 'plink2tkrelated' found. Will be ignored.")

if "bam2plink" in sys.argv:
	os.system("Rscript " + pywd + "/scripts/bam2plink.R " + " ".join(bam2plinkArgs) + " " + pywd)

if "plink2tkrelated" in sys.argv and "bam2plink" in sys.argv:
	os.system("Rscript " + pywd + "/scripts/plink2tkrelated.R " + " ".join(plink2tkrelatedArgs) + " bam2plink" + " " + pywd)
elif "plink2tkrelated" in sys.argv:
	os.system("Rscript " + pywd + "/scripts/plink2tkrelated.R " + " ".join(plink2tkrelatedArgs) + " " + pywd)
