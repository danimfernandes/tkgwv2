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

short_options = "hr:L:e:p:m:b:tf:d:i:"
long_options = ["help", "referenceGenome=", "gwvList=", "bamExtension=", "gwvPlink=", "minMQ=", "minBQ=", "excludeTerminalReadBases", "freqFile=", "dyads=", "ignoreThresh="]

try:
    options, args = getopt.getopt(argList, short_options, long_options)
except:
    print("\n Error: sys.argv error. Please make sure that you include the mandatory arguments and that they include the respective value, where necessary.")
    quit()

bam2plinkArgs = []
plink2tkrelatedArgs = []
for opt, value in options:
	if opt in ['-h', '--help']:
		print("This is gonna be the help")
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

if len(bam2plinkArgs) > 0 and "bam2plink" not in sys.argv:
	print("\n Warning: Arguments for unused 'bam2plink' found. Will be ignored.")
if len(plink2tkrelatedArgs) > 0 and "plink2tkrelated" not in sys.argv:
	print("\n Warning: Arguments for unused 'plink2tkrelated' found. Will be ignored.")

if "bam2plink" in sys.argv:
	os.system("Rscript " + pywd + "/scripts/bam2plink_v0.3_pathToReady.R " + " ".join(bam2plinkArgs) + " " + pywd)

if "plink2tkrelated" in sys.argv and "bam2plink" in sys.argv:
	os.system("Rscript " + pywd + "/scripts/plink2tkrelated_v0.3_pathToReady.R " + " ".join(plink2tkrelatedArgs) + " bam2plink" + " " + pywd)
elif "plink2tkrelated" in sys.argv:
	os.system("Rscript " + pywd + "/scripts/plink2tkrelated_v0.3_pathToReady.R " + " ".join(plink2tkrelatedArgs) + " " + pywd)
