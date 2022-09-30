#!/usr/bin/python3
__author__ = 'Daniel Fernandes'

import os
import sys
import random
import re

cwd = os.getcwd()

exclTerminals = sys.argv[1]

for i in os.listdir(cwd):
	if '.pileupsamtools.gwv.txt' == i[-23:]:
		sample_id, rest = re.split('.pileupsamtools.gwv.txt', i)
		pedout = open(sample_id + '.ped', 'w')
		pileupin = open(i, 'r')		
		snpcount, discar = re.split(' ', os.popen('wc -l ' + i).read())
		snpcount = int(snpcount)
		pedout.write(sample_id + ' ' + sample_id + ' 0 0 0 1 ')
		pileupin = open(i, 'r')
		counter = 1
		countbads = 0
		for line in pileupin:
			if counter <= snpcount:
				line = line.rstrip()
				if line.split("\t")[3] == "0":
						pedout.write('0 0 ')
						countbads += 1
				else:
					chro, coord, refb, cov, rbase, qual = line.split('\t')
					counter += 1
					rbase = rbase.replace(".",refb.capitalize())
					rbase = rbase.replace(",",refb.capitalize())
					if int(cov) == 1:
						if r"*" in rbase or r"+" in rbase or r"-" in rbase:
							pedout.write('0 0 ')
						elif r"$" in rbase:
							rbase = rbase[0]
							if exclTerminals == 'True': 
								pedout.write('0 0 ')
							else:
								pedout.write(rbase.capitalize() + ' ' + rbase.capitalize() + ' ')

						elif r"^" in rbase:
							rbase = rbase[2]
							if exclTerminals == 'True':
								pedout.write('0 0 ')
							else:
								pedout.write(rbase.capitalize() + ' ' + rbase.capitalize() + ' ')
						else:
							pedout.write(rbase.capitalize() + ' ' + rbase.capitalize() + ' ')
					
					## If there is more than one call for this position
					elif int(cov) > 1:
						listBases = list(rbase)
						## Remove indels
						while r"-" in listBases:
							for i in listBases:
								if i == r"-":
									indexi = listBases.index(r"-")
									delL = listBases[indexi+1]
									del listBases[indexi:(indexi+int(delL)+2)] # delete indel annotation
									del listBases[indexi-1] # delete indel
						while r"+" in listBases:
							for i in listBases:
								if i == r"+":
									indexi = listBases.index(r"+")
									delL = listBases[indexi+1]
									del listBases[indexi:(indexi+int(delL)+2)] # delete indel annotation
									del listBases[indexi-1] # delete indel

						## Remove start and end markers (^ and $) and read calls if requested 
						while r"^" in listBases:
							for i in listBases:
								if i == r"^":
									indexi = listBases.index(r"^")
									del listBases[indexi] # delete ^
									del listBases[indexi] # delete the base quality after the ^
									if exclTerminals == 'True' and len(listBases) > 0:
										del listBases[indexi-1] # delete the actual variant 
						while r"$" in listBases:
							for i in listBases:
								if i == r"$":
									indexi = listBases.index(r"$")
									del listBases[indexi] # delete $
									if exclTerminals == 'True' and len(listBases) > 0:
										del listBases[indexi-1] # delete the actual variant 

						if len(listBases) != 0:
							chosenbase = random.choice(listBases)
							pedout.write(chosenbase.capitalize() + ' ' + chosenbase.capitalize() + ' ')
						else:
							countbads += 1
							pedout.write('0 0 ')
			else:
				pass
		pedout.close()

