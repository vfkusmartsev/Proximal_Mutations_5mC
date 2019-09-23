#!/usr/bin/python
# -*- coding: latin-1 -*-
from itertools import dropwhile
import pandas as pd

		#########################################################
		#							#
		#	GOAL: Preprocess population genetics dataset 	#
		#	into smaller, more manageable format		#
		#	Returns file in format: 			#
		#	[CHROM, START, END, REFERENCE, ALTERNATIVE, ALLELE_COUNT, ALLELE_NUMBER]
		#							#
		#########################################################

path = '/Users/vkusmartsev/Documents/PhD_v2/input/arabidopsis/'
filename = path + '1001genomes_snp-short-indel_only_ACGTN.vcf'
data = []
with open(filename) as f:
	#########################################################################
	#									#
	#		Drops the commented lines that begin with a '#'		#
	for line in dropwhile(lambda line: line.startswith('#'), file(filename)):

		###
		#
		#	Takes tab-separated line and splits it
		line = line.split("\t")

		###
		#
		#	Collect info from splitted line
		chrom_no = 'chr' + line[0]
		pos = line[1]
		ref = line[3]
		alt = line[4]

		#################################################################################################
		#												#
		#	Only select single nucleotide variants, ignoring indels and multiple variants 		#
		if len(ref)==1 and len(alt)==1:
			#####
			#	
			#	Metadata contains genotypes of all samples.
			#	Reference allele is 0, alternative allele is 1
			#	Summate number of homozygous, and heterozygous samples
			#	
			metadata = line[9:]
			ref_homo = 0
			alt_homo = 0
			hetero = 0
			for i in metadata:
				test = i.split(':')
				if test[0] == '0|0':
					ref_homo += 2
				elif test[0] == '1|1':
					alt_homo += 2
				elif test[0] == '1|0' or test[0] == '0|1':
					hetero += 1

			###	There are a lot of positions with ./. for missing alleles

			# add homozygotes and heterozygotes for alternative allele
			AC = alt_homo + hetero
			AN = alt_homo + ref_homo + 2*hetero
			data.append([chrom_no, pos, str(int(pos) + 1), ref, alt, AC, AN])


		
		#########################################################
		#							#
		#	A small feedback to check progress: 		#
		#	counts up in 1000s up to total variant number	#
		if len(data) % 1000 == 0:
			print len(data)

###
#
#	Save output
final = pd.DataFrame(data, columns = ['chr','start','end','ref','alt','AC','AN'])
final.to_csv(path + 'arabdpsis_variants.txt', sep = '\t', header = False)










