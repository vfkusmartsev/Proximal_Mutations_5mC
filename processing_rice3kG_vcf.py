#!/usr/bin/python
# -*- coding: latin-1 -*-
from itertools import dropwhile
import pandas as pd
import sys

		#########################################################
		#							#
		#	GOAL: Preprocess population genetics dataset 	#
		#	into smaller, more manageable format		#
		#	Returns file in format: 			#
		#	[CHROM, START, END, REFERENCE, ALTERNATIVE, ALLELE_COUNT, ALLELE_NUMBER]
		#							#
		#########################################################

count = 0
##
#
#	File too large to keep data in RAM
#	so preopen output file to write data
#	to.
with open('/Users/vkusmartsev/Documents/PhD_v2/input/rice/rice_variants/3kRG.txt', 'a+') as final:
	filename = '/Users/vkusmartsev/Documents/PhD_v2/input/rice/rice_variants/NB_final_snp.vcf.vcf'
	with open(filename) as f:
		#########################################################################
		#									#
		#		Drops the commented lines that begin with a '#'		#
		#	Also reads line by line to save RAM (since too large to open) 	#
		for line in dropwhile(lambda line: line.startswith('#'), file(filename)):
			count += 1
			###
			#
			#	Takes tab-separated line, removes unnecessary
			#	whitespace right-hand side of line, and splits it
			#	by tabs
			line = line.rstrip()
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
					if i == '0/0':
						ref_homo += 2
					elif i == '1/1':
						alt_homo += 2
					elif i == '1/0' or i == '0/1':
						hetero += 1
				# add homozygotes and heterozygotes for alternative allele
				AC = alt_homo + hetero
				AN = alt_homo + ref_homo + hetero*2

				###
				#
				#	Write SNP to output file
				item = [chrom_no, pos, int(pos) + 1, ref, alt, AC, AN]
				final.write('\t'.join([str(o) for o in item]) +'\n')
			
			#########################################################
			#							#
			#	A small feedback to check progress: 		#
			#	counts up in 10000s up to total variant number	#
			if count % 10000 == 0:
				print count

















