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
		#	[CHROM, START, END, GENE_NAME, REFERENCE, ALTERNATIVE, ALLELE_FREQUENCY, ALLELE_COUNT, ALLELE_NUMBER, VARIANT_TYPE]
		#							#
		#########################################################


path = '/Users/vkusmartsev/Documents/PhD_v2/input/human/vcf/'
count = 0
##
#
#	File too large to keep data in RAM
#	so preopen output file to write data
#	to.
with open(path + '1000_genomes/1kG.bed', 'a+') as final:
	filename = path + '1000_genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf'
	with open(filename) as f:
		#########################################################################
		#									#
		#		Drops the commented lines that begin with a '#'		#
		#	Also reads line by line to save RAM (since too large to open) 	#
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
				###
				#
				#	Splitted line item 7 contains all metadata
				#	Pull out AC, AF, AN, and CSQ (consequence)
				#	information
				info_temp = line[7].split(';')
				AC_index = [z for z in info_temp if z.startswith('AC=')]
				AC = AC_index[0].replace('AC=','').split(',')[0]
				AF_index = [z for z in info_temp if z.startswith('AF=')]
				AF = AF_index[0].replace('AF=','').split(',')[0]
				AN_index = [z for z in info_temp if z.startswith('AN=')]
				AN = AN_index[0].replace('AN=','').split(',')[0]

				###
				#
				#	If consequence information doesn't exist,
				#	return NA in said column
				CSQ = [i.split('|') for i in info_temp if i.startswith('CSQ=')]
				if CSQ:
					CSQ = CSQ[0]
					gene_name = CSQ[4]
					variant = CSQ[1]
				else:
					gene_name = 'NA'
					variant = 'NA'
				###
				#
				#	Write SNP to output file
				item = [chrom_no, pos, str(int(pos) + 1), gene_name, ref, alt, AF, AC, AN, variant]
				final.write('\t'.join([str(o) for o in item]) +'\n')
			
			#########################################################
			#							#
			#	A small feedback to check progress: 		#
			#	counts up in 1000s up to total variant number	#
			#	Fills terminal with '#' signs 			#
			if count % 10000 == 0:
				sys.stdout.write('#')
				sys.stdout.flush()





