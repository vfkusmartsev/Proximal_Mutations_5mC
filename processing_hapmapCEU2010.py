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
		#	[CHROM, START, END, STRAND, ALLELES]		#
		#							#
		#########################################################


path = '/Users/vkusmartsev/Documents/PhD_v2/input/human/vcf/'


###
#
#	Remove SNPs where the only allele known is 'N'
#		i.e. unknown
def remove_only_Ns(df, col):
	###
	#	df.alleles in format 'X>Y'
	#	so bases are positions 0 and 2
	allele1, allele2 = df.alleles[0], df.alleles[2]
	totalstring = ''
	###
	#
	#	Iterates through all samples (NA_some_number)
	#	and concatenates the alleles into long string.
	#	If string only contains N, drop.
	for i in [z for z in col if z.startswith('NA')]:
		totalstring += df[i]
	if set(totalstring) == {'N'}:
		return False
	else:
		return True


collect = []
for chrom in ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']:
	print chrom
	filename = path + 'HapMap_2010/genotypes_%s_CEU_r28_nr.b36_fwd.txt' %chrom
	tmp = pd.read_table(filename, sep = ' ')
	columns = list(tmp.columns)
	###
	#
	#	Select SNPs that have known alleles
	tmp = tmp[tmp.apply(remove_only_Ns, col = columns, axis = 1)]
	tmp.drop([i for i in columns if i not in ['chrom','pos','strand','alleles']], axis = 1, inplace = True)
	##
	#	Add end position for bed format
	tmp['end'] = tmp.pos + 1
	tmp = tmp[['chrom','pos','end','strand','alleles']]
	collect.append(tmp)


final = pd.concat(collect)
final.to_csv(path + 'CEU_hapmap_2010_hg18_withAF.bed', sep = '\t',header = False, index = False)




