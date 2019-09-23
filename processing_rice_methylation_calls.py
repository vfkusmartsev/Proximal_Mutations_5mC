#!/usr/bin/python
# -*- coding: latin-1 -*-
from useful_funcs import *
import pandas as pd
import pybedtools
import numpy as np
pd.set_option('display.expand_frame_repr', False)
pd.options.display.max_rows = 15


		#########################################################
		#							#
		#	GOAL: Preprocess methylation basecalls to	#
		#	have adequate read coverage (>=10 here) 	#
		#	and separate into unmethylated and methylated 	#
		#	sites 						#
		#	Returns files in format: 			#
		#	['CHROM, START, END, STRAND, METHYLATION']	#
		#							#
		#########################################################


path = '/Users/vkusmartsev/Documents/PhD_v2/input/rice/rice_methylome/'
output = '/Users/vkusmartsev/Documents/PhD_v2/input/rice/rice_methylome/'
read_threshold = 10

def read_upper_threshold(df, threshold = 0.7):

	
basecalls = pd.read_table(path+'GSM1039487_sample01_BSseq.txt', sep = ' ', names = ['chrom','start','strand','context','nC','nT'])
###
#
#	Calculate the read coverage
#	and the methylation stoichiometry
for context in ['CG','CHG','CHH']:
	for chrom in range(1,13):
		print chrom
		tmp = basecalls[basecalls.chrom == chrom]
		###
		#
		#	Subset basecalls to required context
		#	and with sufficient read coverage
		CGs_adequate = tmp[(tmp.context == context) & ((tmp.nC + tmp.nT) >= read_threshold)]

		####
		#
		#	Some sites are duplicated, with one row only having the C reads
		#	and one only having T reads, yet being at the same genomic
		#	position. Do two sets of pd.merge to combine
		#	
		CGs_adequate = pd.merge(CGs_adequate, CGs_adequate.groupby('start').nC.sum().reset_index(), on = 'start')
		CGs_adequate = pd.merge(CGs_adequate, CGs_adequate.groupby('start').nT.sum().reset_index(), on = 'start')
		CGs_adequate = CGs_adequate.drop_duplicates('start')
		CGs_adequate['methylation'] = CGs_adequate.nC_y/(CGs_adequate.nC_y + CGs_adequate.nT_y)

		###
		#
		#	Select methylated CpX
		methylated = CGs_adequate[CGs_adequate.methylation > 0.7].reset_index(drop = True)
		methylated['end'] = methylated['start'] + 1
		mCG = methylated[['chrom','start','end','strand','methylation']]
		mCG.to_csv(output+'m%s_rice_10reads_%s'%(context, chrom), sep = '\t', header = False, index = False)

		###
		#
		#	Select unmethylated CpX
		unmethylated = CGs_adequate[CGs_adequate.methylation < 0.2].reset_index(drop = True)
		unmethylated['end'] = unmethylated['start'] + 1
		CG = unmethylated[['chrom','start','end','strand','methylation']]
		CG.to_csv(output+'%s_rice_10reads_%s'%(context, chrom), sep = '\t', header = False, index = False)


# now need to liftover to IRGSP, then good to go!







