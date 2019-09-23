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

path = '/Users/vkusmartsev/Documents/PhD_v2/input/human/5mC/H1_c_basecalls/'
output = '/Users/vkusmartsev/Documents/PhD_v2/input/human/5mC/me_cytosines/'
read_threshold = 10


def read_upper_threshold(df, threshold = 0.7):
	# Filters basecalls that have methylation present at greater than threshold reads
	# Default threshold is 70%
	meth_fraction = float(df.nC)/(df.nC + df.nT)
	if meth_fraction > threshold:
		return True
	else:
		return False

def read_lower_threshold(df, threshold = 0.2):
	# Filters basecalls that have methylation present at greater than threshold reads
	# Default threshold is 70%
	meth_fraction = float(df.nC)/(df.nC + df.nT)
	if meth_fraction < threshold:
		return True
	else:
		return False

for context in ['CG']:
	for chrom in ['chr'+str(i) for i in range(1,23)]:
		print chrom
		basecalls = pd.read_table(path+'%s_h1_c_basecalls.txt' %chrom, sep = '\t', names = ['chrom','start','strand','context','nC','nT','nAll'])
		###
		#
		#	Subset basecalls to required context,
		#	and to required read coverage
		#	and calculated methylation stoichiometry
		CGs_adequate = basecalls[(basecalls.context == context) & ((basecalls.nC + basecalls.nT) >= read_threshold)]
		CGs_adequate['methylation'] = CGs_adequate.nC/(CGs_adequate.nC + CGs_adequate.nT)

		###
		#
		#	Select methylated CpX
		upper = CGs_adequate.apply(read_upper_threshold, axis = 1)
		methylated = CGs_adequate[upper]
		methylated['end'] = methylated['start'] + 1
		methylated['chrom'] = chrom
		mCG = methylated[['chrom','start','end','strand','methylation']]
		mCG.to_csv(output+'m%s_H1_%s.txt'%(context, chrom), sep = '\t', header = False, index = False)

		###
		#
		#	Select unmethylated CpX
		lower = CGs_adequate.apply(read_lower_threshold, axis = 1)
		unmethylated = CGs_adequate[lower]
		unmethylated['end'] = unmethylated['start'] + 1
		unmethylated['chrom'] = chrom
		CG = unmethylated[['chrom','start','end','strand','methylation']]
		CG.to_csv(output+'%s_H1_%s.txt'%(context, chrom), sep = '\t', header = False, index = False)


# now need to liftover to hg19, then good to go!








