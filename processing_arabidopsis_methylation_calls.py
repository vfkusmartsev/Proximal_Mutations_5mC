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


path = '/Users/vkusmartsev/Documents/PhD_v2/input/arabidopsis/arabidopsis_methylation/'
output = '/Users/vkusmartsev/Documents/PhD_v2/input/arabidopsis/ara_mC/'
read_threshold = 10


def read_upper_threshold(df, threshold = 0.7):
	# Filters basecalls that have methylation present at greater than threshold reads
	# Default threshold is 70%
	meth_fraction = float(df.nC)/(df.nTot)
	if meth_fraction > threshold:
		return True
	else:
		return False

def read_lower_threshold(df, threshold = 0.2):
	# Filters basecalls that have methylation present at lower than threshold reads
	# Default threshold is 20%
	meth_fraction = float(df.nC)/(df.nTot)
	if meth_fraction < threshold:
		return True
	else:
		return False

basecalls = pd.read_table(path + 'At.Ws0.wm.glob.seed.bsseq.chr1-5.perC.txt', names = ['chr','strand','position','context','nC','nT','ratio'])
###
#
#	Calculate the read coverage
#	and the methylation stoichiometry
basecalls['nTot'] = basecalls['nC'] + basecalls['nT']
basecalls['chr'] = ['chr'+str(i) for i in basecalls['chr']]
basecalls['methylation'] = basecalls.nC/basecalls.nTot
for context in ['CG','CHG','CHH']:
	print context
	###
	#
	#	Subset basecalls to required context
	sub = basecalls[basecalls.context == context].reset_index(drop = True)
	for chrom in ['chr'+str(i) for i in range(1,6)]:
		print chrom
		temp = sub[sub.chr == chrom]
		###
		#
		#	Select sites with sufficient read coverage
		CGs_adequate = temp[temp.nTot >= read_threshold]

		###
		#
		#	Select methylated CpX
		upper = CGs_adequate.apply(read_upper_threshold, axis = 1)
		methylated = CGs_adequate[upper]
		methylated['start'] = methylated['position'] - 1
		methylated['chrom'] = chrom
		mCG = methylated[['chrom','start','position','strand','methylation']]
		print len(mCG)
		mCG.to_csv(output+'tair10_m%s_%s_%sreads.txt'%(context, chrom, read_threshold), sep = '\t', header = False, index = False)

		###
		#
		#	Select unmethylated CpX
		lower = CGs_adequate.apply(read_lower_threshold, axis = 1)
		unmethylated = CGs_adequate[lower]
		unmethylated['start'] = unmethylated['position'] - 1
		unmethylated['chrom'] = chrom
		CG = unmethylated[['chrom','start','position','strand','methylation']]
		print len(CG)
		CG.to_csv(output+'tair10_%s_%s_%sreads.txt'%(context, chrom, read_threshold), sep = '\t', header = False, index = False)










