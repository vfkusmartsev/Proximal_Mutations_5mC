#!/usr/bin/python
# -*- coding: latin-1 -*-
from useful_funcs import *
import pandas as pd
import pybedtools as pbt
import os
import sys
import itertools
import string
pd.set_option('display.expand_frame_repr', False)
pd.options.display.max_rows = 10



		#########################################################################
		#									#
		#	GOAL: Pair every location of a methylated CpX with an 		#
		#	unmethylated CpX that has the same ±knt sequence and same	#
		#	chromatin context. Filter out sites containing a nearby CpG 	#
		#	other than the focal CpX, and ones of unknown chromatin state	#
		#									#
		#########################################################################


input_path = '/rds/general/user/vfk16/home/SCRATCH/ax4/input/'
output_path = '/rds/general/user/vfk16/home/SCRATCH/ax4/output/'



replace_chroms = {'chr01':'chr1',
 'chr02':'chr2',
 'chr03':'chr3',
 'chr04':'chr4',
 'chr05':'chr5',
 'chr06':'chr6',
 'chr07':'chr7',
 'chr08':'chr8',
 'chr09':'chr9',
 'chr10':'chr10',
 'chr11':'chr11',
 'chr12':'chr12'}


####################
#
# Whole script works on single chromosome -> allows crude parallelisation over the chromosomes.
#
combos = [z for z in itertools.product(['chr'+str(i) for i in range(1,13)],['none','repeat_pass','repeat_fail'],['CG','CHG','CHH'])]
chromosome = combos[int(sys.argv[1])][0]
mask_str = combos[int(sys.argv[1])][1]
context = combos[int(sys.argv[1])][2]


# repeatmask contains regions of repetitive elements
#	intersect with repeatmask = repetitive elements
mask = pbt.BedTool(input_path + 'repeatmasker.IRGSP.bed')

###################
#
#	Local functions for masking dataframes
#
def masker(df, mask):
	columns = list(df.columns)
	dfbed = pbt.BedTool.from_dataframe(df)
	dfbed = dfbed.intersect(mask)
	return pd.read_table(dfbed.fn, names = columns)

def mask_failed(df, mask):
	columns = list(df.columns)
	dfbed = pbt.BedTool.from_dataframe(df)
	dfbed = dfbed.subtract(mask)
	return pd.read_table(dfbed.fn, names = columns)

def combo_masker(df, mask_str):
	mask_function = {'repeat_pass':mask_failed(df, mask = mask),
		'repeat_fail':masker(df, mask = mask)}
	return mask_function[mask_str]


####################
#
# Local variable for finding reverse complement by sequence.translate(tab)[::-1]
#
old_chars = "ACGT"
replace_chars = "TGCA"
tab = string.maketrans(old_chars, replace_chars)


# load the chromatin file as a BedTool for overlapping
chromatin = pbt.BedTool(input_path+'rice_chromatin.bed')

###########
#
#	Load CpX information in format CHROM, START, END, STRAND, METHYLATION
#		Where methylation is the ratio between methylated
#		and unmethylated reads in the BS-seq data.

# file already preprocessed to include coverage information (only sites >10 reads)
mCG = pd.read_table(input_path+'m%s_rice_10reads_%s_lifted.txt' %(context, chromosome[3:]), names = ['chrom','start','end','strand','methylation'])
CG = pd.read_table(input_path+'%s_rice_10reads_%s_lifted.txt' %(context, chromosome[3:]), names = ['chrom','start','end','strand','methylation'])

# chromosome names corrected for use with 
mCG['chrom'] = ['chr'+str(i) for i in mCG['chrom']]
CG['chrom'] = ['chr'+str(i) for i in CG['chrom']]
mCG.sort_values('start', inplace = True)
CG.sort_values('start', inplace = True)

# correct start site for zero-based counting
mCG['start'] = mCG['start'] - 1
CG['start'] = CG['start'] - 1


###
#
#	For CpGs, combine symmetrically methylated positions
#	to avoid excess duplicates
if context == 'CG':
	CG = combine_reverse(CG)
	mCG = combine_reverse(mCG)

#######
#
#	Convert to BedTool to overlap with chromatin information
mCG = pbt.BedTool.from_dataframe(mCG)
CG = pbt.BedTool.from_dataframe(CG)

CGstates = CG.window(chromatin, w = 0)
mCstates = mCG.window(chromatin, w = 0)

CGdf = pd.read_table(CGstates.fn, names = ['chrom','start','end','strand','methylation','chrom2','state_start','state_end','type','histone'], index_col = False)
mCGdf = pd.read_table(mCstates.fn, names = ['chrom','start','end','strand','methylation','chrom2','state_start','state_end','type','histone'], index_col = False)


########
#
#	Apply genomic masks if necessary
CGdf = combo_masker(CGdf, mask_str)
mCGdf = combo_masker(mCGdf, mask_str)

########
#
#	Convert to numpy array for speed
#	of iterating
CG_array = np.array(CGdf)
mCG_array = np.array(mCGdf)


# open chromosome sequence to get 
# sequence for pairing
anc = open_irgsp_chrom(chromosome[3:])

####
#
#	Can run for multiple k
#
#	±k for size of sequence to pair (including the focal CpX)
for k in [5]:
	# collect the motif and chromatin state of each CpX
	# for both unmethylated and methylated sites
	motif = []
	start = []
	strand = []
	state = []
	methylation = []
	for i in CG_array:
		x = int(i[1])
		d = i[3]
		meth = i[4]
		s = i[8]
		# check strandedness for motif orientation
		if d == '-':
			m = anc[x-k:x+k].translate(tab)[::-1]
		else:
			m = anc[x-k+1:x+k+1]
		# remove sites with unknown/unconfident (lowercase) sequence
		if any(l in m for l in ['.','-','N']) or lowercase_check(m) or m[k-1] != 'C':
			continue
		# remove sites with nearby adjacent CpG
		# that would skew mutational pattern
		if 'CG' in m[:k] or 'CG' in m[k:]:
			continue
		else:
			motif.append(m)
			start.append(x)
			strand.append(d)
			state.append(s)
			methylation.append(meth)

	CG_df = pd.DataFrame({'motif':motif, 'start':start,'strand':strand,'CGstatus':'unmethylated','methylation':methylation,'type':state},index = range(len(motif)))


	methylation = []
	motif = []
	start = []
	strand = []
	state = []
	for i in mCG_array:
		x = int(i[1])
		d = i[3]
		meth = i[4]
		s = i[8]
		# check strandedness for motif orientation
		if d == '-':
			m = anc[x-k:x+k].translate(tab)[::-1]
		else:
			m = anc[x-k+1:x+k+1]
		# remove sites with unknown/unconfident (lowercase) sequence
		if any(l in m for l in ['.','-','N']) or lowercase_check(m) or m[k-1] != 'C':
			continue
		# remove sites with nearby adjacent CpG
		# that would skew mutational pattern
		if 'CG' in m[:k] or 'CG' in m[k:]:
			continue

		else:
			motif.append(m)
			start.append(x)
			strand.append(d)
			state.append(s)
			methylation.append(meth)

	mCG_df = pd.DataFrame({'motif':motif, 'start':start,'strand':strand, 'CGstatus':'methylated','methylation':methylation,'type':state},index = range(len(motif)))


	# compile the sites into single dataframe
	known = CG_df.append(mCG_df, ignore_index = True)
	known.sort_values(['motif','start'], inplace = True)
	known['paired'] = False

	# get list of chromatin states to iterate over
	states = list(set(known['type']))


	#######
	#
	#	PAIRING METHOD
	#
	#	For each state and motif:
	#		Take unmethylated position
	#		Find nearest methylated sites
	#		Since data sorted, it is nearest in dataframe
	#		above or below
	paired = []
	for state in states:
		pre_temp = known[known['type'] == state]
		for motif in set(pre_temp.motif):
			# subset to required chromatin and sequence motif
			temp = np.array(pre_temp[pre_temp.motif == motif])

			# check if exist unmethylated site in this context
			if set(temp[:,0]) == {'methylated'}:
				continue
			else:
				maximum = len(temp)
				# if temp is empty, skip
				if maximum == 0:
					continue

				# in temp, find unmethylated sites,
				#	then radiate out above and below until nearest
				#	methylated site that has not yet been paired
				#		if already paired, temp[,6] == True
				for i in range(maximum):
					if temp[i,0] == 'unmethylated':
						pos = temp[i,3]
						j = 1
						####
						#
						#	Define edge cases when i = 0 (start of array of sites)
						#	or when the algorithm has hit the end of the array
						#	(i + j = maximum)
						#
						# 	When at the edge cases, can only search in one direction:
						#	Either above 0 (i + j), or below the maximum (i - j)
						#	j is the distance away from the unmethylated site being
						#	paired
						if i == 0:
							while i + j != maximum:
								if temp[i+j,0] == 'methylated' and temp[i+j,6] == False:
									paired.append(temp[i])
									paired.append(temp[i+j])
									temp[i,6] = True
									temp[i+j,6] = True
									break
								j+=1
						elif i + j == maximum:
							while i > j:
								if temp[i-j,0] == 'methylated' and temp[i-j,6] == False:
									paired.append(temp[i])
									paired.append(temp[i-j])
									temp[i,6] = True
									temp[i-j,6] = True
									break
								j+=1
						#####
						#
						# 	General algorithm:
						#		Iterate j (distance from site to be paired)
						#		Check if site i+j or i-j satisfies criteria
						#			i.e. if methylated and not yet paired
						else:
							while i + j != maximum and i > j:
								if (temp[i+j,0] == 'methylated' and temp[i+j,6] == False) or (temp[i-j,0] == 'methylated' and temp[i-j,6] == False):
									break
								j+=1
							if i + j == maximum:
								if (temp[i-j,0] == 'methylated' and temp[i-j,6]==False):
									paired.append(temp[i])
									paired.append(temp[i-j])
									temp[i,6] = True
									temp[i-j,6] = True
							#	If both i+j and i-j satisfy
							#		calculate distance from i to both
							#		select shortest distance
							elif (temp[i+j,0] == 'methylated' and temp[i+j,6] == False) and (temp[i-j,0] == 'methylated' and temp[i-j,6]==False):
								behind = temp[i-j,3]
								ahead = temp[i+j,3]
								if (pos - behind) < (ahead - pos):
									paired.append(temp[i])
									paired.append(temp[i-j])
									temp[i,6] = True
									temp[i-j,6] = True
								elif (pos - behind) > (ahead - pos):
									paired.append(temp[i])
									paired.append(temp[i+j])
									temp[i,6] = True
									temp[i+j,6] = True
							elif temp[i+j,0] == 'methylated' and temp[i+j,6] == False:
								paired.append(temp[i])
								paired.append(temp[i+j])
								temp[i,6] = True
								temp[i+j,6] = True
							elif temp[i-j,0] == 'methylated' and temp[i-j,6] == False:
								paired.append(temp[i])
								paired.append(temp[i-j])
								temp[i,6] = True
								temp[i-j,6] = True



	######
	#
	#	Convert paired sites into Dataframe
	paired_all = pd.DataFrame(paired, columns = ['CGstatus','methylation','motif','start','strand','type','paired'])
	paired_all.drop(['paired'],axis = 1, inplace = True)
	# make sign variable for downstream analysis
	paired_all['sign'] = paired_all['strand']+'1'
	paired_all['sign'] = pd.to_numeric(paired_all['sign'])


	# save output file
	paired_all.to_csv(output_path + 'rice_paired_10reads_%s_%s_k%s_%s.txt' %(mask_str, context, k, chromosome), sep = '\t', header = True)


















