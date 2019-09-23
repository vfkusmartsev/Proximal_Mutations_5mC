#!/usr/bin/python
# -*- coding: latin-1 -*-
from useful_funcs import *
import pandas as pd
import os
import sys
import pybedtools as pbt
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
####################
#
# Whole script works on single chromosome -> allows crude parallelisation over the chromosomes.
#
combos = [z for z in itertools.product(['chr'+str(i) for i in range(1,23)],['none','mask_pass','mask_fail','repeat_pass','repeat_fail','both_pass','both_fail'],['CG','CHG','CHH'])]
chromosome = combos[int(sys.argv[1])][0]
mask_str = combos[int(sys.argv[1])][1]
context = combos[int(sys.argv[1])][2]


# mask contains regions of sufficient read quality in the reference mapping
#	intersect with mask = robust genomic regions
mask = pbt.BedTool(input_path + '20140520.strict_mask.autosomes.bed')
# repeatmask contains regions of repetitive elements
#	intersect with repeatmask = repetitive elements
repeatmask = pbt.BedTool(input_path + '/hg19_repeatmasker.bed')

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
	mask_function = {'mask_pass':masker(df, mask = mask),
		'mask_fail':mask_failed(df, mask = mask),
		'repeat_pass':mask_failed(df, mask = repeatmask),
		'repeat_fail':masker(df, mask = repeatmask),
		'both_pass':masker(mask_failed(df, mask = repeatmask),mask = mask),
		'both_fail':mask_failed(masker(df, mask = repeatmask), mask = mask)
	}
	return mask_function[mask_str]

####################
#
# Local variable for finding reverse complement by sequence.translate(tab)[::-1]
#
old_chars = "ACGT"
replace_chars = "TGCA"
tab = string.maketrans(old_chars, replace_chars)

# load the chromatin file as a BedTool for overlapping
chromatin = pbt.BedTool(input_path+'iHMM.M1K16.human_H1.bed')


###########
#
#	Load CpX information in format CHROM, START, END, METHYLATION
#		Where methylation is the ratio between methylated
#		and unmethylated reads in the BS-seq data.

# file already preprocessed to include coverage information (only sites >10 reads)
CGall = pd.read_table(input_path+'human_sperm_CpG_methylation_hg19_10reads.bed', names = ['chrom','start','end','methylation'])

#
# assign methylation state as 'methylated' or 'unmethylated'
CGall['status'] = CGall.apply(lambda x: 'methylated' if x.methylation > 0.7 else ('unmethylated' if x.methylation < 0.2 else 'drop'), axis = 1)
CGall = CGall[(CGall.status != 'drop') & (CGall.chrom == chromosome)].reset_index(drop = True)

# select chromosome and sort
CGall = CGall[CGall.chrom == chromosome].reset_index(drop = True)
CGall = CGall.sort_values(['start'])

CG = CGall[CGall.status == 'unmethylated'].reset_index(drop = True)
mCG = CGall[CGall.status == 'methylated'].reset_index(drop = True)


#######
#
#	Convert to BedTool to overlap with chromatin information
CG = pbt.BedTool.from_dataframe(CG)
mCG = pbt.BedTool.from_dataframe(mCG)

CGstates = CG.window(chromatin, w = 0)
mCstates = mCG.window(chromatin, w = 0)

CGdf = pd.read_table(CGstates.fn, names = ['chrom','start','end','methylation','status','chrom2','state_start','state_end','type','score','strand','thickStart','thickEnd','rgb'], index_col = False)
mCGdf = pd.read_table(mCstates.fn, names = ['chrom','start','end','methylation','status','chrom2','state_start','state_end','type','score','strand','thickStart','thickEnd','rgb'], index_col = False)

########
#
#	Remove sites mapping to unknown chromatin state
CGknown = CGdf[CGdf.type != '17_Unmap']
mCGknown = mCGdf[mCGdf.type != '17_Unmap']


########
#
#	Apply genomic masks if necessary
if mask_str != 'none':
	CGknown = combo_masker(CGknown, mask_str)
	mCGknown = combo_masker(mCGknown, mask_str)


########
#
#	Convert to numpy array for speed
#	of iterating
CG_array = np.array(CGknown)
mCG_array = np.array(mCGknown)


# open ancestral chromosome sequence to get 
# sequence for pairing
anc = open_anc_chrom(chromosome[3:])


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
	state = []
	methylation = []
	for i in CG_array:
		x = int(i[1])
		s = i[8]
		meth = i[3]
		m = anc[x-k+1:x+k+1]
		# remove sites with unknown/unconfident (lowercase) sequence
		if any(l in m for l in ['.','-','N']) or lowercase_check(m) or m[k-1:k+1] != 'CG':
			continue
		# remove sites with nearby adjacent CpG
		# that would skew mutational pattern
		if 'CG' in m[:k] or 'CG' in m[k:]:
			continue
		else:
			motif.append(m)
			start.append(x)
			state.append(s)
			methylation.append(meth)

	CG_df = pd.DataFrame({'motif':motif, 'start':start,'CGstatus':'unmethylated','type':state,'methylation':methylation},index = range(len(motif)))



	motif = []
	start = []
	state = []
	methylation = []
	for i in mCG_array:
		x = int(i[1])
		s = i[8]
		meth = i[3]
		m = anc[x-k+1:x+k+1]
		# remove sites with unknown/unconfident (lowercase) sequence
		if any(l in m for l in ['.','-','N']) or lowercase_check(m) or m[k-1:k+1] != 'CG':
			continue
		# remove sites with nearby adjacent CpG
		# that would skew mutational pattern
		if 'CG' in m[:k] or 'CG' in m[k:]:
			continue

		else:
			motif.append(m)
			start.append(x)
			state.append(s)
			methylation.append(meth)

	mCG_df = pd.DataFrame({'motif':motif, 'start':start,'CGstatus':'methylated','type':state, 'methylation':methylation},index = range(len(motif)))


	# compile the sites into single dataframe
	known = CG_df.append(mCG_df, ignore_index = True)
	known.sort_values(['motif','start'], inplace = True)
	known['paired'] = False

	# get list of chromatin states to iterate over
	states = list(set(known.type))

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
		pre_temp = known[known.type == state]
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
				#		if already paired, temp[,5] == True
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
								if temp[i+j,0] == 'methylated' and temp[i+j,5] == False:
									paired.append(temp[i])
									paired.append(temp[i+j])
									temp[i,5] = True
									temp[i+j,5] = True
									break
								j+=1
						elif i + j == maximum:
							while i > j:
								if temp[i-j,0] == 'methylated' and temp[i-j,5] == False:
									paired.append(temp[i])
									paired.append(temp[i-j])
									temp[i,5] = True
									temp[i-j,5] = True
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
								if (temp[i+j,0] == 'methylated' and temp[i+j,5] == False) or (temp[i-j,0] == 'methylated' and temp[i-j,5] == False):
									break
								j+=1
							if i + j == maximum:
								if (temp[i-j,0] == 'methylated' and temp[i-j,5]==False):
									paired.append(temp[i])
									paired.append(temp[i-j])
									temp[i,5] = True
									temp[i-j,5] = True
							#	If both i+j and i-j satisfy
							#		calculate distance from i to both
							#		select shortest distance
							elif (temp[i+j,0] == 'methylated' and temp[i+j,5] == False) and (temp[i-j,0] == 'methylated' and temp[i-j,5]==False):
								behind = temp[i-j,3]
								ahead = temp[i+j,3]
								if (pos - behind) < (ahead - pos):
									paired.append(temp[i])
									paired.append(temp[i-j])
									temp[i,5] = True
									temp[i-j,5] = True
								elif (pos - behind) > (ahead - pos):
									paired.append(temp[i])
									paired.append(temp[i+j])
									temp[i,5] = True
									temp[i+j,5] = True
							elif temp[i+j,0] == 'methylated' and temp[i+j,5] == False:
								paired.append(temp[i])
								paired.append(temp[i+j])
								temp[i,5] = True
								temp[i+j,5] = True
							elif temp[i-j,0] == 'methylated' and temp[i-j,5] == False:
								paired.append(temp[i])
								paired.append(temp[i-j])
								temp[i,5] = True
								temp[i-j,5] = True



	######
	#
	#	Convert paired sites into Dataframe
	paired_all = pd.DataFrame(paired, columns = ['CGstatus','methylation','motif','start','type','paired'])
	paired_all.drop(['paired'],axis = 1, inplace = True)

	# save output file
	paired_all.to_csv(output_path+'human_sperm_chromatin_paired_%s_%s_k%s_%s.txt' %(mask_str, context, k, chromosome), sep = '\t', header = True)

















