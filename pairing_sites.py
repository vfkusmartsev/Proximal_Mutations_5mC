#!/usr/bin/python
# -*- coding: latin-1 -*-
from useful_funcs import *
import pandas as pd
import os
import sys
import pybedtools as pbt
import itertools
import string
pd.set_option('display.expand_frame_repr', False)
pd.options.display.max_rows = 10


# Run using chromatin_master.sh

		#########################################################################
		#									#
		#	GOAL: Pair every location of a methylated CpG with an 		#
		#	unmethylated CpG that has the same ±knt sequence. If 		#
		#	there is a CpG nearby, filter these and consider them		#
		#	separately. If the latter case, make sure the nearby 		#
		#	CpG context is known and also matched.	Then, for each 		#
		#	set of pairs, calculate the mutation rate by looking at 	#
		#	the occurrence of derived singletons at each position 		#
		#	in the motif. 							#
		#									#
		#########################################################################


#########
#
#	matching exactly to k-1, can analyse mutation rates at k-2.
#	^ for reference about what k means in the file naming convention

input_path = '/rds/general/user/vfk16/home/SCRATCH/ax4/input/'
output_path = '/rds/general/user/vfk16/home/SCRATCH/ax4/output/'
####################
#
# Whole script works on single chromosome -> allows crude parallelisation over chromosomes.
#

##########
#
#	Species sub-data
#
# species = [NAME, CHROMATIN_FILE, ]

humanH1 = ['human_H1',
	input_path + 'iHMM.M1K16.human_H1.bed',
	]
sperm = ['human_sperm',
	input_path + 'iHMM.M1K16.human_H1.bed',
	]
arabidopsis = ['arabidopsis',
	]
rice = ['rice',
	]





##
#
# For selecting mask and chromosome condition when parallelising using qsub on computing cluster

# creates every combo of chromosome, mask, and context
combos = [z for z in itertools.product(['chr'+str(i) for i in range(1,23)],['none','mask_pass','mask_fail','repeat_pass','repeat_fail','both_pass','both_fail'],['CG','CHG','CHH'],[humanH1,sperm,arabidopsis,rice])]
chromosome = combos[int(sys.argv[1])][0]
mask_str = combos[int(sys.argv[1])][1]
context = combos[int(sys.argv[1])][2]
species = combos[int(sys.argv[1])][3]
# mask contains regions of sufficient read quality in the reference mapping
#	intersect with mask = robust genomic regions
mask = pbt.BedTool(input_path + '20140520.strict_mask.autosomes.bed')
# repeatmask contains regions of repetitive elements
#	intersect with repeatmask = repetitive elements
repeatmask = pbt.BedTool(input_path + '/hg19_repeatmasker.bed')

###################
#
#	Functions for masking dataframes
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

######################
#
# A few useful functions that should probably be in the 
# global 'useful_funcs' import
#
def lowercase_check(string):
	#########################################
	#					#
	#	Takes string, returns True 	#
	#	if contains lowercase letters	#
	#					#
	if True in [w.islower() for w in string]:
		return True
	else:
		return False


def variant_checker(df, anc = 'anc', ignore_strand = False):
	#################################################################
	#								#
	#	Takes a dataframe with reference and alternative	#
	# 	alleles and the allele count, compares it to the	#
	#	hg37 ancestral sequence, and returns the derived	#
	#	alternative allele and the correct count.		#
	#								#
	if anc == 'anc':
		anc = df.anc
	if not ignore_strand:
		if df.strand == '+':
			if df.ref == anc:
				return df.alt, df.allele_count, df.freq
			elif df.alt == anc:
				return df.ref, df.allele_number - df.allele_count, 1 - df.freq
			else:
				return 'NA', 'NA', 0.
		else:
			if df.ref.translate(tab) == anc:
				return df.alt.translate(tab), df.allele_count, df.freq
			elif df.alt.translate(tab) == anc:
				return df.ref.translate(tab), df.allele_number - df.allele_count, 1 - df.freq
			else:
				return 'NA', 'NA', 0.
	else:
		if df.ref == anc:
			return df.alt, df.allele_count, df.freq
		elif df.alt == anc:
			return df.ref, df.allele_number - df.allele_count, 1 - df.freq
		else:
			return 'NA', 'NA', 0.

def mutation_checker(df):
	transition_dict = {'C':'T','G':'A','T':'C','A':'G'}
	GC_AT_dict = {'C':'G','G':'C','T':'A','A':'T'}
	AC_GT_dict = {'C':'A','G':'T','A':'C','T':'G'}

	if df.true_alt == transition_dict[df.anc]:
		return 'transition'
	elif df.true_alt == GC_AT_dict[df.anc]:
		return 'GC_AT_transversion'
	elif df.true_alt == AC_GT_dict[df.anc]:
		return 'AC_GT_transversion'

###################
#
#	Pre open all the data files to avoid loop
#

# write separate loading functions

chromatin = pbt.BedTool(species[1])

CG = pbt.BedTool('/rds/general/user/vfk16/home/SCRATCH/ax4/input/me_cytosines/%s_H1_%s_lifted.txt' %(context, chromosome))
mCG = pbt.BedTool('/rds/general/user/vfk16/home/SCRATCH/ax4/input/me_cytosines/m%s_H1_%s_lifted.txt' %(context, chromosome))

CGstates = CG.window(chromatin, w = 0)
mCstates = mCG.window(chromatin, w = 0)

CGdf = pd.read_table(CGstates.fn, names = ['chrom','start','end','strand','chrom2','state_start','state_end','type','score','thickStart','thickEnd','rgb'], index_col = False)
mCGdf = pd.read_table(mCstates.fn, names = ['chrom','start','end','strand','chrom2','state_start','state_end','type','score','thickStart','thickEnd','rgb'], index_col = False)

CGknown = CGdf[(CGdf.type != '17_Unmap') & (CGdf.chrom == chromosome)]
mCGknown = mCGdf[(mCGdf.type != '17_Unmap') & (mCGdf.chrom == chromosome)]


CGknown = combo_masker(CGknown, mask_str)
mCGknown = combo_masker(mCGknown, mask_str)

CG_array = np.array(CGknown)
mCG_array = np.array(mCGknown)



anc = open_anc_chrom(chromosome[3:])

# k is ±length of motif to pair by (including the focal C)
for k in [5]:
	# collect the motif and chromatin state of each CpX
	# for both unmethylated and methylated sites
	motif = []
	start = []
	strand = []
	state = []
	for i in CG_array:
		x = int(i[1])
		d = i[3]
		s = i[7]
		if d == '-':
			m = anc[x-k:x+k].translate(tab)[::-1]
		else:
			m = anc[x-k+1:x+k+1]
		if any(l in m for l in ['.','-','N']) or m[k-1] != 'C':
			continue
		if 'CG' in m[:k].upper() or 'CG' in m[k:].upper():
			continue
		else:
			motif.append(m)
			start.append(x)
			strand.append(d)
			state.append(s)
	CG_df = pd.DataFrame({'motif':motif, 'start':start,'strand':strand,'CGstatus':'unmethylated','type':state},index = range(len(motif)))


	motif = []
	start = []
	strand = []
	state = []
	for i in mCG_array:
		x = int(i[1])
		d = i[3]
		s = i[7]
		if d == '-':
			m = anc[x-k:x+k].translate(tab)[::-1]
		else:
			m = anc[x-k+1:x+k+1]
		if any(l in m for l in ['.','-','N']) or lowercase_check(m) or m[k-1] != 'C':
			continue
		if 'CG' in m[:k] or 'CG' in m[k-1:]:
			continue

		else:
			motif.append(m)
			start.append(x)
			strand.append(d)
			state.append(s)
	mCG_df = pd.DataFrame({'motif':motif, 'start':start,'strand':strand, 'CGstatus':'methylated','type':state},index = range(len(motif)))


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
			temp = np.array(pre_temp[pre_temp.motif == motif])
			if set(temp[:,0]) == {'methylated'}:
				continue
			else:
				maximum = len(temp)
				if maximum == 0:
					continue
				for i in range(maximum):
					if temp[i,0] == 'unmethylated':
						pos = temp[i,2]
						j = 1
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
							elif (temp[i+j,0] == 'methylated' and temp[i+j,5] == False) and (temp[i-j,0] == 'methylated' and temp[i-j,5]==False):
								behind = temp[i-j,2]
								ahead = temp[i+j,2]
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




	paired_all = pd.DataFrame(paired, columns = ['CGstatus','motif','start','strand','type','paired'])
	paired_all.drop(['paired'],axis = 1, inplace = True)
	# make sign variable for downstream analysis
	paired_all['sign'] = paired_all['strand']+'1'
	paired_all['sign'] = pd.to_numeric(paired_all['sign'])


	# save output file
	paired_all.to_csv('/rds/general/user/vfk16/home/SCRATCH/ax4/output/human_chromatin_paired_%s_%s_k%s_%s.txt' %(mask_str, context, k,chromosome), sep = '\t', header = True)

