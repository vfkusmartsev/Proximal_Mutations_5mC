#!/usr/bin/python
# -*- coding: latin-1 -*-
from useful_funcs import *
import pandas as pd
import pybedtools as pbt
import numpy as np
import sys
import itertools
pd.set_option('display.expand_frame_repr', False)
pd.options.display.max_rows = 15

old_chars = "ACGT"
replace_chars = "TGCA"
tab = string.maketrans(old_chars, replace_chars)


input_path = '/rds/general/user/vfk16/home/SCRATCH/ax4/input/'
output_path = '/rds/general/user/vfk16/home/SCRATCH/ax4/output/'

		#########################################################################
		#									#
		#	GOAL: Pair every location of a methylated CpX with an 		#
		#	unmethylated CpX that has the same ±knt sequence and same	#
		#	DDM1 activity. Filter out sites containing a nearby CpG 	#
		#	other than the focal CpX.					#
		#									#
		#########################################################################


######################
#
# A few useful functions that best defined locally
#


def read_splitter(reads):
	####
	#	Splits the reads columns into reads
	#	containing C and reads contain T
	#	at the focal C position
	#	(from BS-seq data)
	####
	tmp = reads.split(';')
	Cs = int(tmp[0][2:])
	Ts = int(tmp[1][2:])
	return Cs, Ts



def file_processor(mutant_df):
	######
	#	Removes unnecessary columns and formats
	#	mutant methylation dataframe to contain
	#	[CHROM, START] as the index to uniquely
	#	identify the methylation of each CpX in
	#	the mutant.
	#######
	mutant_df.drop(['lab','context','dot1','dot2'], axis = 1, inplace = True)
	mutant_df = mutant_df[~mutant_df.reads.isna()]
	mutant_df['C'], mutant_df['T'] = zip(*mutant_df['reads'].map(read_splitter))
	mutant_df['N'] = mutant_df['C'] + mutant_df['T']
	mutant_df['methylation'] = mutant_df.C/mutant_df.N
	mutant_df = mutant_df[mutant_df.N >= 10].reset_index(drop = True)
	mutant_df['start'] = mutant_df['start'].astype(int)
	mutant_df['start'] = mutant_df.start - 1
	mutant_df['end'] = mutant_df['end'].astype(int)
	mutant_df.drop('reads',axis = 1, inplace = True)
	mutant_df.set_index(['chrom','start'], inplace = True)
	return mutant_df



######
#
#	Set of 3 masking functions
#	that given dataframe, overlap
#	with mask using bedtools
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





def motif_finder(df, k, ref):
	########
	#	For each CpX returns the ±k nt motif
	#	up and downstream of the focal C
	#	Also removes unknown sequences, and
	#	sequences containing nearby CpGs
	#	that will skew the result
	#######
	# check strandedness for motif orientation
	if df.strand == '-':
		m = ref[df.start-k:df.start+k].translate(tab)[::-1]
	else:
		m = ref[df.start-k+1: df.start+k+1]

	# remove sites with unknown/unconfident (lowercase) sequence
	if any(l in m for l in ['.','-','N']) or lowercase_check(m) or m[k-1] != 'C':
		return 'NA'
	# remove sites with nearby adjacent CpG
	# that would skew mutational pattern
	if 'CG' in m[:k] or 'CG' in m[k:]:
		return 'NA'
	else:
		return m


def zone_distributor(df):
	######
	#	Function to check whether a sites requires
	#	DDM1 for methylation or is methylated solely
	#	by DRD and RdDM
	######
	zone = ''
	if (df.DDM > 0.7) & (df.DRD < 0.7):
		return 'DRD_target'
	elif (df.DRD > 0.7) & (df.DDM < 0.7):
		return 'DDM_target'
	elif (df.DDM > 0.7) & (df.DRD > 0.7):
		return 'neither'
	else:
		return 'both'


def strand_selector(df, reference):
	#####
	#	If strand data unknown, returns
	#	strand based on whether the base
	#	is a C or G (of the focal CpX)
	#####
	base = reference[df.start].upper()
	if base == 'C':
		return '+'
	elif base == 'G':
		return '-'
	else:
		return 'NA'



k = 5
mask = pbt.BedTool(input_path+'repeatmasker_v2.tair10.bed')
####################
#
# Whole script works on single chromosome -> allows crude parallelisation over the 5 chromosomes.
#
combos = [z for z in itertools.product(['chr'+str(i) for i in range(1,6)],['repeat_pass'],['CG','CHG','CHH'])]
chromosome = combos[int(sys.argv[1])][0]
mask_str = combos[int(sys.argv[1])][1]
context = combos[int(sys.argv[1])][2]



######
#
#	Load the mutant methylation data. The columns named [LAB, CONTEXT, DOT1, DOT2] are not needed and so placeholder names
#
#	NOTE tair9 and tair10 have the same genomic locations - tair10 only an annotation release, not a new genome assembly
DDM = pd.read_table(input_path + 'tair9_GSM1014117_ddm1-2_bs-seq_%s_single-c.gff' %(context), names = ['chrom','lab','context','start','end','methylation','dot1','dot2','reads'])
DRD = pd.read_table(input_path + 'tair9_GSM1014120_drd1-7_bs-seq_%s_single-c.gff' %(context), names = ['chrom','lab','context','start','end','methylation','dot1','dot2','reads'])

######
#
#	Process the files into format where
#	index is [CHROM, START] to unique identify
#	sites, and the rows contain methylation values
DDM = file_processor(DDM)
DRD = file_processor(DRD)


###########
#
#	Load CpX information in format CHROM, START, END, STRAND, METHYLATION
#		Where methylation is the ratio between methylated
#		and unmethylated reads in the BS-seq data.

# file already preprocessed to include coverage information (only sites >10 reads)
mCG = pd.read_table(input_path+'tair10_m%s_%s_10reads.txt' %(context, chromosome), names = ['chrom','start','end','strand','methylation'])
CG = pd.read_table(input_path+'tair10_%s_%s_10reads.txt' %(context, chromosome), names = ['chrom','start','end','strand','methylation'])
mCG.sort_values('start', inplace = True)
CG.sort_values('start', inplace = True)

####
#
#	Apply genomic masks if necessary
mCG = combo_masker(mCG, mask_str)
CG = combo_masker(CG, mask_str)

anc = open_tair10_chrom(chromosome[3:])

####
#	
#	Add the sequence motif ±k for each CpX
#	Drop sites with unknown sequence or
#	sequence containing other CpGs besides
#	the focal one
mCG['motif'] = mCG.apply(motif_finder, k = k, ref = anc, axis = 1)
mCG = mCG[mCG.motif != 'NA'].reset_index(drop = True)
mCG['CGstatus'] = 'methylated'
CG['motif'] = CG.apply(motif_finder, k = k, ref = anc, axis = 1)
CG = CG[CG.motif != 'NA'].reset_index(drop = True)
CG['CGstatus'] = 'unmethylated'

####
#
#	Set index to [CHROM, START]
#	for indexing the mutant methylation
#	values
mCG.set_index(['chrom','start'], inplace = True)
CG.set_index(['chrom','start'], inplace = True)

mCG['DDM'] = DDM.methylation
mCG['DRD'] = DRD.methylation
mCG.dropna(inplace = True)
CG['DDM'] = DDM.methylation
CG['DRD'] = DRD.methylation

######
#
#	Since unmethylated sites aren't methylated
#	can not apply DDM1 or DRD zones to them
#	because not applicable.
#	Instead, pairing mCG to closest unmethylated
#	CpX sites, and giving it the same zone for 
#	downstream normalisation purposes.
mCG['zone'] = mCG.apply(zone_distributor, axis = 1)
mCG.reset_index(inplace = True)
CG.reset_index(inplace = True)


#####
#
#	Filter out nearby CpGs in the paired sequences
#	to avoid skewing results
if context != 'CG':
	mCG = mCG[~mCG.motif.str.contains('CG')].reset_index(drop = True)
	CG = CG[~CG.motif.str.contains('CG')].reset_index(drop = True)


cols = list(mCG.columns)
cols.append('paired')
states = list(set(mCG.zone))
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
	mtemp = mCG[mCG.zone == state].reset_index(drop = True)
	# give unmethylated CpX position same zone for
	# downstream purposes
	CG['zone'] = state
	known = CG.append(mtemp, ignore_index = True)
	known.sort_values(['motif','start'], inplace = True)
	known['paired'] = False
	for motif in set(mtemp.motif):
		# subset to required sequence motif
		temp = np.array(known[known.motif == motif])

		# check if exist unmethylated site in this context
		if set(temp[:,6]) == {'methylated'}:
			continue
		else:
			maximum = len(temp)
			# in temp, find unmethylated sites,
			#	then radiate out above and below until nearest
			#	methylated site that has not yet been paired
			#		if already paired, temp[,10] == True
			for i in range(maximum):
				if temp[i,6] == 'unmethylated':
					pos = temp[i,1]
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
							if temp[i+j,6] == 'methylated' and temp[i+j,10] == False:
								paired.append(temp[i])
								paired.append(temp[i+j])
								temp[i,10] = True
								temp[i+j,10] = True
								break
							j+=1
					elif i + j == maximum:
						while i > j:
							if temp[i-j,6] == 'methylated' and temp[i-j,10] == False:
								paired.append(temp[i])
								paired.append(temp[i-j])
								temp[i,10] = True
								temp[i-j,10] = True
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
							if (temp[i+j,6] == 'methylated' and temp[i+j,10] == False) or (temp[i-j,6] == 'methylated' and temp[i-j,10] == False):
								break
							j+=1
						if i + j == maximum:
							if (temp[i-j,6] == 'methylated' and temp[i-j,10]==False):
								paired.append(temp[i])
								paired.append(temp[i-j])
								temp[i,10] = True
								temp[i-j,10] = True
						#	If both i+j and i-j satisfy
						#		calculate distance from i to both
						#		select shortest distance
						elif (temp[i+j,6] == 'methylated' and temp[i+j,10] == False) and (temp[i-j,6] == 'methylated' and temp[i-j,10]==False):
							behind = temp[i-j,1]
							ahead = temp[i+j,1]
							if (pos - behind) < (ahead - pos):
								paired.append(temp[i])
								paired.append(temp[i-j])
								temp[i,10] = True
								temp[i-j,10] = True
							elif (pos - behind) > (ahead - pos):
								paired.append(temp[i])
								paired.append(temp[i+j])
								temp[i,10] = True
								temp[i+j,10] = True
						elif temp[i+j,6] == 'methylated' and temp[i+j,10] == False:
							paired.append(temp[i])
							paired.append(temp[i+j])
							temp[i,10] = True
							temp[i+j,10] = True
						elif temp[i-j,6] == 'methylated' and temp[i-j,10] == False:
							paired.append(temp[i])
							paired.append(temp[i-j])
							temp[i,10] = True
							temp[i-j,10] = True



######
#
#	Convert paired sites into Dataframe
paired_all = pd.DataFrame(paired, columns = cols)
paired_all.drop(['paired'],axis = 1, inplace = True)
# make sign variable for downstream analysis
paired_all['sign'] = paired_all['strand']+'1'
paired_all['sign'] = pd.to_numeric(paired_all['sign'])

# save output file
paired_all.to_csv(output_path + 'ara_paired_byDDMDRD_%s_%s_k%s_%s.txt' %(mask_str, context, k,chromosome), sep = '\t', header = True)


























