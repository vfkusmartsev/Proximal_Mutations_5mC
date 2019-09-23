#!/usr/bin/python
# -*- coding: latin-1 -*-
from useful_funcs import *
import pandas as pd
import os
import sys
import pybedtools as pbt
import numpy as np
import itertools
pd.set_option('display.expand_frame_repr', False)
pd.options.display.max_rows = 10


		#########################################################################
		#									#
		#	GOAL: Replicate Qu et al., 2012, protocol. 			#
		#	Take CpGs from either sperm or H1 data, combine sites 		#
		#	within 10 bp from each other into CpG site blocks. 		#
		#	Then, overlap these blocks with SNPs from HapMap CEU, 		#
		#	1kG, or gnomAD, whilst removing any SNPs at CpG sites 		#
		#	themselves, and at TpG or CpA sites that would result 		#
		#	in a CpG. 							#
		#	Then compare the SNP incidence in methylated (>0.8) and  	#
		#	unmethylated (<0.2) site blocks. 				#
		#									#
		#########################################################################

input_path = '/rds/general/user/vfk16/home/SCRATCH/ax4/input/'
output_path = '/rds/general/user/vfk16/home/SCRATCH/ax4/output/'
####
#
#	Local dictionaries for different population genetics files
#	pop_gen[0] = snp_file
#	pop_gen[1] = snp_file columns
#	pop_gen[2] = snp_file_with_3nt_contexts
#	pop_gen[3] = dataset name
HapMap = [input_path + 'CEU_hapmap_2010_hg18_withAF.bed',
	['chrom','start','end','strand','alleles'],
	input_path + 'human_CEU_hapmap_hg18_withcontexts.txt',
	'HapMap']

kG = [input_path + '1kG.bed',
	['chrom','start','end','gene','ref','alt','AF','AC','AN','variant_type'],
	input_path + 'human_1kG_nomask_all_withcontexts.txt',
	'1kG']

gnomAD = [input_path + 'gnomad.txt',
	['chrom','start','end','gene','ref','alt','AF','AC','AN','variant_type'],
	input_path + 'gnomad_nomask_all_withcontexts.txt',
	'gnomAD']





####
#
#	Sets of in-silico experimental conditions to test all combinations of:
#		population genetics
#		methylation datasets
#		common vs rare variants (for 1kG and gnomAD)
#		frequency cut-offs (for 1kG and gnomAD)
#
combos = [z for z in itertools.product([kG, gnomAD],['H1','sperm'],['high','low'],[0.01,0.05])] + [j for j in itertools.product([HapMap],['H1','sperm'])]
###
#
#	Using PBS array jobs to crudely
#	parallelise protocol by condition
pop_gen = combos[int(sys.argv[1])-1][0]
meth_map = combos[int(sys.argv[1])-1][1]
if pop_gen != HapMap:
	frequency = combos[int(sys.argv[1])-1][2]
	threshold = combos[int(sys.argv[1])-1][3]


###
#
#	SNP filter function different for HapMap
#	because mutation labelled differently 
#	(C/T or A/G instead of 'transition')
if pop_gen == HapMap:
	def filter_SNPS_v2(df):
		####
		#
		#	Function takes dataframe
		#	Removes any SNP at CpG dinucleotide
		#	Removes any transitions SNP at the T in a TpG
		#	or A in a CpA, since would become a CpG,
		#	as discussed with Wei Qu
		if any(x in df.context for x in ('CG', 'Cg', 'cG', 'cg')):
			return False
		elif any(x in df.context[1:] for x in ('TG','Tg','tG','tg')) & (df.mutation == 'C/T'):
			return False
		elif any(x in df.context[:2] for x in ('CA','Ca','cA','ca')) & (df.mutation == 'A/G'):
			return False
		else:
			return True

else:
	def filter_SNPS_v2(df):
		####
		#
		#	Function takes dataframe
		#	Removes any SNP at CpG dinucleotide
		#	Removes any transitions SNP at the T in a TpG
		#	or A in a CpA, since would become a CpG,
		#	as discussed with Wei Qu
		if any(x in df.context for x in ('CG', 'Cg', 'cG', 'cg')):
			return False
		elif any(x in df.context[1:] for x in ('TG','Tg','tG','tg')) & (df.mutation_class == 'transition'):
			return False
		elif any(x in df.context[:2] for x in ('CA','Ca','cA','ca')) & (df.mutation_class == 'transition'):
			return False
		else:
			return True

def calculate_methylation_v2(df, limit = 5):
	####
	#
	#	Function takes dataframe of CpG site blocks
	#	For each blocks, it takes the coverage and
	#	methylation of each CpG in the blocks.
	#	Returns the average methylation of every CpG
	#	that has sufficient coverage.
	#	As discussed with Wei Qu
	values = [float(i) for i in df.coverage.split(',')]
	meth_vals = [float(i) for i in df.methylation.split(',')]
	meths = []
	for i in range(len(values)):
		if values[i] >= limit:
			meths.append(meth_vals[i])
	if len(meths) == 0:
		return -1.0
	else:
		return np.mean(meths)


def clean_h1_file(df):
	###
	#
	#	Takes the original H1 methylation file
	#	and returns the methylation and coverage.
	#	If the coverage is 0, then methylation
	#	is labelled -1.0, as discussed with Wei
	#	Qu.
	metadata = df.basecalls.split('_')
	nC = float(metadata[1])
	nT = float(metadata[2])
	coverage = nC + nT
	try:
		methylation = float(nC)/coverage
		return pd.Series([methylation, coverage])
	except ZeroDivisionError:
		return pd.Series([-1.0, coverage])




#####
#
#	Load up the original H1 and sperm basecalls
#	and coverage files. They have different formats
#	hence the if;else
#
if meth_map == 'sperm':
	#
	#	Hapmap uses hg18 reference genome,
	#	else hg19 version
	if pop_gen[3] == 'HapMap':
		sperm = pd.read_table(input_path + 'GSE30340_human_sperm_CpG_methylation_hg18.bedgraph',names = ['chrom','start','end','methylation'])
		coverage = pd.read_table(input_path + 'GSE30340_human_sperm_CpG_coverage_hg18.bedgraph', names = ['chrom','start','end','coverage'])
	else:
		sperm = pd.read_table(input_path + 'human_sperm_CpG_methylation_hg19.bed',names = ['chrom','start','end','methylation'])
		coverage = pd.read_table(input_path + 'human_sperm_CpG_coverage_hg19.bed', names = ['chrom','start','end','coverage'])
	sperm['coverage'] = coverage.coverage
	#
	#	Set up block edges
	sperm['start'] = sperm.start - 9
	sperm['end'] = sperm.end + 11
	#
	#	In case the CpG is right at the start of the genome
	#	(Can't have -ve positions)
	sperm['start'] = [0 if i < 0 else i for i in sperm.start]
	sperm.sort_values(['chrom','start'], inplace = True)
	sbed = pbt.BedTool.from_dataframe(sperm)
	#
	#	Overlap the blocks, keeping methylation and coverage information
	blocks = sbed.merge(c = (4,5), o = 'collapse', d = 1)
	blocks = pd.read_table(blocks.fn, names = ['chrom','start','end','methylation','coverage'])
	#
	#	Get the block sizes
	blocks['distance'] = blocks.end - blocks.start + 1
	#
	#	Calculate methylation as described by Wei Qu
	blocks['methylation2'] = blocks.apply(calculate_methylation_v2, axis = 1)
	valid = blocks[blocks.methylation2 >= 0.0].reset_index(drop = True)
	valid.drop(['methylation','coverage'], axis = 1, inplace = True)
	#
	#	Collect blocks with valid methylation and coverage,
	#	convert to bed for next step
	valid_bed = pbt.BedTool.from_dataframe(valid)


elif meth_map == 'H1':
	tmp = []
	for i in range(1,23):
		#
		#	Hapmap uses hg18 reference genome,
		#	else hg19 version
		if pop_gen[3] == 'HapMap':
			CG = pd.read_table(input_path + 'H1_c_basecalls/%s_H1_CG_basecalls' %('chr'+str(i)), names = ['chrom','start','end','strand','basecalls','removable'])
		else:
			CG = pd.read_table(input_path + 'H1_c_basecalls/%s_H1_CG_basecalls_hg19.bed' %('chr'+str(i)), names = ['chrom','start','end','strand','basecalls','removable'])
		CG[['methylation','coverage']] = CG.apply(clean_h1_file, axis = 1)
		CG.drop(['basecalls','removable','strand'], axis = 1, inplace = True)
		tmp.append(CG)
	H1 = pd.concat(tmp)
	H1.drop_duplicates(inplace = True)
	H1.sort_values(['chrom','start'], inplace = True)
	#
	#	Set up block edges
	H1['start'] = H1.start - 9
	H1['end'] = H1.end + 11
	#
	#	In case the CpG is right at the start of the genome
	#	(Can't have -ve positions)
	H1['start'] = [0 if i < 0 else i for i in H1.start]
	H1.sort_values(['chrom','start'],inplace = True)
	hbed = pbt.BedTool.from_dataframe(H1)
	#
	#	Overlap the blocks, keeping methylation and coverage information
	blocks = hbed.merge(c = (4,5), o = 'collapse', d = 1)
	blocks = pd.read_table(blocks.fn, names = ['chrom','start','end','methylation','coverage'])
	#
	#	Get the block sizes
	blocks['distance'] = blocks.end - blocks.start + 1
	#
	#	Calculate methylation as described by Wei Qu
	blocks['methylation2'] = blocks.apply(calculate_methylation_v2, axis = 1)
	valid = blocks[blocks.methylation2 >= 0.0].reset_index(drop = True)
	valid.drop(['methylation','coverage'], axis = 1, inplace = True)
	#
	#	Collect blocks with valid methylation and coverage,
	#	convert to bed for next step
	valid_bed = pbt.BedTool.from_dataframe(valid)




####
#
#	Load SNP files
snp_file = pd.read_table(pop_gen[0], names = pop_gen[1])
# for not HapMap also I haven't included the X and Y chromosomes NB
if pop_gen != HapMap:
	snp_file = snp_file[snp_file.allele_count != 0].reset_index(drop = True)
	if frequency == 'high':
		snp_file = snp_file[snp_file.freq >= threshold]
	elif frequency == 'low':
		snp_file = snp_file[snp_file.freq < threshold]


###
#
#	To each SNP add the 3nt context 
#	to be able to filter out CpG/TpG/CpA sites.
#	Very computationally slow, I personally ran this
#	in a separate terminal window and saved the
#	output before doing the rest.
#
if pop_gen != HapMap:
	snp_file['mutation_class'] = snp_file.apply(mutation_classifier, axis = 1)


snp_file['context'] = 'xxx'
for i in [str(z) for z in range(1,23)]+['X','Y']:
	print i
	ref = open_hg19_chrom(i)
	snp_file['context'] = snp_file.apply(get_context, chromosome = 'chr'+str(i), ref = ref, axis = 1)

###
#
#	Save this snp file for faster access
#	as context addition is slow
#snp_file.to_csv(output_path + '%s_nomask_all_withcontexts.txt' %species[3], sep = '\t', index = False)



###
#
#	Collect results from all the
#	mask and no mask conditions
cols = list(snp_file.columns)
results = []


######
#
#	0. REPLICATING EXACTLY WEI QU'S RESULT

##
#	SNP BED FILE
snp_file_bed = pbt.BedTool.from_dataframe(snp_file)
###
#
#	Overlap SNPs with correct CpG blocks
snps_blocks = snp_file_bed.window(valid_bed, w = 0)
snps_in_blocks = pd.read_table(snps_blocks.fn, names = cols + ['c2','bstart','bend','distance','methylation'])
##
#	Filter out SNPs as described in methods
sub = snps_in_blocks[snps_in_blocks.apply(filter_SNPS_v2, axis = 1)]

####
#
#	Collect data for visualisations
# MASK CONDITION
mask_string = 'no mask'

# SIZE OF METHYLATED CPG BLOCKS
MN = sum(valid[valid.methylation2 >= 0.8].distance)

# SIZE OF UNMETHYLATED CPG BLOCKS
UN = sum(valid[valid.methylation2 < 0.2].distance)

# NUMBER OF SNPS IN METHYLATED BLOCKS
MX = len(snps_in_blocks[snps_in_blocks.methylation >= 0.8])

# NUMBER OF SNPS IN UNMETHYLATED BLOCKS
UX = len(snps_in_blocks[snps_in_blocks.methylation < 0.2])

# NUMBR OF SNPS PASSING FILTER IN METHYLATED BLOCKS
MX_filter = len(sub[sub.methylation >= 0.8])

# NUMBER OF SNPS PASSING FILTER IN UNMETHYLATED BLOCKS
UX_filter = len(sub[sub.methylation < 0.2])

# SNPS PER BASE IN BLOCKS
M_rate = float(MX_filter)/MN
U_rate = float(UX_filter)/UN

# RATIO BETWEEN METHYLATED AND UNMETHYLATED SNP RATES
RelativeRisk = M_rate/U_rate


results.append([meth_map, pop_gen[3], mask_string, MN, UN, MX, UX, MX_filter, UX_filter, M_rate, U_rate, RelativeRisk])


###
#	MASK TIME
#
if pop_gen == HapMap:
	mask = pbt.BedTool(input_path + '20140520.strict_mask.autosomes_hg18.bed')
	repeatmask = pbt.BedTool(input_path + 'hg18_repeatmasker.bed')
else:
	mask = pbt.BedTool(input_path + '20140520.strict_mask.autosomes.bed')
	repeatmask = pbt.BedTool(input_path + 'hg19_repeatmasker.bed')



# use valid_bed as the original block file to mask

####
#
#	1) mask = Repeatmasker success

##
#	Apply mask
masked_bed = valid_bed.subtract(repeatmask)
masked_df = pd.read_table(masked_bed.fn, names = ['chrom','start','end','distance','methylation2'])
masked_df['end'] = masked_df.end - 1
masked_df['distance'] = masked_df.end - masked_df.start + 1
###
#
#	Overlap SNPs with masked CpG blocks
snps_blocks = snp_file_bed.window(masked_bed, w = 0)
snps_in_blocks = pd.read_table(snps_blocks.fn, names = cols + ['c2','bstart','bend','distance','methylation'])
##
#	Filter out SNPs as described in methods
sub = snps_in_blocks[snps_in_blocks.apply(filter_SNPS_v2, axis = 1)]


####
#
#	Collect data for visualisations
# MASK CONDITION
mask_string = 'repeatmasker pass'

# SIZE OF METHYLATED CPG BLOCKS
MN = sum(masked_df[masked_df.methylation2 >= 0.8].distance)

# SIZE OF UNMETHYLATED CPG BLOCKS
UN = sum(masked_df[masked_df.methylation2 < 0.2].distance)

# NUMBER OF SNPS IN METHYLATED BLOCKS
MX = len(snps_in_blocks[snps_in_blocks.methylation >= 0.8])

# NUMBER OF SNPS IN UNMETHYLATED BLOCKS
UX = len(snps_in_blocks[snps_in_blocks.methylation < 0.2])

# NUMBR OF SNPS PASSING FILTER IN METHYLATED BLOCKS
MX_filter = len(sub[sub.methylation >= 0.8])

# NUMBER OF SNPS PASSING FILTER IN UNMETHYLATED BLOCKS
UX_filter = len(sub[sub.methylation < 0.2])

# SNPS PER BASE IN BLOCKS
M_rate = float(MX_filter)/MN
U_rate = float(UX_filter)/UN

# RATIO BETWEEN METHYLATED AND UNMETHYLATED SNP RATES
RelativeRisk = M_rate/U_rate

results.append([meth_map, pop_gen[3], mask_string, MN, UN, MX, UX, MX_filter, UX_filter, M_rate, U_rate, RelativeRisk])



####
#
#	2) mask = Repeatmasker fail

##
#	Apply mask
masked_bed = valid_bed.intersect(repeatmask)
masked_df = pd.read_table(masked_bed.fn, names = ['chrom','start','end','distance','methylation2'])
masked_df['distance'] = masked_df.end - masked_df.start + 1
###
#
#	Overlap SNPs with masked CpG blocks
snps_blocks = snp_file_bed.window(masked_bed, w = 0)
snps_in_blocks = pd.read_table(snps_blocks.fn, names = cols + ['c2','bstart','bend','distance','methylation'])
##
#	Filter out SNPs as described in methods
sub = snps_in_blocks[snps_in_blocks.apply(filter_SNPS_v2, axis = 1)]

####
#
#	Collect data for visualisations
# MASK CONDITION
mask_string = 'repeatmasker fail'

# SIZE OF METHYLATED CPG BLOCKS
MN = sum(masked_df[masked_df.methylation2 >= 0.8].distance)

# SIZE OF UNMETHYLATED CPG BLOCKS
UN = sum(masked_df[masked_df.methylation2 < 0.2].distance)

# NUMBER OF SNPS IN METHYLATED BLOCKS
MX = len(snps_in_blocks[snps_in_blocks.methylation >= 0.8])

# NUMBER OF SNPS IN UNMETHYLATED BLOCKS
UX = len(snps_in_blocks[snps_in_blocks.methylation < 0.2])

# NUMBR OF SNPS PASSING FILTER IN METHYLATED BLOCKS
MX_filter = len(sub[sub.methylation >= 0.8])

# NUMBER OF SNPS PASSING FILTER IN UNMETHYLATED BLOCKS
UX_filter = len(sub[sub.methylation < 0.2])

# SNPS PER BASE IN BLOCKS
M_rate = float(MX_filter)/MN
U_rate = float(UX_filter)/UN

# RATIO BETWEEN METHYLATED AND UNMETHYLATED SNP RATES
RelativeRisk = M_rate/U_rate

results.append([meth_map, pop_gen[3], mask_string, MN, UN, MX, UX, MX_filter, UX_filter, M_rate, U_rate, RelativeRisk])



####
#
#	3) mask = Mappability success

##
#	Apply mask
masked_bed = valid_bed.intersect(mask)
masked_df = pd.read_table(masked_bed.fn, names = ['chrom','start','end','distance','methylation2'])
masked_df['distance'] = masked_df.end - masked_df.start + 1
###
#
#	Overlap SNPs with masked CpG blocks
snps_blocks = snp_file_bed.window(masked_bed, w = 0)
snps_in_blocks = pd.read_table(snps_blocks.fn, names = cols + ['c2','bstart','bend','distance','methylation'])
##
#	Filter out SNPs as described in methods
sub = snps_in_blocks[snps_in_blocks.apply(filter_SNPS_v2, axis = 1)]

####
#
#	Collect data for visualisations
# MASK CONDITION
mask_string = 'mappability pass'

# SIZE OF METHYLATED CPG BLOCKS
MN = sum(masked_df[masked_df.methylation2 >= 0.8].distance)

# SIZE OF UNMETHYLATED CPG BLOCKS
UN = sum(masked_df[masked_df.methylation2 < 0.2].distance)

# NUMBER OF SNPS IN METHYLATED BLOCKS
MX = len(snps_in_blocks[snps_in_blocks.methylation >= 0.8])

# NUMBER OF SNPS IN UNMETHYLATED BLOCKS
UX = len(snps_in_blocks[snps_in_blocks.methylation < 0.2])

# NUMBR OF SNPS PASSING FILTER IN METHYLATED BLOCKS
MX_filter = len(sub[sub.methylation >= 0.8])

# NUMBER OF SNPS PASSING FILTER IN UNMETHYLATED BLOCKS
UX_filter = len(sub[sub.methylation < 0.2])

# SNPS PER BASE IN BLOCKS
M_rate = float(MX_filter)/MN
U_rate = float(UX_filter)/UN

# RATIO BETWEEN METHYLATED AND UNMETHYLATED SNP RATES
RelativeRisk = M_rate/U_rate


results.append([meth_map, pop_gen[3], mask_string, MN, UN, MX, UX, MX_filter, UX_filter, M_rate, U_rate, RelativeRisk])



####
#
#	4) mask = Mappability success

##
#	Apply mask
masked_bed = valid_bed.subtract(mask)
masked_df = pd.read_table(masked_bed.fn, names = ['chrom','start','end','distance','methylation2'])
masked_df['end'] = masked_df.end - 1
masked_df['distance'] = masked_df.end - masked_df.start + 1
###
#
#	Overlap SNPs with masked CpG blocks
snps_blocks = snp_file_bed.window(masked_bed, w = 0)
snps_in_blocks = pd.read_table(snps_blocks.fn, names = cols + ['c2','bstart','bend','distance','methylation'])
##
#	Filter out SNPs as described in methods
sub = snps_in_blocks[snps_in_blocks.apply(filter_SNPS_v2, axis = 1)]

####
#
#	Collect data for visualisations
# MASK CONDITION
mask_string = 'mappability fail'

# SIZE OF METHYLATED CPG BLOCKS
MN = sum(masked_df[masked_df.methylation2 >= 0.8].distance)

# SIZE OF UNMETHYLATED CPG BLOCKS
UN = sum(masked_df[masked_df.methylation2 < 0.2].distance)

# NUMBER OF SNPS IN METHYLATED BLOCKS
MX = len(snps_in_blocks[snps_in_blocks.methylation >= 0.8])

# NUMBER OF SNPS IN UNMETHYLATED BLOCKS
UX = len(snps_in_blocks[snps_in_blocks.methylation < 0.2])

# NUMBR OF SNPS PASSING FILTER IN METHYLATED BLOCKS
MX_filter = len(sub[sub.methylation >= 0.8])

# NUMBER OF SNPS PASSING FILTER IN UNMETHYLATED BLOCKS
UX_filter = len(sub[sub.methylation < 0.2])

# SNPS PER BASE IN BLOCKS
M_rate = float(MX_filter)/MN
U_rate = float(UX_filter)/UN

# RATIO BETWEEN METHYLATED AND UNMETHYLATED SNP RATES
RelativeRisk = M_rate/U_rate

results.append([meth_map, pop_gen[3], mask_string, MN, UN, MX, UX, MX_filter, UX_filter, M_rate, U_rate, RelativeRisk])


####
#
#	Save output file
result_df = pd.DataFrame(results, columns = ['meth_map','pop_gen','mask','MN','UN','MX','UX','MX_noCpG','UX_noCpG','M_rate','U_rate','RelativeRisk'])
result_df.to_csv(output_path + '%s_%s_%sfreq_%sthresh_qureproduced_Aug2019.txt' %(pop_gen[1], meth_map, frequency, str(threshold)[3]), sep = '\t', header = True, index = False)







