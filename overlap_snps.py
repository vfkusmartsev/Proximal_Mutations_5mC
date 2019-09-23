#!/usr/bin/python
# -*- coding: latin-1 -*-
import pandas as pd
import pybedtools as pbt
import numpy as np
import matplotlib.pyplot as plt
import string
import seaborn as sns
import sys
pd.set_option('display.expand_frame_repr', False)

####################
#
# Local variable for finding reverse complement by sequence.translate(tab)[::-1]
#
old_chars = "ACGT"
replace_chars = "TGCA"
tab = string.maketrans(old_chars, replace_chars)

		#################################################################################
		#										#
		#	GOAL: Take the paired sites from scripts 'paired_sites_$SPECIES.py' 	#
		#	and overlap the corresponding singleton SNPs from population genetics 	#
		#	data, to estimate the mutation rate. Remove any SNPs that occur at 	#
		#	TpG or CpA dinucleotides, since may be cryptic CpGs hidden by a 	#
		#	polarization error. 							#
		#										#
		#################################################################################

input_path = '/rds/general/user/vfk16/home/SCRATCH/ax4/input/'
output_path = '/rds/general/user/vfk16/home/SCRATCH/ax4/output/'

###
#
#	Select the size k of the paired motifs
k = 5

#####
#
#	Species-specific information at the start of the script
#
#	Legend of column names:
#		'chrom' = Chromosome the SNP is on
#		'start' = Start position of the SNP
#		'end' 	= End position of the SNP
#		'gene' 	= ENCODE Genename the SNP is in if that information is available
#		'ref' 	= Reference/ancetral allele
#		'alt' 	= Alternative/derived allele
#		'AF' 	= Allele frequency of the SNP
#		'AC' 	= Allele count of the SNP in the population
#		'AN' 	= Allele number of sites in population
#		'variant_type' = Consequence of allele if information is available
#		'c2' 	= Chromosome the CpX is on
#		's2' 	= Start position of the CpX motif
#		'e2' 	= End position of the CpX motif
#		'strand' = Strand the CpX motif is on
#		'CGstatus' = If CpX is methylated or unmethylated
#		'methylation' = Stoichiometry of methylated to unmethylated reads of the focal CpX
#		'motif' = Motif the CpX is in
#		'focalC' = Position of focal CpX
#		'state' = Chromatin state the site is in
#		'sign' 	= Numerical representation of strand
#
#	Species information list is:
#	[SNP_FILE, COLUMNS_FOR_SNP_FILE, PAIRED_SITES_FILE, SPECIES_NAME, NUMBER_OF_CHROMOSOMES, ALLELE_COUNT_FOR_SINGLETONS, COLUMNS_FOR_OVERLAPS]

human = ['gnomad.txt',
	['chrom', 'start','end','gene','ref','alt','AF','AC','AN','variant_type'],
	'human_chromatin_paired_',
	'human',
	22,
	1,
	['chrom', 'start','end','gene','ref','alt','AF','AC','AN','variant_type','c2','s2','e2','strand','CGstatus','methylation','motif','focalC','state','sign']]

human10 = ['gnomad.txt',
	['chrom', 'start','end','gene','ref','alt','AF','AC','AN','variant_type'],
	'human_chromatin_paired_10_',
	'human_10',
	22,
	1,
	['chrom', 'start','end','gene','ref','alt','AF','AC','AN','variant_type','c2','s2','e2','strand','CGstatus','methylation','motif','focalC','state','sign']]

sperm = ['gnomad.txt',
	['chrom', 'start','end','gene','ref','alt','AF','AC','AN','variant_type'],
	'human_chromatin_sperm_paired_',
	'sperm',
	22,
	1,
	['chrom', 'start','end','gene','ref','alt','AF','AC','AN','variant_type','c2','s2','e2','CGstatus','methylation','motif','focalC','state']]

arabidopsis = ['arabdpsis_variants.txt',
	['chrom', 'start','end','ref','alt','AC','AN'],
	'ara_paired_10reads_',
	'arabidopsis',
	5,
	2,
	['chrom', 'start','end','ref','alt','AC','AN','c2','s2','e2','strand','CGstatus','methylation','motif','focalC','state','sign']]

rice10 = ['3kRG.txt',
	['chrom','start','end','ref','alt','AC','AN'],
	'/rds/general/user/vfk16/home/SCRATCH/ax4/output/v3/rice_paired_10reads_',
	'rice_10',
	12,
	1,
	['chrom','start','end','ref','alt','AC','AN','c2','s2','e2','strand','CGstatus','methlyation','motif','focalC','state','sign']]

aramechanism = ['arabdpsis_variants.txt',
	['chrom', 'start','end','ref','alt','AC','AN'],
	'/rds/general/user/vfk16/home/SCRATCH/ax4/output/RdDM_pairs/ara_paired_bymechanism_',
	'arabidopsis_bymechanism',
	5,
	2,
	['chrom', 'start','end','ref','alt','AC','AN','c2','s2','e2','strand','focalC','focalG','WT','dros1','dnrpd1','type2','zone','motif','CGstatus','sign']]

aramethylation = ['arabdpsis_variants.txt',
	['chrom', 'start','end','ref','alt','AC','AN'],
	'/rds/general/user/vfk16/home/SCRATCH/ax4/output/meth_mutants/ara_paired_byDDMDRD_',
	'arabidopsis_DDMDRD',
	5,
	2,
	['chrom', 'start','end','ref','alt','AC','AN','c2','s2','e2','strand','focalC','focalG','WT','motif','CGstatus','DDM','DRD','zone','sign']]

aramethhistone = ['arabdpsis_variants.txt',
	['chrom', 'start','end','ref','alt','AC','AN'],
	'/rds/general/user/vfk16/home/SCRATCH/ax4/output/meth_mutants/ara_paired_DDMDRD_chromatin_',
	'arabidopsis_DDMDRD_chromatin',
	5,
	2,
	['chrom', 'start','end','ref','alt','AC','AN','c2','s2','e2','strand','focalC','focalG','WT','state','histone_mark','motif','CGstatus','DDM','DRD','zone','sign']]

araRdDMhistone = ['arabdpsis_variants.txt',
	['chrom', 'start','end','ref','alt','AC','AN'],
	'/rds/general/user/vfk16/home/SCRATCH/ax4/output/meth_mutants/ara_paired_ROS1POL4_chromatin_',
	'arabidopsis_ROS1POL4_chromatin',
	5,
	2,
	['chrom', 'start','end','ref','alt','AC','AN','c2','s2','e2','strand','focalC','focalG','WT','state','histone_mark','motif','CGstatus','ROS1','NRPD','DKO','zone','sign']]

araros = ['arabdpsis_variants.txt',
	['chrom', 'start','end','ref','alt','AC','AN'],
	'/rds/general/user/vfk16/home/SCRATCH/ax4/output/RdDM_pairs/ara_paired_byROS_',
	'arabidopsis_byROS',
	5,
	2,
	['chrom', 'start','end','ref','alt','AC','AN','c2','s2','e2','strand','focalC','focalG','methylation','dROS','zone','motif','CGstatus','sign']]

ricemech = ['3kRG.txt',
	['chrom','start','end','ref','alt','AC','AN'],
	'/rds/general/user/vfk16/home/SCRATCH/ax4/output/RdDM_pairs/rice_paired_bymechanism_',
	'rice mechanism',
	12,
	1,
	['chrom','start','end','ref','alt','AC','AN','c2','s2','e2','strand','focalC','focalG','methylation','Cs','Ts','Ns','dROS','zone','motif','CGstatus','sign']]

ricemechv2 = ['3kRG.txt',
	['chrom','start','end','ref','alt','AC','AN'],
	'/rds/general/user/vfk16/home/SCRATCH/ax4/output/RdDM_pairs/rice_pairedv2_bymechanism_',
	'rice_mechanism_v2',
	12,
	1,
	['chrom','start','end','ref','alt','AC','AN','c2','s2','e2','strand','focalC','focalG','methylation','Cs','Ts','Ns','dROS','zone','CGstatus','motif','sign']]

def combine_reverse(df):
	forward = df[df.strand == '+'].reset_index()
	reverse = df[df.strand == '-'].reset_index()
	forward['start_reverse'] = forward.start+ 1
	temp = pd.merge(forward, reverse, left_on = 'start_reverse', right_on = 'start')
	temp['check'] = temp.apply(lambda x: True if x.motif_x == x.motif_y.translate(tab)[::-1] else False, axis = 1)
	complements = temp[temp.check == True].reset_index()
	reverse_strand_starts = list(complements.start_y)
	df = df[~df.start.isin(reverse_strand_starts)].reset_index(drop = True)
	return df

def TpG_CpA_polarised_checker(df):
	base = df.Position + 5
	if (df.motif[base:base+2] == 'TG') or (df.motif[base-1:base+1] == 'CA'):
		return False
	else:
		return True



####
#	NOTE for sperm no strand
#


####
#
#	Use sys.argv[1] to select species from command line
#	or when using PBS array jobs
species = [human10, arabidopsis, rice10][int(sys.argv[1])-1]



###
#
#	Load SNP file, as processed by processing scripts
snp_file = pd.read_table(input_path+species[0], names = species[1])
###
#	Correct SNP start site for 0 based counting
#	Make end == start for overlapping purposes
snp_file['start'] = snp_file['start'] - 1
snp_file['end'] = snp_file['end'] - 2
###
#	Remove alleles with 0 allele number
snp_file = snp_file[snp_file.AN != 0].reset_index(drop = True)
snp_file = snp_file.fillna('missing')
#
#	Make BedTool for overlappig
snp_bed = pbt.BedTool.from_dataframe(snp_file)

###
#
#	Can also do the sys.argv method above to select
#	different masks or contexts based on PBS array
#	if necessary (see pairing scripts)
mask = 'repeat_pass'
for context in ['CG','CHG']:
	####
	#
	#	Load the paired site information
	temp = []
	for i in range(1,species[4]+1):
		tmp = pd.read_table(input_path + species[2] + mask + '_'+context+'_k'+str(k)+'_chr'+str(i)+'.txt', index_col = 0) 
		tmp['chrom'] = 'chr'+str(i)
		temp.append(tmp)
	pairs = pd.concat(temp)

	#####
	#
	#	Get the start positions of the sequence motif
	#	for overlapping with SNP file
	pairs['rstart'] = pairs['start'] - (k-1)
	pairs['rend'] = pairs['start'] + k



	cols = list(pairs)

	####
	#
	#	Rearrange the columns into bedfile format
	#
	#	Again no strand in sperm, so done differently
	if species[3] = 'sperm':
		for i in ['rend','rstart','chrom']:
			cols.insert(0, cols.pop(cols.index(i)))
	else:
		for i in ['strand','rend','rstart','chrom']:
			cols.insert(0, cols.pop(cols.index(i)))

	pairs = pairs[cols]
	pairs.sort_values(['chrom','rstart'], inplace = True)
	print 'got pairs'
	pbed = pbt.BedTool.from_dataframe(pairs)

	####
	#
	#	overlap the snp file with the paired sites
	#	Use BedTool.window to keep SNP information
	ind_SNPS = snp_bed.window(pbed, w = 0)
	SNPS = pd.read_table(ind_SNPS.fn, names = species[6])

	###
	#
	#	Add AF if not already in dataframe
	#	Not necessary step
	if species[3] != 'human':
		SNPS['AF'] = SNPS.AC/SNPS.AN

	####
	#
	#	Add position of SNP w.r.t focal CpX
	if species[3] == 'sperm':
		SNPS['Position'] = SNPS.start - SNPS.focalC
	else:
		SNPS['Position'] = ((SNPS.start - SNPS.focalC)*SNPS.sign)


	####
	#
	#	Remove TpG and CpA sites
	SNPS['check'] = SNPS.apply(TpG_CpA_polarised_checker, axis = 1)
	if context != 'CG':
		# requires the check if at position 0 for CpHpX contexts
		SNPS = SNPS[(SNPS.check == True)|(SNPS.Position == 0)].reset_index(drop = True)
	else:	
		SNPS = SNPS[SNPS.check == True].reset_index(drop = True)
	

	####
	#
	#	Filter out the edge positions that do not
	#	have control of outer context. Separate 
	#	entries in dictionary incase wanted to filter 
	#	out specific sites for specific CpX contexts
	context_ignored_positions = {
		'CG':[-k,-(k-1),k],
		'CHG':[-k,-(k-1),k],
		'CHH':[-k,-(k-1),k]
	}
	SNPS = SNPS[~SNPS.Position.isin(context_ignored_positions[context])].reset_index(drop = True)

	###
	#
	#	Save file with all allele frequencies, as a check
	#	for allele frequency on effect size
	SNPS.to_csv(output_path + '%s_%s_%s_allSNPs_noCGTGCA_noreverse_allpos_k%s.txt' %(species[3], context, mask, k), sep = '\t', index = False)

	###
	#
	#	Filter only singletons, and resave as separate file
	SNPS = SNPS[SNPS.AC == species[5]].reset_index(drop = True)
	SNPS.to_csv(output_path + '%s_%s_%s_singletonSNPs_noCGTGCA_noreverse_allpos_k%s.txt' %(species[3], context, mask, k), sep = '\t', index = False)









