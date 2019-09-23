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

######################
#
# A few useful functions that should probably be in the 
# global 'useful_funcs' import
#


def read_splitter(reads):
	tmp = reads.split(';')
	Cs = int(tmp[0][2:])
	Ts = int(tmp[1][2:])
	return Cs, Ts

def file_processor(mutant_df):
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

def add_chromatin(df, chromatin):
	columns = list(df.columns)
	dfbed = pbt.BedTool.from_dataframe(df)
	dfbed = dfbed.window(chromatin, w = 0)
	tmp = pd.read_table(dfbed.fn, names = columns + ['c2','s2','e2','state'])
	tmp.drop(['c2','s2','e2'], axis = 1, inplace = True)
	return tmp

def motif_finder(df, k, ref):
	if df.strand == '-':
		m = ref[df.start-k:df.start+k].translate(tab)[::-1]
	else:
		m = ref[df.start-k+1: df.start+k+1]

	if any(l in m for l in ['.','-','N']) or lowercase_check(m) or m[k-1] != 'C':
		return 'NA'
	if 'CG' in m[:k] or 'CG' in m[k:]:
		return 'NA'
	else:
		return m

def zone_distributor(df):
	zone = ''
	if (df.DDM > 0.7) & (df.DRD < 0.7):
		return 'DRD_target'
	elif (df.DRD > 0.7) & (df.DDM < 0.7):
		return 'DDM_target'
	elif (df.DDM > 0.7) & (df.DRD > 0.7):
		return 'both'
	else:
		return 'other'


def strand_selector(df, reference):
	base = reference[df.start].upper()
	if base == 'C':
		return '+'
	elif base == 'G':
		return '-'
	else:
		return 'NA'

arabidopsis_states_reduced = {
'E1':"H3.3",
'E2':"H3.3",
'E3':'H3K4me',
'E4':'H3K4me',
'E5':'H3K4me',
'E6':'H3K4me',
'E7':'H3K4me',
'E8':'H3K4me',
'E9':'H3K4me',
'E10':'H2A.Z',
'E11':'H3K27me3',
'E12':'H3K27me3',
'E13':'H3K27me3',
'E14':'H3K27me3',
'E15':'H3K27me3',
'E16':'Accessible',
'E17':'Accessible',
'E18':'Accessible',
'E19':'Accessible',
'E20':'Accessible',
'E21':'Accessible',
'E22':'H3K4me3',
'E23':'H3K4me3',
'E24':'H3K4me3',
'E25':'H3K4me3',
'E26':'H3K4me3',
'E27':'H3K4me3',
'E28':'H3K4me3',
'E29':'low signal',
'E30':'low signal',
'E31':'H3K9me2',
'E32':'H3K9me2',
'E33':'H3K9me2',
'E34':'H3K9me2',
'E35':'H3K9me2',
'E36':'H3K9me2'
}

def add_chromatin(df, chromatin):
	columns = list(df.columns)
	dfbed = pbt.BedTool.from_dataframe(df)
	dfbed = dfbed.window(chromatin, w = 0)
	tmp = pd.read_table(dfbed.fn, names = columns + ['c2','s2','e2','state'])
	tmp.drop(['c2','s2','e2'], axis = 1, inplace = True)
	return tmp

k = 5
mask_str = 'repeat_pass'
#chromosome = 'chr'+ str(int(sys.argv[1]))
mask = pbt.BedTool('/rds/general/user/vfk16/home/SCRATCH/ax4/input/repeatmasker_v2.tair10.bed')
#context = 'CG'
combos = [z for z in itertools.product(['chr'+str(i) for i in range(1,6)],['CG','CHG','CHH'])]
chromosome = combos[int(sys.argv[1])-1][0]
context = combos[int(sys.argv[1])-1][1]

path = '/rds/general/user/vfk16/home/SCRATCH/ax4/input/ara_meth_mutants/'


chromatin = pbt.BedTool('/rds/general/user/vfk16/home/SCRATCH/ax4/input/arabidopsis_chromatin.bed')

##
#	tair9 has same locations as tair10
#		where tair10 is only a genome annotation release, not a genome assembly
DDM = pd.read_table(path + 'tair9_GSM1014117_ddm1-2_bs-seq_%s_single-c.gff' %(context), names = ['chrom','lab','context','start','end','methylation','dot1','dot2','reads'])
DRD = pd.read_table(path + 'tair9_GSM1014120_drd1-7_bs-seq_%s_single-c.gff' %(context), names = ['chrom','lab','context','start','end','methylation','dot1','dot2','reads'])


DDM = file_processor(DDM)
DRD = file_processor(DRD)



mCG = pd.read_table('/rds/general/user/vfk16/home/SCRATCH/ax4/input/ara_mC/tair10_m%s_%s_10reads.txt' %(context, chromosome), names = ['chrom','start','end','strand','methylation'])
CG = pd.read_table('/rds/general/user/vfk16/home/SCRATCH/ax4/input/ara_mC/tair10_%s_%s_10reads.txt' %(context, chromosome), names = ['chrom','start','end','strand','methylation'])

anc = open_tair10_chrom(chromosome[3:])
mCG.sort_values('start', inplace = True)
CG.sort_values('start', inplace = True)

mCG = combo_masker(mCG, mask_str)
CG = combo_masker(CG, mask_str)

mCG = add_chromatin(mCG, chromatin)
CG = add_chromatin(CG, chromatin)

mCG['histone_mark'] = [arabidopsis_states_reduced[i] for i in mCG.state]
CG['histone_mark'] = [arabidopsis_states_reduced[i] for i in CG.state]


mCG['motif'] = mCG.apply(motif_finder, k = k, ref = anc, axis = 1)
mCG = mCG[mCG.motif != 'NA'].reset_index(drop = True)
mCG['CGstatus'] = 'methylated'
CG['motif'] = CG.apply(motif_finder, k = k, ref = anc, axis = 1)
CG = CG[CG.motif != 'NA'].reset_index(drop = True)
CG['CGstatus'] = 'unmethylated'

mCG.set_index(['chrom','start'], inplace = True)
CG.set_index(['chrom','start'], inplace = True)


mCG['DDM'] = DDM.methylation
mCG['DRD'] = DRD.methylation
mCG.dropna(inplace = True)
CG['DDM'] = DDM.methylation
CG['DRD'] = DRD.methylation

mCG['zone'] = mCG.apply(zone_distributor, axis = 1)
mCG.reset_index(inplace = True)
CG.reset_index(inplace = True)



if context != 'CG':
	mCG = mCG[~mCG.motif.str.contains('CG')].reset_index(drop = True)
	CG = CG[~CG.motif.str.contains('CG')].reset_index(drop = True)


cols = list(mCG.columns)
cols.append('paired')
states = list(set(mCG.zone))
histones = list(set(mCG.histone_mark))

paired = []
for state in states:
	for histone in histones:
		mtemp = mCG[(mCG.zone == state)&(mCG.histone_mark == histone)].reset_index(drop = True)
		ctemp = CG[CG.histone_mark == histone].reset_index(drop = True)
		ctemp['zone'] = state
		known = ctemp.append(mtemp, ignore_index = True)
		known.sort_values(['motif','start'], inplace = True)
		known['paired'] = False
		for motif in set(mtemp.motif):
			temp = np.array(known[known.motif == motif])
			if set(temp[:,8]) == {'methylated'}:
				continue
			else:
				maximum = len(temp)
				for i in range(maximum):
					if temp[i,8] == 'unmethylated':
						pos = temp[i,1]
						j = 1
						if i == 0:
							while i + j != maximum:
								if temp[i+j,8] == 'methylated' and temp[i+j,12] == False:
									paired.append(temp[i])
									paired.append(temp[i+j])
									temp[i,12] = True
									temp[i+j,12] = True
									break
								j+=1
						elif i + j == maximum:
							while i > j:
								if temp[i-j,8] == 'methylated' and temp[i-j,12] == False:
									paired.append(temp[i])
									paired.append(temp[i-j])
									temp[i,12] = True
									temp[i-j,12] = True
									break
								j+=1
						else:
							while i + j != maximum and i > j:
								if (temp[i+j,8] == 'methylated' and temp[i+j,12] == False) or (temp[i-j,8] == 'methylated' and temp[i-j,12] == False):
									break
								j+=1
							if i + j == maximum:
								if (temp[i-j,8] == 'methylated' and temp[i-j,12]==False):
									paired.append(temp[i])
									paired.append(temp[i-j])
									temp[i,12] = True
									temp[i-j,12] = True
							elif (temp[i+j,8] == 'methylated' and temp[i+j,12] == False) and (temp[i-j,8] == 'methylated' and temp[i-j,12]==False):
								behind = temp[i-j,1]
								ahead = temp[i+j,1]
								if (pos - behind) < (ahead - pos):
									paired.append(temp[i])
									paired.append(temp[i-j])
									temp[i,12] = True
									temp[i-j,12] = True
								elif (pos - behind) > (ahead - pos):
									paired.append(temp[i])
									paired.append(temp[i+j])
									temp[i,12] = True
									temp[i+j,12] = True
							elif temp[i+j,8] == 'methylated' and temp[i+j,12] == False:
								paired.append(temp[i])
								paired.append(temp[i+j])
								temp[i,12] = True
								temp[i+j,12] = True
							elif temp[i-j,8] == 'methylated' and temp[i-j,12] == False:
								paired.append(temp[i])
								paired.append(temp[i-j])
								temp[i,12] = True
								temp[i-j,12] = True




paired_all = pd.DataFrame(paired, columns = cols)
paired_all.drop(['paired'],axis = 1, inplace = True)
paired_all['sign'] = paired_all['strand']+'1'
paired_all['sign'] = pd.to_numeric(paired_all['sign'])


paired_all.to_csv('/rds/general/user/vfk16/home/SCRATCH/ax4/output/meth_mutants/ara_paired_DDMDRD_chromatin_%s_%s_k%s_%s.txt' %(mask_str, context, k,chromosome), sep = '\t', header = True)


























