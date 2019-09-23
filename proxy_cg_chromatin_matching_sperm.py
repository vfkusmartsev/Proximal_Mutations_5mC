#!/usr/bin/python
# -*- coding: latin-1 -*-
from useful_funcs import *
import pandas as pd
import os
import sys
import pybedtools as pbt
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



####################
#
# Whole script works on single chromosome -> allows crude parallelisation over the 22 chromosomes.
#
chromosome = 'chr'+ str(sys.argv[1])



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
	 ########################################
	#\					/#
	#\	Takes string, returns True 	/#
	#\	if contains lowercase letters	/#
	#\					/#
	if True in [w.islower() for w in string]:
		return True
	else:
		return False


def variant_checker(df, anc = 'anc', ignore_strand = True):
	 ################################################################
	#\								/#
	#\	Takes a dataframe with reference and alternative	/#
	#\ 	alleles and the allele count, compares it to the	/#
	#\	hg37 ancestral sequence, and returns the derived	/#
	#\	alternative allele and the correct count.		/#
	#\								/#
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
"""
print 'loading 1kG snp data....'
ogkg = pd.read_table('/rds/general/user/vfk16/home/SCRATCH/ax4/input/gnomad.txt', names = ['chrom', 'start','end','gene','ref','alt','freq','allele_count','allele_number','variant_type'])
kgsub = ogkg[(ogkg['chrom']==chromosome) & (ogkg.allele_count != 0)].reset_index()
kgsub['start'] = kgsub.start - 1
kgsub['freq'] = kgsub.freq.apply(float)
kgsub['allele_count'] = kgsub.allele_count.apply(int)
kgsub['allele_number'] = kgsub.allele_number.apply(int)
"""

chromatin = pbt.BedTool('/rds/general/user/vfk16/home/SCRATCH/ax4/input/iHMM.M1K16.human_H1.bed')

CGall = pd.read_table('/rds/general/user/vfk16/home/SCRATCH/ax4/input/human_sperm_CpG_methylation_hg19_10reads.bed', names = ['chrom','start','end','methylation'])

CGall['status'] = CGall.apply(lambda x: 'methylated' if x.methylation > 0.7 else ('unmethylated' if x.methylation < 0.2 else 'drop'), axis = 1)
CGall = CGall[(CGall.status != 'drop') & (CGall.chrom == chromosome)].reset_index(drop = True)


CGall = CGall[CGall.chrom == chromosome].reset_index(drop = True)
CGall = CGall.sort_values(['start'])

CG = CGall[CGall.status == 'unmethylated'].reset_index(drop = True)
mCG = CGall[CGall.status == 'methylated'].reset_index(drop = True)

CG = pbt.BedTool.from_dataframe(CG)
mCG = pbt.BedTool.from_dataframe(mCG)

CGstates = CG.window(chromatin, w = 0)
mCstates = mCG.window(chromatin, w = 0)

CGdf = pd.read_table(CGstates.fn, names = ['chrom','start','end','methylation','status','chrom2','state_start','state_end','type','score','strand','thickStart','thickEnd','rgb'], index_col = False)
mCGdf = pd.read_table(mCstates.fn, names = ['chrom','start','end','methylation','status','chrom2','state_start','state_end','type','score','strand','thickStart','thickEnd','rgb'], index_col = False)

CGknown = CGdf[CGdf.type != '17_Unmap']
mCGknown = mCGdf[mCGdf.type != '17_Unmap']


CG_array = np.array(CGknown)
mCG_array = np.array(mCGknown)


anc = open_anc_chrom(chromosome[3:])


for k in [5]:
	######################
	#
	# File for storing the numbers of things
	#




	#for every unmethylated CpG, find nearest methylated CpG
	# match sequence and context
	# look at nearby mutation rates


	motif = []
	start = []
	state = []
	methylation = []
	for i in CG_array:
		x = int(i[1])
		s = i[8]
		meth = i[3]
		m = anc[x-k+1:x+k+1]
		if any(l in m for l in ['.','-','N']) or lowercase_check(m) or m[k-1:k+1] != 'CG':
			continue
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
		if any(l in m for l in ['.','-','N']) or lowercase_check(m) or m[k-1:k+1] != 'CG':
			continue
		if 'CG' in m[:k] or 'CG' in m[k:]:
			continue

		else:
			motif.append(m)
			start.append(x)
			state.append(s)
			methylation.append(meth)

	mCG_df = pd.DataFrame({'motif':motif, 'start':start,'CGstatus':'methylated','type':state, 'methylation':methylation},index = range(len(motif)))



	known = CG_df.append(mCG_df, ignore_index = True)
	known.sort_values(['motif','start'], inplace = True)
	known['paired'] = False

	states = list(set(known.type))


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
						pos = temp[i,3]
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




	paired_all = pd.DataFrame(paired, columns = ['CGstatus','methylation','motif','start','type','paired'])
	paired_all.drop(['paired'],axis = 1, inplace = True)


	paired_all.to_csv('/rds/general/user/vfk16/home/SCRATCH/ax4/output/human_chromatin_sperm_paired_k%s_%s.txt' %(k,chromosome), sep = '\t', header = True)

"""
	#paired_all = pd.read_table('/scratch/vfk16/output/cg_chromatin_fullypaired_k%s_%s.txt' %(k, chromosome), sep = '\t', index_col = 0)

	Mutation_class = []
	Position = []
	Context = []
	X = []
	label = []
	N = ([len(paired_all[paired_all.CGstatus == 'methylated'])]*4 + [len(paired_all[paired_all.CGstatus == 'unmethylated'])]*4 ) * (k-2)*2 * 15

	for p in range(1,16):
		for l in range(1,k-1)+range(k+1,k*2-1):
			print 'checking position',l
			paired_all['mutpos'] = paired_all['start'] + (l-(k-1))
			paired_all['anc'] = [i[l] for i in paired_all['motif']]
			paired_all['cpg_check'] = ['nah' if (i[l:l+2] == 'TG' or i[l-1:l+1] == 'CA') else 'ok' for i in paired_all['motif']]
			Position += [str(l-(k-1))]*8
			Context += ['mCG']*4
			Context += ['CG']*4
			Mutation_class += (['transitions'] + ['GC_AT_transversions'] + ['AC_GT_transversions'] + ['all_transversion'])*2
			
			sub_snps = pd.merge(paired_all, kgsub, left_on = 'mutpos', right_on = 'start')
			


			if not sub_snps.empty:
				container = np.array([[i,j,f] for i,j,f in sub_snps.apply(variant_checker, axis = 1)])
				sub_snps['true_alt'] = container[:,0]
				sub_snps['true_AC'] = container[:,1]
				sub_snps['true_freq'] = container[:,2]
				sub_snps = sub_snps[(sub_snps.cpg_check == 'ok') & (sub_snps.true_AC != 'NA')]
				sub_snps['true_AC'] = sub_snps.true_AC.apply(int)
				sub_snps['true_freq'] = sub_snps.true_freq.apply(float)
				if p < 6:
					label += [float(p)/max(sub_snps.allele_number)]*8
					sub_snps = sub_snps[(sub_snps.true_AC == p)]
				else:
					p_lower = np.percentile(kgsub[kgsub.allele_count > 5].freq, (p-6)*10)
					p_upper = np.percentile(kgsub[kgsub.allele_count > 5].freq, ((p-6)*10)+10)
					label += [float(np.mean([p_lower,p_upper]))] * 8
					sub_snps = sub_snps[(sub_snps.true_freq > p_lower)&(sub_snps.true_freq < p_upper)]
					
				if not sub_snps.empty:
					sub_snps['mutation_type'] = sub_snps.apply(mutation_checker, axis = 1)

					m_snps = sub_snps[sub_snps.CGstatus == 'methylated']
					u_snps = sub_snps[sub_snps.CGstatus == 'unmethylated']
					X += [len(m_snps[m_snps.mutation_type == 'transition'])]
					mGC_AT = len(m_snps[m_snps.mutation_type == 'GC_AT_transversion'])
					mAC_GT = len(m_snps[m_snps.mutation_type == 'AC_GT_transversion'])
					X += [mGC_AT]
					X += [mAC_GT]
					X += [mGC_AT + mAC_GT]

					X += [len(u_snps[u_snps.mutation_type == 'transition'])]
					uGC_AT = len(u_snps[u_snps.mutation_type == 'GC_AT_transversion'])
					uAC_GT = len(u_snps[u_snps.mutation_type == 'AC_GT_transversion'])
					X += [uGC_AT]
					X += [uAC_GT]
					X += [uGC_AT + uAC_GT]
				else:
					X += [0]*8
			else:
				if p < 6:
					label += [float(p)/max(kgsub.allele_number)]*8
				else:
					p_lower = np.percentile(kgsub[kgsub.allele_count > 10].freq, (p-6)*10)
					p_upper = np.percentile(kgsub[kgsub.allele_count > 10].freq, ((p-6)*10)+10)
					label += [float(np.mean([p_lower,p_upper]))] * 8
				X += [0]*8

		

	mutation_data = pd.DataFrame({
		'Position':Position,
		'Mutation_class':Mutation_class,
		'Context':Context,
		'X':X,
		'N':N,
		'freq_bin':label},
		index = range(len(Position)))

	mutation_data.to_csv('/scratch/vfk16/output/gnomad_chromatin_muts_sperm_freqbinned_k%s_%s.txt' %(k, chromosome), sep = '\t', header = True)


"""


""" ##### This next bit is for looking at individual nucleotide motifs


		for nucleotide in ['G','A','T','C']:
			print 'checking context', nucleotide
			subset = paired_all[[i[l]==nucleotide for i in paired_all.motif]]
			if nucleotide == 'C':
				tmp = subset[[i[l+1]=='G' for i in subset.motif]]
				tmp['proxy_pos'] = tmp['start'] + (l-3)*subset['sign']
				nearby_CG.append(tmp)
				subset = subset[[i[l+1]!='G' for i in subset.motif]]
			if nucleotide == 'G' and l != 0:
				subset = subset[[i[l-1]!='C' for i in subset.motif]]
			subset['mutpos'] = subset['start'] + (l-3)*subset['sign']
			Position += [str(l-3)]*6
			Nt += [nucleotide]*6
			Context += ['mCG']*3
			Context += ['CG']*3
			Mutation_class += (['transitions'] + ['symmetric_transversions'] + ['transversions'])*2

			sub_snps = pd.merge(subset, kgsub, left_on = 'mutpos', right_on = 'start')
			

			N += [len(subset[subset.CGstatus == 'methylated'])]*3
			N += [len(subset[subset.CGstatus == 'unmethylated'])]*3 #need to filter by methylation
			if not sub_snps.empty:
				container = np.array([[i,j] for i,j in sub_snps.apply(variant_checker, anc=nucleotide, axis = 1)])
				sub_snps['true_alt'] = container[:,0]
				sub_snps['true_AC'] = container[:,1]

				sub_snps = sub_snps[(sub_snps['true_alt'] != 'NA') & (sub_snps.true_AC == '1')]
				m_snps = sub_snps[sub_snps.CGstatus == 'methylated']
				u_snps = sub_snps[sub_snps.CGstatus == 'unmethylated']
				X += [len(m_snps[m_snps.true_alt == transition_dict[nucleotide]])]
				X += [len(m_snps[m_snps.true_alt == sym_version_dict[nucleotide]])]
				X += [len(m_snps[m_snps.true_alt == transversion_dict[nucleotide]])]

				X += [len(u_snps[u_snps.true_alt == transition_dict[nucleotide]])]
				X += [len(u_snps[u_snps.true_alt == sym_version_dict[nucleotide]])]
				X += [len(u_snps[u_snps.true_alt == transversion_dict[nucleotide]])]
			else:
				X += [0]*6



	mutation_data = pd.DataFrame({
		'Position':Position,
		'Mutation_class':Mutation_class,
		'Nucleotide':Nt,
		'Context':Context,
		'X':X,
		'N':N},
		index = range(len(Position)))

	mutation_data.to_csv('/scratch/vfk16/output/proxy_cg_fullypaired_%s_%s.txt' %(k,chromosome), sep = '\t', header = True)

	
	nearby_CG = pd.concat(nearby_CG)




	mCmC = nearby_CG[(nearby_CG.CGstatus == 'methylated') & (nearby_CG.proxy_pos.isin(mCG.start))]
	mCmC['proxy_status'] = 'methylated'

	CmC = nearby_CG[(nearby_CG.CGstatus == 'unmethylated') & (nearby_CG.proxy_pos.isin(mCG.start))]
	CmC['proxy_status'] = 'methylated'

	CC = nearby_CG[(nearby_CG.CGstatus == 'unmethylated') & (nearby_CG.proxy_pos.isin(CG.start))]
	CC['proxy_status'] = 'unmethylated'

	mCC = nearby_CG[(nearby_CG.CGstatus == 'methylated') & (nearby_CG.proxy_pos.isin(CG.start))]
	mCC['proxy_status'] = 'unmethylated'

	known_nearby = pd.concat([mCmC, CmC, CC, mCC])


	Mutation_class = []
	Position = []
	Context = []
	proxy_status = []
	X = []
	N = []
	for l in range(2)+range(5,7):
		print 'checking position',l
		for nucleotide in ['CG']:
			subset = known_nearby[[i[l:l+2]==nucleotide for i in known_nearby.motif]]
			Position += [str(l-3)]*12
			Context += ['mCG']*6
			proxy_status += ['methylated']* 3 + ['unmethylated'] * 3
			Context += ['CG']*6
			proxy_status += ['methylated']* 3 + ['unmethylated'] * 3
			Mutation_class += (['transitions'] + ['symmetric_transversions'] + ['transversions'])*4

			sub_snps = pd.merge(subset, kgsub, left_on = 'proxy_pos', right_on = 'start')
			

			if not sub_snps.empty:
				container = np.array([[i,j] for i,j in sub_snps.apply(variant_checker, anc='C', axis = 1)])
				sub_snps['true_alt'] = container[:,0]
				sub_snps['true_AC'] = container[:,1]

				sub_snps = sub_snps[(sub_snps['true_alt'] != 'NA') & (sub_snps.true_AC == '1')]

				mm_snps = sub_snps[(sub_snps.CGstatus == 'methylated')&(sub_snps.proxy_status == 'methylated')]
				mu_snps = sub_snps[(sub_snps.CGstatus == 'methylated')&(sub_snps.proxy_status == 'unmethylated')]
				um_snps = sub_snps[(sub_snps.CGstatus == 'unmethylated')&(sub_snps.proxy_status == 'methylated')]
				uu_snps = sub_snps[(sub_snps.CGstatus == 'unmethylated')&(sub_snps.proxy_status == 'unmethylated')]

				N += [len(subset[(subset.CGstatus == 'methylated')&(subset.proxy_status == 'methylated')])]*3
				N += [len(subset[(subset.CGstatus == 'methylated')&(subset.proxy_status == 'unmethylated')])]*3
				N += [len(subset[(subset.CGstatus == 'unmethylated')&(subset.proxy_status == 'methylated')])]*3
				N += [len(subset[(subset.CGstatus == 'unmethylated')&(subset.proxy_status == 'unmethylated')])]*3

				X += [len(mm_snps[mm_snps.true_alt == transition_dict['C']])]
				X += [len(mm_snps[mm_snps.true_alt == sym_version_dict['C']])]
				X += [len(mm_snps[mm_snps.true_alt == transversion_dict['C']])]

				X += [len(mu_snps[mu_snps.true_alt == transition_dict['C']])]
				X += [len(mu_snps[mu_snps.true_alt == sym_version_dict['C']])]
				X += [len(mu_snps[mu_snps.true_alt == transversion_dict['C']])]

				X += [len(um_snps[um_snps.true_alt == transition_dict['C']])]
				X += [len(um_snps[um_snps.true_alt == sym_version_dict['C']])]
				X += [len(um_snps[um_snps.true_alt == transversion_dict['C']])]

				X += [len(uu_snps[uu_snps.true_alt == transition_dict['C']])]
				X += [len(uu_snps[uu_snps.true_alt == sym_version_dict['C']])]
				X += [len(uu_snps[uu_snps.true_alt == transversion_dict['C']])]
			else:
				X += [0]*12
				N += [len(subset[(subset.CGstatus == 'methylated')&(subset.proxy_status == 'methylated')])]*3
				N += [len(subset[(subset.CGstatus == 'methylated')&(subset.proxy_status == 'unmethylated')])]*3
				N += [len(subset[(subset.CGstatus == 'unmethylated')&(subset.proxy_status == 'methylated')])]*3
				N += [len(subset[(subset.CGstatus == 'unmethylated')&(subset.proxy_status == 'unmethylated')])]*3




	proxy_data = pd.DataFrame({
		'Position':Position,
		'Context':Context,
		'proxy_status':proxy_status,
		'Mutation_class':Mutation_class,
		'X':X,
		'N':N},
		index = range(len(Position)))

	proxy_data.to_csv('/scratch/vfk16/output/mC_proxy_mC_data_fullypaired_%s.txt' %chromosome, sep = '\t', header = True)
"""


