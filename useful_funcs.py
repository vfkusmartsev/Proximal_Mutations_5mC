#!/usr/bin/python
# -*- coding: latin-1 -*-

import numpy as np 
import re
from itertools import dropwhile
from collections import Counter
import string
#################################################################################
#	Set up code for complementing positive strand nucleotide sequence 	#
#	To get complement, type "sequence".translate(tab)[::-1]			#
#################################################################################
old_chars = "ACGT"
replace_chars = "TGCA"
tab = string.maketrans(old_chars, replace_chars)





def open_mm10_chrom(chromosome):
	 #################################################################
	#\								 /#
	#\	Opens the reference sequence for 'chromosome'		 /#
	#\	Returns a single string of the sequence 		 /#
	#\								 /#
	 #################################################################
	path = "/rds/general/user/vfk16/home/SCRATCH/ax4/input/mm10_chromosomes/"
	filename = path + "chr" + str(chromosome) + ".fa"
	with open(filename) as f:
		lines = f.readlines()[1:]
		reference = ''.join(lines)
		reference = reference.replace("\n","")
		return reference

def open_irgsp_chrom(chromosome):
	 #################################################################
	#\								 /#
	#\	Opens the reference sequence for 'chromosome'		 /#
	#\	Returns a single string of the sequence 		 /#
	#\								 /#
	 #################################################################
	path = "/rds/general/user/vfk16/home/SCRATCH/ax4/input/IRGSP_v1/"
	filename = path + "IRGSP_v1" + str(chromosome-1).zfill(2) + ".fa"
	with open(filename) as f:
		lines = f.readlines()[1:]
		reference = ''.join(lines)
		reference = reference.replace("\n","")
		return reference

		

def open_tair10_chrom(chromosome):
	 #################################################################
	#\								 /#
	#\	Opens the reference sequence for 'chromosome'		 /#
	#\	Returns a single string of the sequence 		 /#
	#\								 /#
	 #################################################################
	path = "/rds/general/user/vfk16/home/SCRATCH/ax4/input/TAIR10/"
	filename = path + "chr" + str(chromosome) + ".fas"
	with open(filename) as f:
		lines = f.readlines()[1:]
		reference = ''.join(lines)
		reference = reference.replace("\n","")
		return reference





def open_anc_chrom(chromosome):
	 #################################################################
	#\								 /#
	#\	Opens the ancestral sequence for 'chromosome'		 /#
	#\	Returns a single string of the sequence 		 /#
	#\								 /#
	 #################################################################
	path = "/rds/general/user/vfk16/home/SCRATCH/ax4/human_ancestor_GRCh37_e59/"
	filename = path + "human_ancestor_" + str(chromosome) + ".fa"
	with open(filename) as f:
		lines = f.readlines()[1:]
		ancestor = ''.join(lines)
		ancestor = ancestor.replace("\n","")
		return ancestor



def open_hg19_chrom(chromosome):
	 #################################################################
	#\								 /#
	#\	Opens the reference sequence for 'chromosome'		 /#
	#\	Returns a single string of the sequence 		 /#
	#\								 /#
	 #################################################################
	path = "/rds/general/user/vfk16/home/SCRATCH/ax4/GRCh37_hg19/"
	filename = path + "chr" + str(chromosome) + ".fa"
	with open(filename) as f:
		lines = f.readlines()[1:]
		reference = ''.join(lines)
		reference = reference.replace("\n","")
		return reference


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


def variant_checker(df, anc = 'anc', ignore_strand = False):
	#################################################################
	#								#
	#	Takes a dataframe with reference and alternative	#
	# 	alleles and the allele count, compares it to the	#
	#	hg37 ancestral sequence, and returns the derived	#
	#	alternative allele and the correct count.		#
	#								#
	#	Returns derived allele; allele count; allele freq 	#
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



def combine_reverse(df):
	#########################################################
	#							#
	#	Function to combine CpGs from opposite strands 	#
	forward = df[df.strand == '+'].reset_index()
	reverse = df[df.strand == '-'].reset_index()
	forward['start_reverse'] = forward.start+ 1
	temp = pd.merge(forward, reverse, left_on = 'start_reverse', right_on = 'start')
	reverse_strand_starts = list(temp.start_y)
	df = df[~df.start.isin(reverse_strand_starts)].reset_index(drop = True)
	return df


def translate(sequence):
	 ################################################################
	#\								/#
	#\	Takes coding strand sequence and translates it 		/#
	#\	into a protein sequence.				/#
	#\								/#
	 ################################################################
	my_seq = sequence.lower()
	
	#setting up the codon table (can import text file)
	CodonTable = {}
	table_file = open("/rds/general/user/vfk16/home/SCRATCH/ax4/trans_tbl")
	for line in table_file:
		(key,value) = line.split()
		CodonTable[key] = value
	table_file.close()

	#selecting the codons and translating
	n = len(my_seq)# - start_pos
	n_amino = n//3
	overhang = n%3
	if overhang != 0:
		print "Sequence not a multiple of 3! There will be an excess of %d nucleotides" %(overhang)
	prot_seq = ""
	for i in range(n_amino):
		codon = my_seq[3*i: (3*i)+3]

		if len(codon) != 3:
			print "not enough nucleotides for next amino acid!"
			break
		if "n" in codon or "." in codon or "-" in codon:
			amino_acid = 'X'
		else:	
			amino_acid = CodonTable.get(codon)
		# if amino_acid == "*":
		# 	print "Reached a stop codon!"
		# 	break
		prot_seq += amino_acid
	return prot_seq







