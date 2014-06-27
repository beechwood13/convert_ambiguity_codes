#!/bin/python
# -- coding: utf-8 -- 

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from random import randint

def eval_ambig_code(base):
	""" Take an ambiguity code and randomly choose a heterozygote base"""
	keto = "GT"
	amino = "AC"
	purine = "AG"
	pyrimidine = "CT"
	strong = "CG"
	weak = "AT"
	
	if base == "K" or base == "k":
		new_base =  keto[randint(0,1)]
	elif base == "M" or base == "m":
		new_base =  amino[randint(0,1)]
	elif base == "R" or base == "r":
		new_base =  purine[randint(0,1)]
	elif base == "Y" or base == "y":
		new_base =  pyrimidine[randint(0,1)]
	elif base ==  "S" or base == "s":
		new_base =  strong[randint(0,1)]
	elif base ==  "W" or base == "w":
		new_base =  weak[randint(0,1)]
	else:
		new_base = base
	
	return new_base
	
def alter_ambig_code(mutable_seq):
	"""Iterate through a sequence and alter the ambiguity codes
	returning a new sequence at the end"""
	for index, base in enumerate(mutable_seq):
 		if base not in "AaCcGgTtNn":
 			mutable_seq[index] = eval_ambig_code(base)
 		
	return mutable_seq

# set up output_handle
output_handle = open(sys.argv[2], "w")

# read in fasta, change seq and write out
for seq_record in SeqIO.parse(sys.argv[1], "fasta", IUPAC.ambiguous_dna):
	input_seq = seq_record.seq
	input_mutable = input_seq.tomutable()
	seq_record.seq = alter_ambig_code(input_mutable)
	SeqIO.write(seq_record, output_handle, "fasta")
	
output_handle.close()