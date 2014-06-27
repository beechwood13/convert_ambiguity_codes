convert_ambiguity_codes
=======================

A simply python script to randomly convert ambiguity codes in a FASTA file.

A range of coalescent modelling methods require FASTA sequences generated from whole-genome requencing data. Heterozygouse base calls can be a bit of headache in these cases. So a simple way around this is simply to generate a consensus fasta file and then randomly choose one of the segregating alleles based on the ambiguity codes. 

This is exactly what this python script does. It's relatively fast and extremely simple. It will only accept FASTA files.
