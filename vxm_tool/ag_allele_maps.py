#! usr/bin/python

import os, re

import reverse_conversion 
from  reverse_conversion import map_single_ag_to_alleles, map_ag_for_proposed_algo


antigen_list = ["A2", "A26", "C05", "C01", "B72", "B27", "DR12", "DR1", "DQ7", "DQ5"] 

locus_ags_list = ["A2", "A26"]



def ags_to_strings(ag1, ag2):
	allele_list1 = map_ag_for_proposed_algo(ag1)
	allele_string1 = "/".join(allele_list1)
	allele_list2 = map_ag_for_proposed_algo(ag2)
	allele_string2 = "/".join(allele_list2)

	genotype_list = allele_string1 + "+" + allele_string2
	return genotype_list





def vxm_proposed_for_uags(donorantigens, donorRace, candidateUAs):
	conflicts = []
	conflict_ag_probs = {}
	donor_allele_freqs = {}
	ag_probs = {}
	bw_prob










# Donor Typing DOnor Alleles with probabilities Conflicting Antigens Candidate's Unacceptab	