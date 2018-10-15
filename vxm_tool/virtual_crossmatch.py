#! usr/bin/python


import os, re
import requests
import vxm_hla 
import itertools
from vxm_hla import allele_truncate, locus_string_geno_list, expand_ac, single_locus_allele_codes_genotype, gl_string_alleles_list, allele_code_to_allele_list

import conversion_functions_for_VXM
from conversion_functions_for_VXM import  gl_string_ags, genotype_ags, allele_code_ags, unosagslist, convert_allele_list_to_ags, allele_freq

import reverse_conversion

from reverse_conversion import map_single_ag_to_alleles

UA_eq_dict = {}



UNOS_UA_eq_filename = "UNOS_UA_ag_equivalencies.csv"
UNOS_UA_eq_file = open(UNOS_UA_eq_filename, 'r')

for row in UNOS_UA_eq_file:
	if row.startswith("Antigen"):
		continue
	else:
		row = row.strip("\n")
		row_split = row.split(",")
		ua_ag = row_split[0]
		ua_ag_eqs = row_split[1:]
		ua_ag_eqs = list(filter(None, ua_ag_eqs))
		UA_eq_dict[ua_ag] = ua_ag_eqs
print(UA_eq_dict)


def vxm_uags(donorags, candidateags):
	conflicts = []
	donor_ags_alleles = []
	UA_list = []
	for ag in candidateags:
		if ag in UA_eq_dict.keys():
			UA_list.append(UA_eq_dict[ag])
		else:
			UA_list.append([ag])	


	recepient_ags = [item for sublist in UA_list for item in sublist]


	for ag in donorags:
		alleles = map_single_ag_to_alleles(ag)
		donor_ags_alleles.append([ag])
		if alleles:
			donor_ags_alleles.append(alleles)
	
	merged_dags_alleles = list(itertools.chain(*donor_ags_alleles))

	for ag in recepient_ags:
		if ag in merged_dags_alleles:
			conflicts.append(ag)

	return (donorags, recepient_ags, conflicts)


def vxm_hIresalleles(donorsAlleleList, candidateags):
	conflicts = []
	donorags = []

	UA_list = []
	

	for ag in candidateags:
		if ag in UA_eq_dict.keys():
			UA_list.append(UA_eq_dict[ag])
		else:
			UA_list.append([ag])	


	recepient_ags = [item for sublist in UA_list for item in sublist]

	donorags = convert_allele_list_to_ags(donorsAlleleList)

	donorags_alleles = donorsAlleleList + donorags
	
	for ag in recepient_ags:
		if ag in donorags_alleles:
			conflicts.append(ag)

	return(donorags, recepient_ags, conflicts)


def vxm_gls(donor_gl_string, donor_ethnicity, recipient_UA_list):
	conflicts = []
	ag_probs = {}
	donor_ags = []
	output = gl_string_ags(donor_gl_string, donor_ethnicity)
	#print("Output")
	#print(output)
	
	for i in output:
		ag_list = i[0].split("+")
		for j in ag_list:
			if j in ag_probs.keys():
				ag_probs[j] += i[1]
			else:
				ag_probs[j] = i[1]
	
	##################################################################################################

	for ag in ag_probs.keys():
		if ag in UA_dict.keys():
			p=0
			for eqs in UA_dict[ag]:
				if eqs in ag_probs.keys():
					p += ag_probs[eqs]
					ag_probs[ag] = p
		else:
			continue

	#############################################################################################
	##print(ag_probs)
	donor_alleles = vxm_hla.gl_string_alleles_list(donor_gl_string)
	#print(donor_alleles)
	donor_allele_freqs = conversion_functions_for_VXM.allele_freq(donor_alleles, donor_ethnicity)
	#print(donor_allele_freqs)

	for k in ag_probs.keys():        #### gets all the mapped donor antigens as a list
		donor_ags.append(k)

	UA_list = []
	for ag in recipient_UA_list:           
		if ag in UA_eq_dict.keys():
			UA_list.append(UA_eq_dict[ag])
		else:
			UA_list.append([ag])	


	recepient_ags = [item for sublist in UA_list for item in sublist]            ### makes a list from list of lists


	donor_alleles_ags = donor_ags + donor_alleles 							######## list of donor alleles from GL string and the antigens mapped as per ALLAN
	
	for ag in recepient_ags:
		if ag in donor_alleles_ags:
			conflicts.append(ag)
		
	conflict_ag_probs = {}

	for i in conflicts:
		if i in ag_probs.keys():
			conflict_ag_probs[i] = round(ag_probs[i], 4)

		elif i in donor_allele_freqs.keys():
			conflict_ag_probs[i] = round(donor_allele_freqs[i], 4)

		else:
			conflict_ag_probs[i] = 0	
	
	for i,j in conflict_ag_probs.items():
		if j > 1.00:
			j = 1.00
			conflict_ag_probs[i] = j

	#print(conflict_ag_probs)


	return(donor_ags, recepient_ags, conflicts, conflict_ag_probs)








def vxm_allele_codes(allele_codes_list, donor_ethnicity, recepient_UA_list):
	conflicts = []
	ag_probs = {}
	donor_ags = []
	output = allele_code_ags(allele_codes_list, donor_ethnicity)

	
	for i in output:
		ag_list = i[0].split("+")
		for j in ag_list:
			if j in ag_probs.keys():
				ag_probs[j] += i[1]
			else:
				ag_probs[j] = i[1]
	#print(ag_probs)
	
	for ag in ag_probs.keys():
		if ag in UA_dict.keys():
			p=0
			for eqs in UA_dict[ag]:
				if eqs in ag_probs.keys():
					p += ag_probs[eqs]
					ag_probs[ag] = p
		else:
			continue
	
	donor_alleles = vxm_hla.allele_code_to_allele_list(allele_codes_list)
	donor_allele_freqs = conversion_functions_for_VXM.allele_freq(donor_alleles, donor_ethnicity)

	for k in ag_probs.keys():
		donor_ags.append(k)

	UA_list = []
	for ag in recepient_UA_list:
		if ag in UA_eq_dict.keys():
			UA_list.append(UA_eq_dict[ag])
		else:
			UA_list.append([ag])	


	recepient_ags = [item for sublist in UA_list for item in sublist]

	donor_alleles_ags = donor_ags + donor_alleles

	for ag in recepient_ags:
		if ag in donor_alleles_ags:
			conflicts.append(ag)
	
	conflict_ag_probs = {}

	
	for i in conflicts:
		if i in ag_probs.keys():
			conflict_ag_probs[i] = round(ag_probs[i], 4)

		elif i in donor_allele_freqs.keys():
			conflict_ag_probs[i] = round(donor_allele_freqs[i], 4)

		else:
			conflict_ag_probs[i] = 0



	for i,j in conflict_ag_probs.items():
		if j > 1.00:
			j = 1.00
			conflict_ag_probs[i] = j

	#print(conflict_ag_probs)


	return(donor_ags, recepient_ags, conflicts, conflict_ag_probs)





