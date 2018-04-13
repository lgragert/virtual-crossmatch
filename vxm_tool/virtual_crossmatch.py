#! usr/bin/python


import os, re
import requests
import hla 
from hla import allele_truncate, locus_string_geno_list, expand_ac, single_locus_allele_codes_genotype

import conversion_functions
from conversion_functions import convert_allele_to_ag, convert_allele_list_to_ags, gl_string_ags, genotype_ags, allele_code_ags


UA_eq_dict = {}



UNOS_UA_eq_filename = "UNOS_4-10_ag_equivalencies.csv"
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
#print(UA_eq_dict)



def vxm_gls(donor_gl_string, donor_ethnicity, recepient_UA_list):
	conflicts = []
	output = gl_string_ags(donor_gl_string, donor_ethnicity)
	
	ags = list(output['UNOS antigens']) ## returns tuple of strings
	ags_list = []

	for i in ags: 						### converts tuple of strings to list
		isplit = i.split(",")
		x = isplit[0]
		ags_list.append(x)
		y = isplit[1]
		ags_list.append(y)



	pub_epis = list(output['Bw4/6 epitopes'])
	epis_list = []
	for i in pub_epis:
		isplit = i.split(",")
		x = isplit[0]
		epis_list.append(x)
		y = isplit[1]
		epis_list.append(y)


	episn = [x for x in epis_list if x!="NA"]
	
	donor_ags = ags_list + episn

	#print("Donor Antigens")
	#print(donor_ags)

	UA_list = []
	for ag in recepient_UA_list:
		if ag in UA_eq_dict.keys():
			UA_list.append(UA_eq_dict[ag])
		else:
			UA_list.append([ag])	


	recepient_ags = [item for sublist in UA_list for item in sublist]
	#print("Recepient Unacceptable Antigens")
	#print(recepient_ags)

	for ag in donor_ags:
		if ag in recepient_ags:
			conflicts.append(ag)
		
			
	if len(conflicts) == 0:
		print("Virtual Crossmatch is negative")

	else:
		print("Virtual Crossmatch is positive and the overlapping antigens are following")	

	return (donor_ags, recepient_ags, conflicts)





def vxm_allele_codes(allele_codes_list, donor_ethnicity, recepient_UA_list):
	conflicts = []
	output = allele_code_ags(allele_codes_list, donor_ethnicity)
	
	ags = list(output['Antigens']) ## returns tuple of strings
	print(ags)
	ags_list = []

	for i in ags: 						### converts tuple of strings to list
		isplit = i.split(",")
		x = isplit[0]
		ags_list.append(x)
		y = isplit[1]
		ags_list.append(y)



	pub_epis = list(output['Bw4/6 epitopes'])
	epis_list = []
	for i in pub_epis:
		isplit = i.split(",")
		x = isplit[0]
		epis_list.append(x)
		y = isplit[1]
		epis_list.append(y)


	episn = [x for x in epis_list if x!="NA"]
	#print(episn)
	
	donor_ags = ags_list + episn

	#print("Donor Antigens")
	#print(donor_ags)

	UA_list = []
	for ag in recepient_UA_list:
		if ag in UA_eq_dict.keys():
			UA_list.append(UA_eq_dict[ag])
		else:
			UA_list.append([ag])	


	recepient_ags = [item for sublist in UA_list for item in sublist]
	#print("Recepient Unacceptable Antigens")
	#print(recepient_ags)

	for ag in donor_ags:
		if ag in recepient_ags:
			conflicts.append(ag)
		
	if len(conflicts) == 0:
		print("Virtual Crossmatch is negative")

	else:
		print("Virtual Crossmatch is positive and the overlapping antigens are following")	

	return (donor_ags, recepient_ags, conflicts)
