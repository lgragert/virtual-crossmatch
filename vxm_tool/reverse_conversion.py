#! usr/bin/python 

import os
import re
import requests   
import hla
from hla import allele_truncate

ag_to_allele_dict = {}
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

UNOS_conversion_table_filename = "conversion_table.csv"
UNOS_conversion_table_file = open(UNOS_conversion_table_filename, 'r')
for row in UNOS_conversion_table_file:
	expression_character = ""
	if row.startswith("Allele"):
		continue 
	else:
		allele = row.split(',')[0]
		allele_4d = hla.allele_truncate(allele)
		antigen = row.split(',')[1]
		rule = row.split(',') [2]
		bw4_6 = row.split(',')[3]
		
		if antigen in ag_to_allele_dict.keys():
			if allele_4d in ag_to_allele_dict[antigen]:
				continue
			else:	

				ag_to_allele_dict[antigen].append(allele_4d)

		else:
			ag_to_allele_dict[antigen] = [allele_4d]	 


#print(ag_to_allele_dict)

final_dict = {}
#ag_list = []

#for ag in ag_to_allele_dict.keys():
	#ag_list.append(ag)
d = {}

for ag in ag_to_allele_dict.keys():
	allele_list = []
	d = {}
	if ag in UA_eq_dict.keys():
		ag_eqs = UA_eq_dict[ag]
		for ages in ag_eqs:
			ages = ages.strip()
			alleles = ag_to_allele_dict[ages]
			d[ages] = [alleles]

			#allele_list.append(alleles)
		
		#allele_list = [item for sublist in allele_list for item in sublist]
		#allele_list = allele_list.sort()
		#print(type(allele_list))
		final_dict[ag] = d

	else:
		final_dict[ag] = ag_to_allele_dict[ag]		
#print(final_dict)

def map_single_ag_to_alleles(antigen):
	allele_list = {}
	allele_only_list = []
	if antigen in final_dict:
		allele_list = final_dict[antigen]
	
	for i in allele_list.values():
		allele_only_list.append(i)

	flat_list = [item for sublist in allele_only_list for item in sublist]	
	flat_list = [item for sublist in flat_list for item in sublist]
	return flat_list
