#! usr/bin/python 

import os
import re
import requests   
import operator
import glob
import vxm_hla
from vxm_hla import allele_truncate, locus_string_geno_list, expand_ac, single_locus_allele_codes_genotype
import decimal
from decimal import Decimal

allele_to_ag_dict = {}
population_allele_frequencies = {}
allele_frequencies = {}
b_bw_dict = {}
bw4_list = []
bw6_list = []
ag_to_allele_dict = {}

unosEQags = []
agbw46 = {}
### Dictionary with alleles and equivalent antigens

UNOS_conversion_table_filename = "UNOS_conversion_table_with_rules.csv"
UNOS_conversion_table_file = open(UNOS_conversion_table_filename, 'r')
for row in UNOS_conversion_table_file:
	expression_character = ""
	if row.startswith("Allele"):
		continue 
	else:
		allele = row.split(',')[0]
		allele_4d = vxm_hla.allele_truncate(allele)
		antigen = row.split(',')[1]
		unosEQags.append(antigen)
		rule = row.split(',') [2]
		bw4_6 = row.split(',')[3]

		if bw4_6 != "NA":
			agbw46[antigen] = bw4_6
		if bw4_6 == "Bw4":
			bw4_list.append(antigen)
		if bw4_6 == "Bw6":
			bw6_list.append(antigen)

		if antigen in ag_to_allele_dict.keys():
			if allele_4d in ag_to_allele_dict[antigen]:
				continue
			else:
				ag_to_allele_dict[antigen].append(allele_4d)

		else:
			ag_to_allele_dict[antigen] = [allele_4d]	 

		
	allele_to_ag_dict[allele] = antigen, rule, bw4_6 
	allele_to_ag_dict[allele_4d] = antigen, rule, bw4_6

bw4_list = list(set(bw4_list))
bw6_list = list(set(bw6_list))

unosagslist = list(set(unosEQags))

b_bw_dict["Bw4"] = bw4_list
b_bw_dict["Bw6"] = bw6_list
	

for file in glob.glob('*.freqs'):
	#print(file)
	pop = file.split(".")[0]
	#print(pop)
	population_allele_frequencies[pop] = {}
	freq_file = open(file, 'r')
	for line in freq_file:
		if line.startswith("Haplo"):
			continue
		else:
			
			line_split = line.split(",")
			allele_list = line_split[0]
			count = line_split[1]
			haplotype_frequency = line_split[2]
			allele_split = allele_list.split("~")
			for allele in allele_split:
				allele = allele.rstrip("g")
				if allele in population_allele_frequencies[pop]:
					population_allele_frequencies[pop][allele] += float(haplotype_frequency)
				else:
					population_allele_frequencies[pop][allele] = float(haplotype_frequency)	

#print(population_allele_frequencies)


def map_ag_to_bw46(ag):
	if ag in agbw46.keys():
		bwe = agbw46[ag]

	else:
		bwe = ""	

def map_ag_to_alleles_pop_specific(antigen, pop):
	alleles = ag_to_allele_dict[antigen]
	common_alleles_in_pop = []
	
	for allele in alleles:
		if allele in population_allele_frequencies[pop]:
			common_alleles_in_pop.append(allele)


	pop_specific_alleles = list(set(common_alleles_in_pop))
	return pop_specific_alleles

def ags_to_strings(ag1, ag2, pop):
	allele_list1 = map_ag_to_alleles_pop_specific(ag1, pop)
	#print(len(allele_list1))
	allele_string1 = "/".join(allele_list1)
	allele_list2 = map_ag_to_alleles_pop_specific(ag2, pop)
	#print(len(allele_list2))x

	allele_string2 = "/".join(allele_list2)

	genotype_list = allele_string1 + "+" + allele_string2
	return genotype_list

def convert_allele_list_to_ags(hla_allele_list):
	
	"""This function can be called if a list of alleles has to be converted to antigens. Input format is a list and 
	the corresponding antigens and rules will be printed out"""
	allele_list_dict = {}
	ag_list = []
	bw4_6_list = []
	for allele in hla_allele_list:
		allele = allele.rstrip("p P g G")
		if allele in allele_to_ag_dict:
			ag = allele_to_ag_dict[allele][0]
			ag_list.append(ag)
			bw4_6 = allele_to_ag_dict[allele][2]
			if bw4_6 != "NA":
				bw4_6_list.append(bw4_6)
			
		else:
			continue
	
	bw46_list = list(set(bw4_6_list))
	ags = ag_list + bw46_list
				
	return ags

def gl_string_ags(gl_string, pop):
	gl_dict = {}
	gl_string = gl_string.replace("HLA-", "")
	locus_split = gl_string.split("^")

	ag_freq_1 = 0.0
	ag_freq_2 = 0.0
	geno_antigen_freq = {}

	ag_list = ""

	One_locus_typing = 0
	Two_locus_typing = 0
	Three_locus_typing = 0
	Four_locus_typing = 0
	Five_locus_typing = 0
	Six_locus_typing = 0

	

	if len(locus_split) == 1:
		One_locus_typing = 1
		#print("One locus typing")
		geno_antigen_freq = {}
		a_locus = locus_split[0]
		a_genotype_list = vxm_hla.locus_string_geno_list(a_locus)
		#print(a_genotype_list)
		a_ags = genotype_ags(a_genotype_list,pop)
		a_alleles = genotype_alleles(a_genotype_list, pop)
		#print(a_alleles)
		ag_list = a_ags  
		allele_list = a_alleles




	if len(locus_split) == 2:
		Two_locus_typing = 1
		print("Two locus typing")
		geno_antigen_freq = {}
		a_locus = locus_split[0]
		a_genotype_list = vxm_hla.locus_string_geno_list(a_locus)
		a_ags = genotype_ags(a_genotype_list,pop)
		a_alleles = genotype_alleles(a_genotype_list, pop)
		geno_antigen_freq = {}
		b_locus = locus_split[1]
		b_genotype_list = vxm_hla.locus_string_geno_list(b_locus)
		b_ags = genotype_ags(b_genotype_list,pop)
		b_alleles = genotype_alleles(b_genotype_list, pop)
		ag_list = a_ags  + b_ags 
		allele_list = a_alleles + b_alleles




	if len(locus_split) == 3:
		Three_locus_typing = 1
		print("Three locus typing")
		geno_antigen_freq = {}
		a_locus = locus_split[0]
		a_genotype_list = vxm_hla.locus_string_geno_list(a_locus)
		a_ags = genotype_ags(a_genotype_list,pop)
		a_alleles = genotype_alleles(a_genotype_list, pop)
		geno_antigen_freq = {}
		b_locus = locus_split[1]
		b_genotype_list = vxm_hla.locus_string_geno_list(b_locus)
		b_ags = genotype_ags(b_genotype_list,pop)
		b_alleles = genotype_alleles(b_genotype_list, pop)
		geno_antigen_freq = {}
		c_locus = locus_split[2]
		c_genotype_list = vxm_hla.locus_string_geno_list(c_locus)
		c_ags = genotype_ags(c_genotype_list,pop)
		c_alleles = genotype_alleles(c_genotype_list, pop)
		ag_list = a_ags  + b_ags + c_ags
		allele_list = a_alleles + b_alleles + c_alleles
	
	
	if len(locus_split) == 4:
		Four_locus_typing = 1
		print("Four locus typing")
		geno_antigen_freq = {}
		a_locus = locus_split[0]
		a_genotype_list = vxm_hla.locus_string_geno_list(a_locus)
		a_ags = genotype_ags(a_genotype_list,pop)
		a_alleles = genotype_alleles(a_genotype_list, pop)
		geno_antigen_freq = {}
		b_locus = locus_split[1]
		b_genotype_list = vxm_hla.locus_string_geno_list(b_locus)
		b_ags = genotype_ags(b_genotype_list,pop)
		b_alleles = genotype_alleles(b_genotype_list, pop)
		geno_antigen_freq = {}
		c_locus = locus_split[2]
		c_genotype_list = vxm_hla.locus_string_geno_list(c_locus)
		c_ags = genotype_ags(c_genotype_list,pop)
		c_alleles = genotype_alleles(c_genotype_list, pop)
		geno_antigen_freq = {}
		dr_locus = locus_split[3]
		dr_genotype_list = vxm_hla.locus_string_geno_list(dr_locus)
		dr_ags = genotype_ags(dr_genotype_list,pop)
		dr_alleles = genotype_alleles(dr_genotype_list, pop)
		ag_list = a_ags  + b_ags  + c_ags  + dr_ags
		allele_list = a_alleles + b_alleles + c_alleles + dr_alleles

	if len(locus_split) == 5:
		Five_locus_typing = 1
		print("Five locus typing")
		geno_antigen_freq = {}
		a_locus = locus_split[0]
		a_genotype_list = vxm_hla.locus_string_geno_list(a_locus)
		a_ags = genotype_ags(a_genotype_list,pop)
		a_alleles = genotype_alleles(a_genotype_list, pop)
		geno_antigen_freq = {}
		b_locus = locus_split[1]
		b_genotype_list = vxm_hla.locus_string_geno_list(b_locus)
		b_ags = genotype_ags(b_genotype_list,pop)
		b_alleles = genotype_alleles(b_genotype_list, pop)
		geno_antigen_freq = {}
		c_locus = locus_split[2]
		c_genotype_list = vxm_hla.locus_string_geno_list(c_locus)
		c_ags = genotype_ags(c_genotype_list,pop)
		c_alleles = genotype_alleles(c_genotype_list, pop)
		geno_antigen_freq = {}
		dr_locus = locus_split[3]
		dr_genotype_list = vxm_hla.locus_string_geno_list(dr_locus)
		dr_ags = genotype_ags(dr_genotype_list,pop)
		dr_alleles = genotype_alleles(dr_genotype_list, pop)
		geno_antigen_freq = {}
		dqb_locus = locus_split[4]
		dqb_genotype_list = vxm_hla.locus_string_geno_list(dqb_locus)
		dqb_ags = genotype_ags(dqb_genotype_list,pop)
		dq_alleles = genotype_alleles(dqb_genotype_list, pop)
		ag_list = a_ags  + b_ags  + c_ags   + dr_ags  + dqb_ags
		allele_list = a_alleles + b_alleles + c_alleles + dr_alleles + dq_alleles

	

	if len(locus_split) == 6:
		Six_locus_typing = 1
		print("Six locus typing")
		geno_antigen_freq = {}
		a_locus = locus_split[0]
		a_genotype_list = vxm_hla.locus_string_geno_list(a_locus)
		a_ags = genotype_ags(a_genotype_list,pop)
		a_alleles = genotype_alleles(a_genotype_list, pop)
		geno_antigen_freq = {}
		b_locus = locus_split[1]
		b_genotype_list = vxm_hla.locus_string_geno_list(b_locus)
		b_ags = genotype_ags(b_genotype_list,pop)
		b_alleles = genotype_alleles(b_genotype_list, pop)
		geno_antigen_freq = {}
		c_locus = locus_split[2]
		c_genotype_list = vxm_hla.locus_string_geno_list(c_locus)
		c_ags = genotype_ags(c_genotype_list,pop)
		c_alleles = genotype_alleles(c_genotype_list, pop)
		geno_antigen_freq = {}
		dr_locus = locus_split[3]
		dr_genotype_list = vxm_hla.locus_string_geno_list(dr_locus)
		dr_ags = genotype_ags(dr_genotype_list,pop)
		dr_alleles = genotype_alleles(dr_genotype_list, pop)
		geno_antigen_freq = {}
		dqb_locus = locus_split[4]
		dqb_genotype_list = vxm_hla.locus_string_geno_list(dqb_locus)
		dqb_ags = genotype_ags(dqb_genotype_list,pop)
		dq_alleles = genotype_alleles(dqb_genotype_list, pop)
		geno_antigen_freq = {}
		dr345_locus = locus_split[5]
		dr345_genotype_list = vxm_hla.locus_string_geno_list(dr345_locus)
		dr345_ags = genotype_ags(dr345_genotype_list,pop)
		dr345_alleles = genotype_alleles(dr345_genotype_list, pop)
		ag_list = a_ags + b_ags  + c_ags   + dr_ags  + dqb_ags  + dr345_ags
		allele_list = a_alleles + b_alleles + c_alleles + dr_alleles + dq_alleles + dr345_alleles
	#print(ag_list)
	#ages = ag_list[0::3]
	#bw46_list = ag_list[1::3]
	#probs = ag_list[2::3]
	#gl_dict = {"GL_string" : gl_string, "UNOS antigens": ages, "Bw4/6 epitopes": bw46_list, "Antigen Probablities": probs }
	#print(ag_list)
	return ag_list, allele_list

	
def genotype_ags(genotype_list, pop):
	ag_freq_1 = 0.0
	ag_freq_2 = 0.0
	
	geno_antigen_freq = {}
	for genotype in genotype_list:
		allele_1 = genotype.split("+")[0]
		allele_1 = allele_1.rstrip("g p P G")
		allele_1 = vxm_hla.allele_truncate(allele_1)

		allele_2 = genotype.split("+")[1]
		allele_2 = allele_2.rstrip("g p P G")
		allele_2 = vxm_hla.allele_truncate(allele_2)

		if allele_1 in allele_to_ag_dict.keys():
			ag_1 = allele_to_ag_dict[allele_1][0]
			bw46_1 = allele_to_ag_dict[allele_1][2]

		if allele_2 in allele_to_ag_dict.keys():	
			ag_2 = allele_to_ag_dict[allele_2][0]
			bw46_2 = allele_to_ag_dict[allele_2][2]
		


		if allele_1 in population_allele_frequencies[pop]:
			ag_freq_1 = population_allele_frequencies[pop][allele_1]
		if allele_2 in population_allele_frequencies[pop]:
			ag_freq_2 = population_allele_frequencies[pop][allele_2]

		gf = 0
		if (ag_1 == ag_2):
			gf = float(ag_freq_1) * float(ag_freq_2)
		else:
			gf = 2 * float(ag_freq_1) * float(ag_freq_2)	

		geno_antigen = ag_1 + "+" + ag_2
		
		if bw46_1 != "NA":

			geno_antigen = geno_antigen + "+" + bw46_1 

		if bw46_2	!= "NA":
			geno_antigen = geno_antigen + "+" + bw46_2	
			
		#print(geno_antigen)

		if geno_antigen in geno_antigen_freq.keys():
			geno_antigen_freq[geno_antigen] += float(gf)
		else:
			geno_antigen_freq[geno_antigen] = float(gf)

	TF = sum(geno_antigen_freq.values())
	if TF == 0.0:
		TF = 1
	else:
		TF = TF	

		
	for i,j in geno_antigen_freq.items():
		ag_probs = j/TF
		geno_antigen_freq[i] = ag_probs

		
	#print(geno_antigen_freq)		
	sorted_gf = sorted(geno_antigen_freq.items(), key = operator.itemgetter(1), reverse = True)
	#print(sorted_gf)
	#if len(sorted_gf) == 1:
		#ag_prob = 1
	#else:
		#antigen_list = sortef_gf[1::2]	
	#print(sorted_gf)
	# top_ag_geno = sorted_gf[0][0]
	# top_gf = sorted_gf[0][1]
	#print(top_ag_geno)	
	# ag_1 = top_ag_geno.split("+")[0]
	# ag_2 = top_ag_geno.split("+")[1]	
	#print(ag_1)
	#print(ag_2)
	# ag_list = ag_1 + "," + ag_2
	#bw46_list = bw46_1	+ "," + bw46_2
	return (sorted_gf)





def genotype_alleles(genotype_list, pop):
	allele_freq_1 = 0.0
	allele_freq_2 = 0.0
	
	geno_allele_freq = {}
	for genotype in genotype_list:
		allele_1 = genotype.split("+")[0]
		allele_1 = allele_1.rstrip("g p P G")
		allele_1 = vxm_hla.allele_truncate(allele_1)

		allele_2 = genotype.split("+")[1]
		allele_2 = allele_2.rstrip("g p P G")
		allele_2 = vxm_hla.allele_truncate(allele_2)

		
		if allele_1 in population_allele_frequencies[pop]:
			allele_freq_1 = population_allele_frequencies[pop][allele_1]
		if allele_2 in population_allele_frequencies[pop]:
			allele_freq_2 = population_allele_frequencies[pop][allele_2]

		gf = 0
		if (allele_1 == allele_2):
			gf = float(allele_freq_1) * float(allele_freq_2)
		else:
			gf = 2 * float(allele_freq_1) * float(allele_freq_2)	

		geno_allele = allele_1 + "+" + allele_2
	

		if geno_allele in geno_allele_freq.keys():
			geno_allele_freq[geno_allele] += float(gf)
		else:
			geno_allele_freq[geno_allele] = float(gf)

	TF = sum(geno_allele_freq.values())
	if TF == 0.0:
		TF = 1
	else:
		TF = TF	

		
	for i,j in geno_allele_freq.items():
		allele_probs = j/TF
		geno_allele_freq[i] = allele_probs

		
	#print(geno_antigen_freq)		
	sorted_gf = sorted(geno_allele_freq.items(), key = operator.itemgetter(1), reverse = True)
	#print(sorted_gf)
	
	return (sorted_gf)


def allele_freq(allele_list, pop):
	allele_pop_freqs = {}
	for i in allele_list:
		if i in population_allele_frequencies[pop]:
			allele_pop_freqs[i] = round(population_allele_frequencies[pop][i],4)
		else: 
			allele_pop_freqs[i] = 0	

	return allele_pop_freqs	


def allele_code_ags(allele_codes_list, pop):

	ag_freq_1 = 0.0
	ag_freq_2 = 0.0
	geno_antigen_freq = {}
	ag_list = ""

	
	if len(allele_codes_list) == 2:
		print("One locus typing")
		One_locus_typing = 1
		A_1_code = allele_codes_list[0]
		A_2_code = allele_codes_list[1]
		A_codes_pair = [A_1_code, A_2_code]
		geno_antigen_freq = {}
		acodes_genotype = vxm_hla.single_locus_allele_codes_genotype(A_codes_pair)
		a_ags = genotype_ags(acodes_genotype, pop)
		a_alleles = genotype_alleles(acodes_genotype, pop)
		ag_list = a_ags  
		allele_list = a_alleles


	if len(allele_codes_list) == 4:
		print("Two locus typing")
		Two_locus_typing = 1
		A_1_code = allele_codes_list[0]
		A_2_code = allele_codes_list[1]
		A_codes_pair = [A_1_code, A_2_code]
		geno_antigen_freq = {}
		acodes_genotype = vxm_hla.single_locus_allele_codes_genotype(A_codes_pair)
		a_ags = genotype_ags(acodes_genotype, pop)
		a_alleles = genotype_alleles(acodes_genotype, pop)

		C_1_code = allele_codes_list[2]
		C_2_code = allele_codes_list[3]
		C_codes_pair = [C_1_code, C_2_code]
		geno_antigen_freq = {}
		ccodes_genotype = vxm_hla.single_locus_allele_codes_genotype(C_codes_pair)
		c_ags = genotype_ags(ccodes_genotype, pop)
		c_alleles = genotype_alleles(ccodes_genotype, pop)

		ag_list = a_ags  + c_ags  

		allele_list = a_alleles + c_alleles


	if len(allele_codes_list) == 6:
		print("Three locus typing")
		Three_locus_typing = 1
		A_1_code = allele_codes_list[0]
		A_2_code = allele_codes_list[1]
		A_codes_pair = [A_1_code, A_2_code]
		geno_antigen_freq = {}
		acodes_genotype = vxm_hla.single_locus_allele_codes_genotype(A_codes_pair)
		a_ags = genotype_ags(acodes_genotype, pop)
		a_alleles = genotype_alleles(acodes_genotype, pop)

		C_1_code = allele_codes_list[2]
		C_2_code = allele_codes_list[3]
		C_codes_pair = [C_1_code, C_2_code]
		geno_antigen_freq = {}
		ccodes_genotype = vxm_hla.single_locus_allele_codes_genotype(C_codes_pair)
		c_ags = genotype_ags(ccodes_genotype, pop)
		c_alleles = genotype_alleles(ccodes_genotype, pop)

		B_1_code = allele_codes_list[4]
		B_2_code = allele_codes_list[5]
		B_codes_pair = [B_1_code, B_2_code]
		geno_antigen_freq = {}
		bcodes_genotype = vxm_hla.single_locus_allele_codes_genotype(B_codes_pair)
		b_ags = genotype_ags(bcodes_genotype, pop)
		b_alleles = genotype_alleles(bcodes_genotype, pop)
		
		ag_list = a_ags  + c_ags + b_ags

		allele_list = a_alleles + c_alleles + b_alleles

	if len(allele_codes_list) == 8:
		print("Four locus typing")
		Four_locus_typing = 1
		A_1_code = allele_codes_list[0]
		A_2_code = allele_codes_list[1]
		A_codes_pair = [A_1_code, A_2_code]
		geno_antigen_freq = {}
		acodes_genotype = vxm_hla.single_locus_allele_codes_genotype(A_codes_pair)
		a_ags = genotype_ags(acodes_genotype, pop)
		a_alleles = genotype_alleles(acodes_genotype, pop)

		C_1_code = allele_codes_list[2]
		C_2_code = allele_codes_list[3]
		C_codes_pair = [C_1_code, C_2_code]
		geno_antigen_freq = {}
		ccodes_genotype = vxm_hla.single_locus_allele_codes_genotype(C_codes_pair)
		c_ags = genotype_ags(ccodes_genotype, pop)
		c_alleles = genotype_alleles(ccodes_genotype, pop)

		B_1_code = allele_codes_list[4]
		B_2_code = allele_codes_list[5]
		B_codes_pair = [B_1_code, B_2_code]
		geno_antigen_freq = {}
		bcodes_genotype = vxm_hla.single_locus_allele_codes_genotype(B_codes_pair)
		b_ags = genotype_ags(bcodes_genotype, pop)
		b_alleles = genotype_alleles(bcodes_genotype, pop)
			
		dr_1_code = allele_codes_list[6]
		dr_2_code = allele_codes_list[7]
		dr_codes_pair = [dr_1_code, dr_2_code]
		geno_antigen_freq = {}
		drcodes_genotype = vxm_hla.single_locus_allele_codes_genotype(dr_codes_pair)
		dr_ags = genotype_ags(drcodes_genotype, pop)
		dr_alleles = genotype_alleles(drcodes_genotype, pop)

		ag_list = a_ags   + c_ags  + b_ags + dr_ags

		allele_list = a_alleles + c_alleles + b_alleles + dr_alleles
	
	if len(allele_codes_list) == 10:
		print("Five locus typing")
		Five_locus_typing = 1
		A_1_code = allele_codes_list[0]
		A_2_code = allele_codes_list[1]
		A_codes_pair = [A_1_code, A_2_code]
		geno_antigen_freq = {}
		acodes_genotype = vxm_hla.single_locus_allele_codes_genotype(A_codes_pair)
		a_ags = genotype_ags(acodes_genotype, pop)
		a_alleles = genotype_alleles(acodes_genotype, pop)

		C_1_code = allele_codes_list[2]
		C_2_code = allele_codes_list[3]
		C_codes_pair = [C_1_code, C_2_code]
		geno_antigen_freq = {}
		ccodes_genotype = vxm_hla.single_locus_allele_codes_genotype(C_codes_pair)
		c_ags = genotype_ags(ccodes_genotype, pop)
		c_alleles = genotype_alleles(ccodes_genotype, pop)

		B_1_code = allele_codes_list[4]
		B_2_code = allele_codes_list[5]
		B_codes_pair = [B_1_code, B_2_code]
		geno_antigen_freq = {}
		bcodes_genotype = vxm_hla.single_locus_allele_codes_genotype(B_codes_pair)
		b_ags = genotype_ags(bcodes_genotype, pop)
		b_alleles = genotype_alleles(bcodes_genotype, pop)
			
		dr_1_code = allele_codes_list[6]
		dr_2_code = allele_codes_list[7]
		dr_codes_pair = [dr_1_code, dr_2_code]
		geno_antigen_freq = {}
		drcodes_genotype = vxm_hla.single_locus_allele_codes_genotype(dr_codes_pair)
		dr_ags = genotype_ags(drcodes_genotype, pop)	
		dr_alleles = genotype_alleles(drcodes_genotype, pop)
		
		dqb_1_code = allele_codes_list[8]
		dqb_2_code = allele_codes_list[9]
		dqb_codes_pair = [dqb_1_code, dqb_2_code]
		geno_antigen_freq = {}
		dqbcodes_genotype = vxm_hla.single_locus_allele_codes_genotype(dqb_codes_pair)
		dqb_ags = genotype_ags(dqbcodes_genotype, pop)
		dqb_alleles = genotype_alleles(dqbcodes_genotype, pop)

		ag_list = a_ags  + c_ags  + b_ags + dr_ags  + dqb_ags

		allele_list = a_alleles + c_alleles + b_alleles + dr_alleles + dqb_alleles
				
	if len(allele_codes_list) == 12:
		print("Six locus typing")
		Six_locus_typing = 1
		A_1_code = allele_codes_list[0]
		A_2_code = allele_codes_list[1]
		A_codes_pair = [A_1_code, A_2_code]
		geno_antigen_freq = {}
		acodes_genotype = vxm_hla.single_locus_allele_codes_genotype(A_codes_pair)
		a_ags = genotype_ags(acodes_genotype, pop)
		a_alleles = genotype_alleles(acodes_genotype, pop)


		C_1_code = allele_codes_list[2]
		C_2_code = allele_codes_list[3]
		C_codes_pair = [C_1_code, C_2_code]
		geno_antigen_freq = {}
		ccodes_genotype = vxm_hla.single_locus_allele_codes_genotype(C_codes_pair)
		c_ags = genotype_ags(ccodes_genotype, pop)
		c_alleles = genotype_alleles(ccodes_genotype, pop)

		B_1_code = allele_codes_list[4]
		B_2_code = allele_codes_list[5]
		B_codes_pair = [B_1_code, B_2_code]
		geno_antigen_freq = {}
		bcodes_genotype = vxm_hla.single_locus_allele_codes_genotype(B_codes_pair)
		b_ags = genotype_ags(bcodes_genotype, pop)
		b_alleles = genotype_alleles(bcodes_genotype, pop)
			
		dr_1_code = allele_codes_list[6]
		dr_2_code = allele_codes_list[7]
		dr_codes_pair = [dr_1_code, dr_2_code]
		geno_antigen_freq = {}
		drcodes_genotype = vxm_hla.single_locus_allele_codes_genotype(dr_codes_pair)
		dr_ags = genotype_ags(drcodes_genotype, pop)
		dr_alleles = genotype_alleles(drcodes_genotype, pop)	
		
		dqb_1_code = allele_codes_list[8]
		dqb_2_code = allele_codes_list[9]
		dqb_codes_pair = [dqb_1_code, dqb_2_code]
		geno_antigen_freq = {}
		dqbcodes_genotype = vxm_hla.single_locus_allele_codes_genotype(dqb_codes_pair)
		dqb_ags = genotype_ags(dqbcodes_genotype, pop)
		dqb_alleles = genotype_alleles(dqbcodes_genotype, pop)
							

		dr345_1_code = allele_codes_list[10]
		dr345_2_code = allele_codes_list[11]
		dr345_codes_pair = [dr345_1_code, dr345_2_code]
		geno_antigen_freq = {}
		dr345codes_genotype = vxm_hla.single_locus_allele_codes_genotype(dr345_codes_pair)
		dr345_ags = genotype_ags(dr345codes_genotype, pop)
		dr345_alleles = genotype_alleles(dr345codes_genotype, pop)			


		ag_list = a_ags  + c_ags  + b_ags + dr_ags  + dqb_ags  + dr345_ags

		allele_list = a_alleles + c_alleles + b_alleles + dr_alleles + dqb_alleles + dr345_alleles


	#ages = ag_list[0::3]
	#bw46_list = ag_list[1::3]
	#probs = ag_list[2::3]
	#al_dict = {"Allele Codes" : allele_codes_list, "Antigens": ages, "Bw4/6 epitopes": bw46_list, "Antigen Probablities": probs}
	#print(ag_list)
	return ag_list, allele_list
	
def convert_ag_list_to_gls(ag_list, pop):
	if len(ag_list) == 2:
		print("One Locus Typing")
		ag1 = ag_list[0]
		ag2 = ag_list[1]
		gls = ags_to_strings(ag1, ag2, pop)

	if len(ag_list) == 4:
		print("Two Locus Typing")
		ag1 = ag_list[0]
		ag2 = ag_list[1]
		gls1 = ags_to_strings(ag1, ag2, pop)

		ag3 = ag_list[2]
		ag4 = ag_list[3]
		gls2 = ags_to_strings(ag3, ag4, pop)

		gls = gls1 + "^" + gls2

	if len(ag_list) == 6:
		print("Three Locus Typing")
		ag1 = ag_list[0]
		ag2 = ag_list[1]
		gls1 = ags_to_strings(ag1, ag2, pop)

		ag3 = ag_list[2]
		ag4 = ag_list[3]
		gls2 = ags_to_strings(ag3, ag4, pop)

		ag5 = ag_list[4]
		ag6 = ag_list[5]
		gls3 = ags_to_strings(ag5, ag6, pop)

		gls = gls1 + "^" + gls2 + "^" + gls3

	if len(ag_list) == 8:
		print("Four Locus Typing")
		ag1 = ag_list[0]
		ag2 = ag_list[1]
		gls1 = ags_to_strings(ag1, ag2, pop)
		
		ag3 = ag_list[2]
		ag4 = ag_list[3]
		gls2 = ags_to_strings(ag3, ag4, pop)
		
		ag5 = ag_list[4]
		ag6 = ag_list[5]
		gls3 = ags_to_strings(ag5, ag6, pop)
		
		
		ag7 = ag_list[6]
		ag8 = ag_list[7]
		gls4 =  ags_to_strings(ag7, ag8, pop)

		gls = gls1 + "^" + gls2 + "^" + gls3 + "^" + gls4 

	if len(ag_list) == 10:
		print("Five locus Typing")	

		ag1 = ag_list[0]
		ag2 = ag_list[1]
		gls1 = ags_to_strings(ag1, ag2, pop)
		
		ag3 = ag_list[2]
		ag4 = ag_list[3]
		gls2 = ags_to_strings(ag3, ag4, pop)
		
		ag5 = ag_list[4]
		ag6 = ag_list[5]
		gls3 = ags_to_strings(ag5, ag6, pop)
		
		
		ag7 = ag_list[6]
		ag8 = ag_list[7]
		gls4 =  ags_to_strings(ag7, ag8, pop)

		ag9 = ag_list[8]
		ag10 = ag_list[9]
		gls5 =  ags_to_strings(ag9, ag10, pop)

		gls = gls1 + "^" + gls2 + "^" + gls3 + "^" + gls4 + "^" + gls5



		
	if len(ag_list) == 12:
		print("Six locus Typing")	

		ag1 = ag_list[0]
		ag2 = ag_list[1]
		gls1 = ags_to_strings(ag1, ag2, pop)
		
		ag3 = ag_list[2]
		ag4 = ag_list[3]
		gls2 = ags_to_strings(ag3, ag4, pop)
		
		ag5 = ag_list[4]
		ag6 = ag_list[5]
		gls3 = ags_to_strings(ag5, ag6, pop)
		
		
		ag7 = ag_list[6]
		ag8 = ag_list[7]
		gls4 =  ags_to_strings(ag7, ag8, pop)

		ag9 = ag_list[8]
		ag10 = ag_list[9]
		gls5 =  ags_to_strings(ag9, ag10, pop)

		ag11 = ag_list[10]
		ag12 = ag_list[11]
		gls6 = ags_to_strings(ag11, ag12, pop)

		gls = gls1 + "^" + gls2 + "^" + gls3 + "^" + gls4 + "^" + gls5	+ "^" + gls6

	#print(gls)

	return gls

