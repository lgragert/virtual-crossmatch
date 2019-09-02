
#! usr/bin/python

##### This script has functions that map the VXM parameters in the high resolution VXM probability table
import os, re
import vxm_hla
import decimal
from decimal import Decimal

regex = re.compile(r'(\d+|\s+)')

#### Making dictionary of the alleles and antigens, bw4/6 from conversion table to map elements in the high res table
allele_to_ag_dict = {}
UNOS_conversion_table_filename = "UNOS_conversion_table_with_rules.csv"
#UNOS_conversion_table_filename = "/opt/bitnami/apps/django/django_projects/Project/conversion_table.csv"
UNOS_conversion_table_file = open(UNOS_conversion_table_filename, 'r')
for row in UNOS_conversion_table_file:
	expression_character = ""
	if row.startswith("Allele"):
		continue 
	else:
		allele = row.split(',')[0]
		allele_4d = vxm_hla.allele_truncate(allele)
		antigen = row.split(',')[1]
		bw4_6 = row.split(',')[3]
	
	allele_to_ag_dict[allele_4d] = [antigen, bw4_6]

regex = re.compile(r'(\d+|\s+)')



def group_serotypes_per_locus(ag_list):
	'''Function to group antigens and alleles as per loci (UNOS antigen equivalents, UAs and OPTN equiavalents)'''
	a_locus_ag = []
	b_locus_ag = []
	bw_locus_ag = []
	drb345_locus_ag = []
	c_locus_ag = []
	dr_locus_ag = []
	dq_locus_ag = []

	for ag in ag_list:
		if "*" in ag:
			locus = ag.split("*")[0]
			if locus == "A":
				a_locus_ag.append(ag)
			if locus == "B":
				b_locus_ag.append(ag)
			if locus == "C":
				c_locus_ag.append(ag)
			if locus == "DRB1":
				dr_locus_ag.append(ag)
			if locus == "DQB1":
				dq_locus_ag.append(ag)
			if locus == "DRB3" or locus == "DRB4" or locus == "DRB5":
				drb345_locus_ag.append(ag)

	
		elif ag == "Bw4" or ag == "Bw6":
			bw_locus_ag.append(ag)
		elif ag == "DR51" or ag == "DR52" or ag == "DR53":
			drb345_locus_ag.append(ag)	
	
		else:
			ag_locus = regex.split(ag)[0]
			if ag_locus == "A":
				a_locus_ag.append(ag)
			if ag_locus == "B":
				b_locus_ag.append(ag)
			if ag_locus == "C":
				c_locus_ag.append(ag)
			if ag_locus == "DR":
				dr_locus_ag.append(ag)	
			if ag_locus == "DQ":
				dq_locus_ag.append(ag)	

			

	locus_sorted_ag_list = [", ".join(sorted(a_locus_ag))] + [", ".join(sorted(b_locus_ag))] + [", ".join(sorted(bw_locus_ag))] + [", ".join(sorted(c_locus_ag))] + [", ".join(sorted(dr_locus_ag))] + [", ".join(sorted(drb345_locus_ag))] + [", ".join(sorted(dq_locus_ag))]


	return locus_sorted_ag_list	


def group_serotypes_per_locus_with_bw(sorted_alleles_dict, ag_list):
	'''Maps the antigens to alleles in each row of the high res table'''
	a_locus_ag = []
	b_locus_ag = []
	bw_locus_ag = []
	c_locus_ag = []
	dr_locus_ag = []
	drb345_locus_ag = []
	dq_locus_ag = []

	for tup in sorted_alleles_dict:
		i = tup[0]
		locus = i.split("*")[0]
		j = tup[1]
		ag = allele_to_ag_dict[i][0]
		bw = allele_to_ag_dict[i][1]
		if locus == "A":
			if ag in ag_list:
				a_locus_ag.append(ag)
			else:
				ag = ""	
				a_locus_ag.append(ag)

		if locus == "B":
			if ag in ag_list:
				b_locus_ag.append(ag)
			else:
				ag = ""	
				b_locus_ag.append(ag)
			
		if locus == "C":
			if ag in ag_list:
				c_locus_ag.append(ag)
			else:
				ag = ""	
				c_locus_ag.append(ag)	

		if locus == "DRB1":
			if ag in ag_list:
				dr_locus_ag.append(ag)
			else:
				ag = ""	
				dr_locus_ag.append(ag)

		if locus == "DRB3" or locus == "DRB4" or locus == "DRB5":
			if ag in ag_list:
				drb345_locus_ag.append(ag)
			else:
				ag = ""	
				drb345_locus_ag.append(ag)		
		
		if locus == "DQB1":
			if ag in ag_list:
				dq_locus_ag.append(ag)
			else:
				ag = ""	
				dq_locus_ag.append(ag)
	for i in ag_list:
		if i == "Bw4" or i == "Bw6":
			bw_locus_ag.append(i)	
		

	locus_sorted_ag_list = [a_locus_ag] + [b_locus_ag] + [bw_locus_ag] + [c_locus_ag] + [dr_locus_ag] + [drb345_locus_ag] + [dq_locus_ag]


	return locus_sorted_ag_list	

def split_gl_string_per_locus(gl_string, donor_bws_string):
	'''The GL string entered by user is chopped based on loci for distinct rows in the high res table'''
	a_string = ""
	b_string = ""
	bw_string = ""
	c_string = ""
	dr_string = ""
	dq_string = ""
	drb345_string = ""

	gl_string_split = gl_string.split("^")	
	

	for string in gl_string_split:
		string_locus = string.split("*")[0]

		if string_locus == "DRB3" or string_locus == "DRB4" or string_locus == "DRB5":
			drb345_string = string.replace("+", " + ") 
			drb345_string = drb345_string.replace("/", "/ ")
		
		if string_locus == "A":
			a_string = string.replace("+", " + ") 
			a_string = a_string.replace("/", "/ ")

		if string_locus == "B":
			b_string = string.replace("+", " + ") 
			b_string = b_string.replace("/", "/ ")

		if string_locus == "C":
			c_string = string.replace("+", " + ")
			c_string = c_string.replace("/", "/ ") 

		if string_locus == "DRB1":
			dr_string = string.replace("+", " + ") 
			dr_string = dr_string.replace("/", "/ ")

		if string_locus == "DQB1":
			dq_string = string.replace("+", " + ") 
			dq_string = dq_string.replace("/", "/ ") 

	bw_string = donor_bws_string

	string_list = [a_string] + [b_string] + [bw_string] +[c_string] + [dr_string] + [drb345_string] + [dq_string]	
	
	return string_list				


##### for alleles #################

def prob_dict_list_of_strings(allele_probs, bw_prob):
	'''The alleles from distinct loci are listed in different rows in high res table'''

	a_alleles = [] 
	bw_epitopes = []
	b_alleles = []
	c_alleles = []
	dr_alleles = []
	drb345_alleles = []
	dq_alleles = []


	for tup in allele_probs:
		i = tup[0]
		j = tup[1]
		locus = i.split("*")[0]
		if locus == "A":
			ji = round_freq_se(j)
			fi = i + ": " + ji 
			a_alleles.append(fi)

		if locus == "B":
			ji = round_freq_se(j)
			fi = i + ": " + ji 
			b_alleles.append(fi)


		if locus == "C":
			ji = round_freq_se(j)
			fi = i + ": " + ji 
			c_alleles.append(fi)

		if locus == "DQB1":
			ji = round_freq_se(j)
			fi = i + ": " + ji 
			dq_alleles.append(fi)       

		if locus == "DRB1":
			ji = round_freq_se(j)
			fi = i + ": " + ji 
			dr_alleles.append(fi)

		if locus == "DRB3" or locus == "DRB4" or locus == "DRB5":
			ji = round_freq_se(j)
			fi = i + ": " + ji 
			drb345_alleles.append(fi)	


	for tup2 in bw_prob:
		i = tup2[0]
		j = tup2[1]
		ji = round_freq_se(j)
		fi = i + ": " + ji
		bw_epitopes.append(fi)

	list_of_allele_probs = [a_alleles] + [b_alleles] + [sorted(bw_epitopes)] + [c_alleles] + [dr_alleles] + [drb345_alleles] + [dq_alleles] 

	return list_of_allele_probs


###### for antigens, conflicting antigens, optne equivalents

def prob_dict_list_of_strings_for_antigens(sorted_alleles_dict, ag_probs):
	'''Antigens are mapped to alleles in the second column in the high res table'''

	a_ags = [] 
	b_ags = []
	c_ags = []
	dr_ags = []
	drb345_ags = []
	dq_ags = []
	bw_ags = []


	for tup in sorted_alleles_dict:
		i = tup[0]
		j = tup[1]

		bw_prob = ""
		ag_prob = ""

		ag = allele_to_ag_dict[i][0]
		bw = allele_to_ag_dict[i][1]
		prob1 = ag_probs[ag]
		ag_prob = ag + ": " + round_freq_se(prob1)

		if bw != "NA":
			prob2 = ag_probs[bw]
			bw_prob = bw + ":" + round_freq_se(prob2)
			bw_ags.append(bw_prob)
		
		if ag == "DR51" or ag == "DR52" or ag == "DR53":
			drb345_ags.append(ag_prob)
		
		elif ag.startswith("A"):
			a_ags.append(ag_prob)

		elif ag.startswith("B"):
			b_ags.append(ag_prob)	


		elif ag.startswith("C"):
			c_ags.append(ag_prob)

		elif ag.startswith("DR"):
			dr_ags.append(ag_prob)

		elif ag.startswith("DQ"):
			dq_ags.append(ag_prob)				


	list_of_ag_probs = [a_ags] + [b_ags] + [sorted(list(set(bw_ags)))] +[c_ags] + [dr_ags] + [drb345_ags] +[dq_ags] 

	return list_of_ag_probs


def group_list_of_alleles_per_locus(donor_typing_list):
	'''The alleles are grouped according to loci'''
	donor_a_alleles = []
	donor_b_alleles = []
	donor_c_alleles = []
	donor_dr_alleles = []
	donor_dr345_alleles = []
	donor_dq_alleles = []
	

	for i in donor_typing_list:
		locus = i.split("*")[0]
		if locus == "A":
			donor_a_alleles.append(i)

		if locus == "B":
			donor_b_alleles.append(i)

		if locus == "C":
			donor_c_alleles.append(i)
		
		if (locus == "DRB3") or (locus == "DRB4") or (locus == "DRB5"):
			donor_dr345_alleles.append(i)
		
		if (locus == "DRB1"):
			donor_dr_alleles.append(i)

		if (locus == "DQB1") or (locus == "DQA1"):
			donor_dq_alleles.append(i)


	final_typing_list = [", ".join(sorted(donor_a_alleles))] + [", ".join(sorted(donor_b_alleles))] + [""] + [", ".join(sorted(donor_c_alleles))] + [", ".join(sorted(donor_dr_alleles))] + [", ".join(sorted(donor_dr345_alleles))] +[", ".join(sorted(donor_dq_alleles))] 

	return final_typing_list


#### Function to group allele codes as per locus

def group_allele_codes_per_locus(donor_typing_list, donor_bws_string):
	'''Allele codes are grouped according to loci'''
	donor_a_alleles = []
	donor_b_alleles = []
	donor_c_alleles = []
	donor_dr_alleles = []
	donor_drb345_alleles = []
	donor_dq_alleles = []
	

	for i in donor_typing_list:
		locus = i.split("*")[0]
		if locus == "A":
			donor_a_alleles.append(i)

		if locus == "B":
			donor_b_alleles.append(i)

		if locus == "C":
			donor_c_alleles.append(i)

		if (locus == "DRB1"): 
			donor_dr_alleles.append(i)

		if (locus == "DRB3") or (locus == "DRB4") or (locus == "DRB5"):	
			donor_drb345_alleles.append(i)

		if (locus == "DQB1") or (locus == "DQA1"):
			donor_dq_alleles.append(i)

	bw_string = donor_bws_string

	final_typing_list = [", ".join(sorted(donor_a_alleles))] + [", ".join(sorted(donor_b_alleles))] + [bw_string] + [", ".join(sorted(donor_c_alleles))] + [", ".join(sorted(donor_dr_alleles))] + [", ".join(sorted(donor_drb345_alleles))] + [", ".join(sorted(donor_dq_alleles))] 

	return final_typing_list



#### Grouping conflicting antigens as per locus and with their probabilities. These are DSA

def conflicts_ags(sorted_allele_dict, conflicts): 
	'''The conflicting antigens are mapped to alleles in the second column'''

	cag_prob_dict = {}

	for tup in sorted_allele_dict:
		i = tup[0]
		j = tup[1]
		ag = allele_to_ag_dict[i][0]
		if i in conflicts:
			#print(i)
			prob = conflicts[i]
			cag_prob_dict[i] = i + ": " + round_freq_se(prob)

		if ag in conflicts:
			prob = conflicts[ag]
			if i in cag_prob_dict.keys():
				cag_prob_dict[i] = cag_prob_dict[i] + " " + ("(" + ag + ": " + round_freq_se(prob)+ ")")

			else:
				cag_prob_dict[i] = 	ag + ": " + round_freq_se(prob)

		if i not in conflicts and ag not in conflicts:
			cag_prob_dict[i] = ""
	
	#print(cag_prob_dict)
	bw_cags = []
	
	for i,j in conflicts.items():
		if i == "Bw4" or i == "Bw6":
			bw_prob_list_element = i + ": " + round_freq_se(j)
			bw_cags.append(bw_prob_list_element)	

	a_cags = []
	b_cags = []
	c_cags = []
	dr_cags = []
	drb345_cags = []
	dq_cags = []
	

	for i,j in cag_prob_dict.items():
		locus = i.split("*")[0]

		if locus == "A":
			a_cags.append(j)

		if locus == "B":
			b_cags.append(j)

		if locus == "C":
			c_cags.append(j)

		if locus == "DRB1":
			dr_cags.append(j)

		if locus == "DRB3" or locus == "DRB4" or locus == "DRB5": 
			drb345_cags.append(j)	

		if locus == "DQB1":
			dq_cags.append(j)			

	list_of_cag_probs = [a_cags] + [b_cags] + [bw_cags] + [c_cags] + [dr_cags] + [drb345_cags] + [dq_cags]
	#print(list_of_cag_probs) 

	return list_of_cag_probs



def mapping_bws_for_gls(gl_string):
	'''Function to map Bw4 and Bw6 epitopes to B alleles in Genotype List String'''
	bws_list = []
	alleles_list = vxm_hla.gl_string_alleles_list(gl_string)

	for allele in alleles_list:
		bw = allele_to_ag_dict[allele][1]
		if bw != "NA":
			bws_list.append(bw)

	final_bws_mapped_to_gls = list(set(bws_list))
	
	return final_bws_mapped_to_gls		




def mapping_bws_for_macs(allele_codes_list):
	'''Function to map Bw4/Bw6 epitopes to B alleles in MACS'''
	bws_list = []
	alleles_list = vxm_hla.allele_code_to_allele_list(allele_codes_list)

	for allele in alleles_list:
		bw = allele_to_ag_dict[allele][1]
		if bw != "NA":
			bws_list.append(bw)

	final_bws_mapped_to_macs = list(set(bws_list))
	
	return final_bws_mapped_to_macs	



def group_unacceptable_antigens_per_locus_with_bw(sorted_allele_list,  ag_list):

	cag_prob_dict = {}
	
	## sorted allele dict is a list of sorted tuples 
	for tup in sorted_allele_list:
		allele = tup[0]
		ag = allele_to_ag_dict[allele][0]
		
		if allele in ag_list:
			cag_prob_dict[allele] = allele 

		if ag in ag_list:
			if allele in cag_prob_dict.keys():
				cag_prob_dict[allele] = cag_prob_dict[allele] + " " + "(" + ag + ")"

			else:
				cag_prob_dict[allele] = 	ag 

		if allele not in ag_list and ag not in ag_list:
			cag_prob_dict[allele] = ""
	
	
	bw_cags = []
	
	for i in ag_list:
		if i == "Bw4" or i == "Bw6":
			bw_prob_list_element = i 
			bw_cags.append(bw_prob_list_element)

	
	

	a_cags = []
	b_cags = []
	c_cags = []
	dr_cags = []
	drb345_cags = []
	dq_cags = []
	

	for i,j in cag_prob_dict.items():
		locus = i.split("*")[0]

		if locus == "A":
			a_cags.append(j)

		if locus == "B":
			b_cags.append(j)

		if locus == "C":
			c_cags.append(j)

		if locus == "DRB1":
			dr_cags.append(j)

		if locus == "DRB3" or locus == "DRB4" or locus == "DRB5":
			drb345_cags.append(j)	

		if locus == "DQB1":
			dq_cags.append(j)			

	list_of_cag_probs = [a_cags] + [b_cags] + [bw_cags] + [c_cags] + [dr_cags] + [drb345_cags] + [dq_cags]


	return list_of_cag_probs


##### Function to represent allele and antigen frequencies lower than 0.01 as scientific notation and rounding to 4 digits

def round_freq_se(frequency):
	'''Rounds the allele/ag frequencies or represents them in scientific notation depending on the length of floating point number'''
	if frequency < 0.01:
		edited_frequency = '%.2E' % Decimal(frequency)
	else: 
		edited_frequency = str(round(frequency, 4))

	return edited_frequency	
