from django.http import HttpResponse
from django.shortcuts import render
import operator
import re
import pandas as pd 
#from pd import *
import conversion_functions_for_VXM
from conversion_functions_for_VXM import gl_string_ags, genotype_ags, allele_code_ags, convert_allele_list_to_ags, convert_ag_list_to_gls

import virtual_crossmatch
from virtual_crossmatch import vxm_gls, vxm_allele_codes, vxm_uags, vxm_hIresalleles, UA_eq_dict, vxm_proposed_for_uags

import ancillary_funcs
from ancillary_funcs import group_serotypes_per_locus, split_gl_string_per_locus, prob_dict_list_of_strings, group_list_of_alleles_per_locus, prob_dict_list_of_strings_for_antigens
from ancillary_funcs import conflicts_ags, mapping_bws_for_gls, group_serotypes_per_locus_with_bw, group_unacceptable_antigens_per_locus_with_bw, mapping_bws_for_macs, group_allele_codes_per_locus
from ancillary_funcs import round_freq_se
import vxm_hla

############################################################To display full race for acronym selected ##############################################################

pop_acro_dict = {"AFA": "African American", "API": "Asia/Pacific Islander", "CAU": "Caucasian", "HIS": "Hispanic", 
"NAM": "Native American", "AAFA": "African American", "AFB": "African Black",  "AINDI": "South Asian Indian", 
"AISC": "American Indian-South or Central American", "ALANAM": "Alaska Native", "AMIND": "North American Indian", 
"CARB": "Caribbean Black", "CARHIS": "Caribbean Hispanic", "CARIBI": "Caribbean Indian", "EURCAU": "European Caucasian", 
"FILII": "Filipino", "HAWI": "Hawaiian or Pacific Islander", "JAPI": "Japanese", "KORI": "Korean", 
"MENAFC": "Middle Eastern or N. Coast of Africa", "MSWHIS": "Mexican or Chicano HIS", "NCHI": "Chinese", 
"SCAHIS": "South or Central AMerican Hispanic", "SCAMB": "South or Central American Black", 
"SCSEAI": "South East Asian", "VIET": "Vietnamese"}

#######################################################################################################################################################################

#########################################################################################################################################


def vxm_home(request):
	return render(request, 'victor/vxmHome.html')

#########################################################################################################################################

def victor_license(request):
    return render(request, 'victor/victor_license.html')

########################################################################################################################################

def unos_ags(request):
	return render(request, 'victor/uagsVxm.html')

########################################################################################################################################

def highres_allele(request):
	return render(request, 'victor/highRESallele.html')	

########################################################################################################################################

def gl_string(request):
    return render(request, 'victor/glsVxm.html')

########################################################################################################################################

def multiple_allele_codes(request):
	return render(request, 'victor/macVxm.html')

##########################################################################################################################################

def proposed_unos_ags(request):
	return render(request, 'victor/proposedUAGS.html')

#################################################################################################################################################################
def match_ags(request):
	output_dict = {}
	donorAgs = request.GET['userinput1'].strip()
	donorAgs = re.split(r'[;,\s]\s*' , donorAgs)
	#print(donorAgs)

	
	recepientAgs = request.GET['userinput2'].strip()
	if len(recepientAgs) == 0:
		recepientAgs = []
	else:
		recepientAgs = re.split(r'[;,\s]\s*' , recepientAgs)
	
	#print(recepientAgs)
	
	vxm_output = vxm_uags(donorAgs, recepientAgs)
	donorAntigens = sorted(list(set(vxm_output[0])))
	dags = ', '.join(donorAntigens)
	
	rags = ', '.join(sorted(vxm_output[1]))
	
	#print(vxm_output[2])
	conflicts = ', '.join(sorted(vxm_output[2]))
	
	

	##### all data required for table is pput in list 
	donor_typing = donorAgs
	recepient_ags = list(set(recepientAgs))
	conflicting_ags = vxm_output[2]
	#print(len(conflicting_ags))
	optn_equis = vxm_output[1]


	###### get list of locus from donor typing

	locus_list = []
	
	for ag in donor_typing:
		if ag[0:2].isalpha():
			locus = ag[0:2]
			locus_list.append(locus)
		else:
			locus = ag[0]
			locus_list.append(locus)	


	#final_locus_list = list(set(locus_list))

	final_locus_list = ["A", "B", "Bw", "C", "DR", "DRB3/4/5", "DQ"]	

		

	donor_ags_locus_list = group_serotypes_per_locus(donorAntigens)
	candidate_uags_locus_list = group_serotypes_per_locus(recepient_ags) 
	conflict_ags_locus_list = group_serotypes_per_locus(conflicting_ags)
	optne_locus_list = group_serotypes_per_locus(optn_equis)

	if len(conflicting_ags) == 0:
		end_result = "Virtual Crossmatch is Negative"
	else:
		end_result = "Virtual Crossmatch is Positive"	
	return render(request, 'victor/vxmUAGSmatch.html', {
		'donor_typing': dags, 'candidate_ags': rags, 'conflicted_ag': conflicts, "output3": end_result,  "zipped_list": zip(final_locus_list, donor_ags_locus_list, 
			candidate_uags_locus_list, optne_locus_list, conflict_ags_locus_list)})

###########################################################################################################################################################################



def match_hi_res_alleles(request):

	donor_typing = request.GET['userinput1'].strip()
	donor_alleles_list = re.split(r'[;,\s]\s*' , donor_typing)
	#print(donor_alleles_list)

	#donorUNOSAgs = conversion_functions_for_VXM.convert_allele_list_to_ags(donor_alleles_list)
	#print(donorUNOSAgs)
	recepientAgs = request.GET['userinput2'].strip()
	if len(recepientAgs) == 0:
		recepientAgs = []
	else:
		recepientAgs = re.split(r'[;,\s]\s*' , recepientAgs) 	

	vxm_output = vxm_hIresalleles(donor_alleles_list, recepientAgs)
	donorUNOSAgs= (vxm_output)[0]

	
	

	##############################################Locus list to be printed in the table rows ##############################################################################
	
	final_locus_list = ["A", "B", "Bw", "C", "DR", "DRB3/4/5", "DQ"]	
	final_alleles_list = group_list_of_alleles_per_locus(donor_alleles_list)	
	unos_eq_locus_list = group_serotypes_per_locus(donorUNOSAgs)
	ua_locus_list = group_serotypes_per_locus(recepientAgs)
	optne_locus_list = group_serotypes_per_locus(vxm_output[1])
	conflict_ags_locus_list = group_serotypes_per_locus(vxm_output[2])
	
	##################################################################################################################################################################################

	dags = ', '.join(sorted(vxm_output[0]))                            ##### String from list of donor antigens
	rags = ', '.join(sorted(vxm_output[1]))							   ##### String from list of recepient antigens

	conflicts = ', '.join(sorted(vxm_output[2]))

	if len(vxm_output[2]) == 0:
		end_result = "Virtual Crossmatch is Negative"
	else:
		end_result = "Virtual Crossmatch is Positive"	


	return render(request, 'victor/vxmAllelematch.html', {'donor_alleles': donor_typing,
		'donor_antigens': dags, 'candidate_ags': rags, 'conflicted_ag': conflicts, "output3": end_result, 
		'zipped_output': zip(final_locus_list, final_alleles_list, unos_eq_locus_list, ua_locus_list, optne_locus_list, conflict_ags_locus_list)})



#############################################################################################################################################

def match_gl(request):
	donorTyping = request.GET['userinput1']
	donorTyping = donorTyping.strip()
	popSpec = request.GET['userinput2']
	popSpecFul = pop_acro_dict[popSpec]
	recepientAntigens = request.GET['userinput3'].strip()
	pbTh = float(request.GET['userinput4'])

	

	if len(recepientAntigens) == 0:
		recepientAntigens = []
	else:
		recepientAntigens = re.split(r'[;,\s]\s*' , recepientAntigens)
	#print(recepientAntigens)
	
	entered_recepient_antigens = ", ".join(sorted(list(set(filter(None,recepientAntigens)))))
	#print(entered_recepient_antigens)
	vxm_output = vxm_gls(donorTyping, popSpec, recepientAntigens)
	donorAgs = sorted(vxm_output[0])
	donor_ags = ', '.join(donorAgs)
	candags = vxm_output[1]
	ag_probabilities = vxm_output[3]
	allele_probs = vxm_output[4]
	sorted_allele_probs = sorted(allele_probs.items(), key=operator.itemgetter(1), reverse=True)
	#print(sorted_allele_probs)
	antigen_probs = vxm_output[5]

	bw_prob = vxm_output[6]
	sorted_bw_probs = sorted(bw_prob.items(), key=operator.itemgetter(1), reverse=True)

	
	ua_locus_list = group_unacceptable_antigens_per_locus_with_bw(sorted_allele_probs, recepientAntigens)
	

	
	optne_locus_list = group_unacceptable_antigens_per_locus_with_bw(sorted_allele_probs, candags)
	recepient_ags = ', '.join(sorted(list(set(candags))))
	conflicted_ag = ', '.join(sorted(list(set(vxm_output[2]))))
	
	donor_bws = mapping_bws_for_gls(donorTyping)
	#print(donor_bws)
	donor_bws_string = "+ ".join(donor_bws)
	
	
	allele_list_with_probs = prob_dict_list_of_strings(sorted_allele_probs, sorted_bw_probs)
	
	antigen_list_with_probs = prob_dict_list_of_strings_for_antigens(sorted_allele_probs, antigen_probs)
	
	new_ag_probs = {}
	for ag, pp in ag_probabilities.items():
		if pp >= pbTh:
			new_ag_probs[ag] = pp
	#print(new_ag_probs)
	cags = []
	cag_probs = []
	for i, k in sorted(new_ag_probs.items()):
		ki = round_freq_se(k)
		cags.append(i)
		cag_probs.append(ki)
	
	#print(cags)
	#cag_list_above_th_locus_sorted = prob_dict_list_of_strings(new_ag_probs)
	cag_list_above_th_locus_sorted = conflicts_ags(sorted_allele_probs, new_ag_probs)
	#print(cag_list_above_th_locus_sorted)
	afterThcags = ", ".join(sorted(cags))
	final_locus_list = ["A", "B", "Bw", "C", "DR", "DRB3/4/5", "DQ"]	
	donor_strings = split_gl_string_per_locus(donorTyping, donor_bws_string)
	

	if len(cags) == 0:
		end_result = "Virtual Crossmatch is Negative"
	else:
		end_result = "Virtual Crossmatch is Positive"	
	

	gls_output_zipped_list = zip(final_locus_list, donor_strings, allele_list_with_probs, antigen_list_with_probs, ua_locus_list, optne_locus_list, cag_list_above_th_locus_sorted)
	
	return render(request, 'victor/vxmGlsmatch.html', {'donor_ags': donor_ags, 
		'entered_recepient_antigens': entered_recepient_antigens, 'recepient_ags': recepient_ags, 'output1': donorTyping, 'ethinicity': popSpecFul, 
		'output3': end_result, "conflicts": afterThcags, 'zipped_list': zip(cags, cag_probs), "output_zipped_list": 
		gls_output_zipped_list})


###################################################################################################################################################################

def match_proposed_uags(request):
	donorTyping = request.GET['userinput1'].strip()
	donorTyping = re.split(r'[;,\s]\s*' , donorTyping)
	#print(donorTyping)
	popSpec = request.GET['userinput2']
	popSpecFul = pop_acro_dict[popSpec]
	recepientAntigens = request.GET['userinput3']
	pbTh = float(request.GET['userinput4'])

	donor_gls = convert_ag_list_to_gls(donorTyping, popSpec)

	if len(recepientAntigens) == 0:
		recepientAntigens = []
	else:
		recepientAntigens = re.split(r'[;,\s]\s*' , recepientAntigens)
	
	#recepientAntigens = filter(None, recepientAntigens)
	
	entered_recepient_antigens = ", ".join(sorted(list(set(filter(None,recepientAntigens)))))
	#print(entered_recepient_antigens)
	vxm_output = vxm_proposed_for_uags(donorTyping, popSpec, recepientAntigens)
	donorAgs = sorted(vxm_output[0])
	donor_ags = ', '.join(donorAgs)
	candags = vxm_output[1]
	ag_probabilities = vxm_output[3]
	allele_probs = vxm_output[4]
	sorted_allele_probs = sorted(allele_probs.items(), key=operator.itemgetter(1), reverse=True)
	
	antigen_probs = vxm_output[5]
	bw_prob = vxm_output[6]
	sorted_bw_probs = sorted(bw_prob.items(), key=operator.itemgetter(1), reverse=True)
	
	ua_locus_list = group_unacceptable_antigens_per_locus_with_bw(sorted_allele_probs, recepientAntigens)
	
	
	
	
	
	optne_locus_list = group_unacceptable_antigens_per_locus_with_bw(sorted_allele_probs, candags)
	
	recepient_ags = ', '.join(sorted(list(set(candags))))
	conflicted_ag = ', '.join(sorted(list(set(vxm_output[2]))))
	
	donor_bws = mapping_bws_for_gls(donor_gls)
	#print(donor_bws)
	donor_bws_string = "+ ".join(donor_bws)
	
	
	allele_list_with_probs = prob_dict_list_of_strings(sorted_allele_probs, sorted_bw_probs)
	
	antigen_list_with_probs = prob_dict_list_of_strings_for_antigens(sorted_allele_probs, antigen_probs)
	
	new_ag_probs = {}
	for ag, pp in ag_probabilities.items():
		if pp >= pbTh:
			new_ag_probs[ag] = pp
	#print(new_ag_probs)
	cags = []
	cag_probs = []
	for i, k in sorted(new_ag_probs.items()):
		ki = round_freq_se(k)
		cags.append(i)
		cag_probs.append(ki)
	

	#cag_list_above_th_locus_sorted = prob_dict_list_of_strings(new_ag_probs)
	cag_list_above_th_locus_sorted = conflicts_ags(sorted_allele_probs, new_ag_probs)
	#print(cag_list_above_th_locus_sorted)
	afterThcags = ", ".join(sorted(cags))
	final_locus_list = ["A", "B", "Bw", "C", "DR", "DRB3/4/5", "DQ"]	
	donor_strings = split_gl_string_per_locus(donor_gls, donor_bws_string)
	

	if len(cags) == 0:
		end_result = "Virtual Crossmatch is Negative"
	else:
		end_result = "Virtual Crossmatch is Positive"	
	
	donor_ags_locus_list = group_serotypes_per_locus(donorTyping)
	gls_output_zipped_list = zip(final_locus_list, donor_ags_locus_list, allele_list_with_probs, antigen_list_with_probs, ua_locus_list, optne_locus_list, cag_list_above_th_locus_sorted)
	
	return render(request, 'victor/proposedUAGSmatch.html', {'donor_ags': donor_ags, 
		'entered_recepient_antigens': entered_recepient_antigens, 'recepient_ags': recepient_ags, 'output1': ", ".join(sorted(donorTyping)), 'ethinicity': popSpecFul, 
		'output3': end_result, "conflicts": afterThcags, 'zipped_list': zip(cags, cag_probs), "output_zipped_list": 
		gls_output_zipped_list})

	
###################################################################################################################################################################


def match_ac(request):
	donorTyping = request.GET['userinput1'].strip()
	donorCodes = re.split(r'[;,\s]\s*' , donorTyping)
	
	popSpec = request.GET['userinput2']
	
	popSpecFul = pop_acro_dict[popSpec]
	recepientAntigens = request.GET['userinput3'].strip()
	pbTh = float(request.GET['userinput4'])

	if len(recepientAntigens) == 0:
		recepientAntigens = []
	else:
		recepientAntigens = re.split(r'[;,\s]\s*' , recepientAntigens)
	
	entered_recepient_antigens = ", ".join(sorted(list(set(filter(None,recepientAntigens)))))

	vxm_output = vxm_allele_codes(donorCodes, popSpec, recepientAntigens)
	donorAgs = sorted(vxm_output[0])
	donor_ags = ', '.join(donorAgs)
	candags = vxm_output[1]
	'''recepient_ags = ', '.join(vxm_output[1])
	conflicted_ag = ', '.join(vxm_output[2])'''
	ag_probabilities = vxm_output[3]
	allele_probs = vxm_output[4]
	sorted_allele_probs = sorted(allele_probs.items(), key=operator.itemgetter(1), reverse=True)
	antigen_probs = vxm_output[5]
	bw_prob = vxm_output[6]
	sorted_bw_probs = sorted(bw_prob.items(), key=operator.itemgetter(1), reverse=True)

	
	ua_locus_list = group_unacceptable_antigens_per_locus_with_bw(sorted_allele_probs, recepientAntigens)
	print(ua_locus_list)


	optne_locus_list = group_unacceptable_antigens_per_locus_with_bw(sorted_allele_probs, candags)

	recepient_ags = ', '.join(sorted(list(set(candags))))
	conflicted_ag = ', '.join(sorted(list(set(vxm_output[2]))))

	#donor_bws = mapping_bws_for_gls(donorTyping)
	donor_bws = sorted(mapping_bws_for_macs(donorCodes))
	donor_bws_string = "+ ".join(donor_bws)

	new_ag_probs = {}
	for ag, pp in ag_probabilities.items():
		if pp >= pbTh:
			new_ag_probs[ag] = pp

	#print(new_ag_probs)
	cags = []
	cag_probs = []
	for i, k in sorted(new_ag_probs.items()):
		ki = round_freq_se(k)
		cags.append(i)
		cag_probs.append(ki)
	#print(cag_probs)

	#donor_strings = split_gl_string_per_locus(donorTyping, donor_bws_string)
	final_alleles_codes_list = group_allele_codes_per_locus(donorCodes, donor_bws_string)
	final_locus_list = ["A", "B", "Bw", "C", "DR", "DRB3/4/5", "DQ"]	
	
	allele_list_with_probs = prob_dict_list_of_strings(sorted_allele_probs, sorted_bw_probs)
	
	antigen_list_with_probs = prob_dict_list_of_strings_for_antigens(sorted_allele_probs, antigen_probs)
	

	
	#cag_list_above_th_locus_sorted = prob_dict_list_of_strings(new_ag_probs)
	cag_list_above_th_locus_sorted = conflicts_ags(sorted_allele_probs, new_ag_probs)

	afterThcags = ", ".join(sorted(cags))
	
	if len(cags) == 0:
		end_result = "Virtual Crossmatch is Negative"
	else:
		end_result = "Virtual Crossmatch is Positive"	
	#print(len(end_result))

	mac_output_zipped_list = zip(final_locus_list, final_alleles_codes_list, allele_list_with_probs, antigen_list_with_probs, ua_locus_list, optne_locus_list, cag_list_above_th_locus_sorted)
	return render(request, 'victor/vxmMACmatch.html', {'donor_ags': donor_ags, 'recepient_ags': recepient_ags, 
		'entered_recepient_antigens': entered_recepient_antigens, 'output1': donorTyping, 'ethinicity': popSpecFul, 
		'output3': end_result, "conflicts": afterThcags,  'zipped_list': zip(cags, cag_probs), "output_zipped_list": 
		mac_output_zipped_list})


#############################################################################################################################################################################	


def gl_string_extended_table(request):
    return render(request, 'victor/extendedTableGLmatch.html', {"output": match_gl[1][output_zipped_list]})
