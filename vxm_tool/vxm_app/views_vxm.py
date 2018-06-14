from django.http import HttpResponse
from django.shortcuts import render

import re

import conversion_functions_for_VXM
from conversion_functions_for_VXM import gl_string_ags, genotype_ags, allele_code_ags

import virtual_crossmatch

from virtual_crossmatch import vxm_gls, vxm_allele_codes, vxm_uags, vxm_hIresalleles

pop_acro_dict = {"AFA": "African American", "API": "Asia/Pacific Islander", "CAU": "Caucasian", "HIS": "Hispanic", 
"NAM": "Native American", "AAFA": "African American", "AFB": "African Black",  "AINDI": "South Asian Indian", 
"AISC": "American Indian-South or Central American", "ALANAM": "Alaska Native", "AMIND": "North American Indian", 
"CARB": "Caribbean Black", "CARHIS": "Caribbean Hispanic", "CARIBI": "Caribbean Indian", "EURACU": "European Caucasian", 
"FILII": "Filipino", "HAWI": "Hawaiian or Pacific Islander", "JAPI": "Japanese", "KORI": "Korean", 
"MENAFC": "Middle Eastern or N. Coast of Africa", "MSWHIS": "Mexican or Chicano HIS", "NCHI": "Chinese", 
"SCAHIS": "South or Central AMerican Hispanic", "SCAMB": "South or Central American Black", 
"SCSEAI": "South East Asian", "VIET": "Vietnamese"}


#########################################################################################################################################


def vxm_home(request):
	return render(request, 'vxm_app/vxmHome.html')

#########################################################################################################################################

def victor_license(request):
    return render(request, 'vxm_app/victor_license.html')

########################################################################################################################################

def unos_ags(request):
	return render(request, 'vxm_app/uagsVxm.html')

########################################################################################################################################

def highres_allele(request):
	return render(request, 'vxm_app/highRESallele.html')	

########################################################################################################################################

def gl_string(request):
    return render(request, 'vxm_app/glsVxm.html')

########################################################################################################################################

def multiple_allele_codes(request):
	return render(request, 'vxm_app/macVxm.html')

##########################################################################################################################################

def match_ags(request):
	donorAgs = request.GET['userinput1'].strip()
	donorAgs = re.split(r'[;,\s]\s*' , donorAgs)

	
	recepientAgs = request.GET['userinput2'].strip()
	if len(recepientAgs) == 0:
		recepientAgs = []
	else:
		recepientAgs = re.split(r'[;,\s]\s*' , recepientAgs)
	

	vxm_output = vxm_uags(donorAgs, recepientAgs)
	dags = ', '.join(sorted(vxm_output[0]))
	rags = ', '.join(sorted(vxm_output[1]))
	
	conflicts = ', '.join(sorted(vxm_output[2]))
	
	if len(conflicts) == 0:
		end_result = "Virtual Crossmatch is Negative"
	else:
		end_result = "Virtual Crossmatch is Positive Because of Following Conflicting Antigens"	
	return render(request, 'vxm_app/vxmUAGSmatch.html', {
		'donor_typing': dags, 'candidate_ags': rags, 'conflicted_ag': conflicts, "output3": end_result})

############################################################################################################################################################



def match_hi_res_alleles(request):

	donor_typing = request.GET['userinput1'].strip()
	donor_alleles_list = re.split(r'[;,\s]\s*' , donor_typing)


	recepientAgs = request.GET['userinput2'].strip()
	if len(recepientAgs) == 0:
		recepientAgs = []
	else:
		recepientAgs = re.split(r'[;,\s]\s*' , recepientAgs) 	

	vxm_output = vxm_hIresalleles(donor_alleles_list, recepientAgs)

	dags = ', '.join(sorted(vxm_output[0]))
	rags = ', '.join(sorted(vxm_output[1]))

	conflicts = ', '.join(sorted(vxm_output[2]))


	if len(conflicts) == 0:
		end_result = "Virtual Crossmatch is Negative"
	else:
		end_result = "Virtual Crossmatch is Positive Because of Following Conflicting Antigens"	


	return render(request, 'vxm_app/vxmAllelematch.html', {'donor_alleles': donor_typing,
		'donor_antigens': dags, 'candidate_ags': rags, 'conflicted_ag': conflicts, "output3": end_result})























#############################################################################################################################################

def match_gl(request):
	donorTyping = request.GET['userinput1']
	donorTyping = donorTyping.strip()
	popSpec = request.GET['userinput2']
	popSpecFul = pop_acro_dict[popSpec]
	recepientAntigens = request.GET['userinput3']
	#print(type(recepientAntigens))

	if len(recepientAntigens) == 0:
		recepientAntigens = []
	else:
		recepientAntigens = re.split(r'[;,\s]\s*' , recepientAntigens)

	vxm_output = vxm_gls(donorTyping, popSpec, recepientAntigens)
	donor_ags = ', '.join(vxm_output[0])
	print(donor_ags)
	recepient_ags = ', '.join(vxm_output[1])
	conflicted_ag = ', '.join(vxm_output[2])
	ag_probabilities = vxm_output[3]
	print(ag_probabilities)
	print(conflicted_ag)
	cags = []
	cag_probs = []
	for i, k in ag_probabilities.items():
		cags.append(i)
		cag_probs.append(k)
	
	if len(conflicted_ag) == 0:
		end_result = "Virtual Crossmatch is Negative"
	else:
		end_result = "Virtual Crossmatch is Positive Because of Following Conflicting Antigens"	
	#print(len(end_result))
	return render(request, 'vxm_app/vxmGlsmatch.html', {'donor_ags': donor_ags, 
		'recepient_ags': recepient_ags, 'output1': donorTyping, 'ethinicity': popSpecFul, 
		'output3': end_result, "conflicts": conflicted_ag, "vxm_probs": ag_probabilities, 'zipped_list': zip(cags, cag_probs)})


###################################################################################################################################################################


def match_ac(request):
	donorTyping = request.GET['userinput1'].strip()
	donorCodes = re.split(r'[;,\s]\s*' , donorTyping)
	
	popSpec = request.GET['userinput2']
	
	popSpecFul = pop_acro_dict[popSpec]
	recepientAntigens = request.GET['userinput3'].strip()

	if len(recepientAntigens) == 0:
		recepientAntigens = []
	else:
		recepientAntigens = re.split(r'[;,\s]\s*' , recepientAntigens)
	
	

	vxm_output = vxm_allele_codes(donorCodes, popSpec, recepientAntigens)
	print(vxm_output)
	donor_ags = ', '.join(vxm_output[0])
	print(donor_ags)
	recepient_ags = ', '.join(vxm_output[1])
	conflicted_ag = ', '.join(vxm_output[2])
	print(conflicted_ag)
	ag_probabilities = vxm_output[3]
	print("antigen probabilities")
	print(ag_probabilities)
	cags = []
	cag_probs = []
	for i, k in ag_probabilities.items():
		cags.append(i)
		cag_probs.append(k)
	#print(cag_probs)
	if len(conflicted_ag) == 0:
		end_result = "Virtual Crossmatch is Negative"
	else:
		end_result = "Virtual Crossmatch is Positive Because of Following Conflicting Antigens"	
	#print(len(end_result))
	return render(request, 'vxm_app/vxmMACmatch.html', {'donor_ags': donor_ags, 
		'recepient_ags': recepient_ags, 'output1': donorTyping, 'ethinicity': popSpecFul, 
		'output3': end_result, "conflicts": conflicted_ag, "vxm_probs": ag_probabilities, 'zipped_list': zip(cags, cag_probs)})
