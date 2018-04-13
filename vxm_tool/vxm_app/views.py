from django.http import HttpResponse
from django.shortcuts import render

import re

import conversion_functions
from conversion_functions import convert_allele_to_ag, convert_allele_list_to_ags, gl_string_ags, genotype_ags, allele_code_ags

import virtual_crossmatch

from virtual_crossmatch import vxm_gls, vxm_allele_codes

pop_acro_dict = {"AFA": "African American", "API": "Asia/Pacific Islander", "CAU": "Caucasian", "HIS": "Hispanic", 
"NAM": "Native American", "AAFA": "African American", "AFB": "African Black",  "AINDI": "South Asian Indian", 
"AISC": "American Indian-South or Central American", "ALANAM": "Alaska Native", "AMIND": "North American Indian", 
"CARB": "Caribbean Black", "CARHIS": "Caribbean Hispanic", "CARIBI": "Caribbean Indian", "EURACU": "European Caucasian", 
"FILII": "Filipino", "HAWI": "Hawaiian or Pacific Islander", "JAPI": "Japanese", "KORI": "Korean", 
"MENAFC": "Middle Eastern or N. Coast of Africa", "MSWHIS": "Mexican or Chicano HIS", "NCHI": "Chinese", 
"SCAHIS": "South or Central AMerican Hispanic", "SCAMB": "South or Central American Black", 
"SCSEAI": "South East Asian", "VIET": "Vietnamese"}


def vxm_home(request):
	return render(request, 'vxmHome.html')

def license(request):
    return render(request, 'license.html')

def gl_string(request):
    return render(request, 'glsVxm.html')


def multiple_allele_codes(request):
	return render(request, 'macVxm.html')

def match_gl(request):
	donorTyping = request.GET['userinput1']
	popSpec = request.GET['userinput2']
	popSpecFul = pop_acro_dict[popSpec]
	recepientAntigens = request.GET['userinput3']
	#print(type(recepientAntigens))
	recepientAntigens = re.split(r'[;,\s]\s*' , recepientAntigens)

	vxm_output = vxm_gls(donorTyping, popSpec, recepientAntigens)
	print(vxm_output)
	return render(request, 'vxmGlsmatch.html', {'donor_ags': vxm_output[0], 'recepient_ags': vxm_output[1], 'output1': donorTyping, 'ethinicity': popSpecFul, 'output3': recepientAntigens})




#def match_ac(request):
	#donorTyping = request.GET['userinput1']
	#popSpec = request.GET['userinput2']
	#recepientAntigens = request.GET['userinput3']
	#vxm_output = vxm_allele_codes(donorTyping, popSpec, recepientAntigens)
	#return render(request, vxmMac.html, {'output': vxm_output})
