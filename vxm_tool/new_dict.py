#!usr/env/bin/python
import operator
import decimal
from decimal import Decimal
import vxm_hla

allele_to_ag_dict = {}
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
		bw4_6 = row.split(',')[3]
	
	allele_to_ag_dict[allele_4d] = [antigen, bw4_6]


alleles_dict = {'A*02:01': 0.99951694966802, 'A*26:01': 0.99951694966802, 'A*02:55': 0.00048305033198005554, 
				'A*26:07': 0.00048305033198005554, 'C*05:01': 1.0, 'C*01:02': 1.0, 'B*15:03': 0.8535156141374103, 
				'B*27:05': 1.0000000000000002, 'B*15:01': 0.1408549708554045, 'B*15:02': 0.005040888092480535, 
				'B*15:04': 0.0005885269147048043, 
				'DRB1*12:01': 1.0, 'DRB1*01:01': 1.0, 'DQB1*03:01': 1.0, 'DQB1*05:01': 1.0}



			


sorted_dict = sorted(alleles_dict.items(), key=operator.itemgetter(1), reverse=True)

#print(sorted_dict)



def round_freq_se(frequency):
	if frequency < 0.01:
		edited_frequency = '%.2E' % Decimal(frequency)
	else: 
		edited_frequency = str(round(frequency, 4))

	return edited_frequency	







a_alleles = [] 
bw_epitopes = []
b_alleles = []
c_alleles = []
dr_alleles = []
dq_alleles = []


for tup in sorted_dict:
	i = tup[0]
	j = tup[1]


	if i.startswith("A"):
		ji = round_freq_se(j)
		fi = i + ": " + ji 
		a_alleles.append(fi)

	if i.startswith("B"):
		ji = round_freq_se(j)
		fi = i + ": " + ji 
		b_alleles.append(fi)


	if i.startswith("C"):
		ji = round_freq_se(j)
		fi = i + ": " + ji 
		c_alleles.append(fi)

	if i.startswith("DQ"):
		ji = round_freq_se(j)
		fi = i + ": " + ji 
		dq_alleles.append(fi)       

	if i.startswith("DR"):
		ji = round_freq_se(j)
		fi = i + ": " + ji 
		dr_alleles.append(fi)


	
list_of_allele_probs = [a_alleles] + [b_alleles]  + [c_alleles] + [dr_alleles] + [dq_alleles] 

#print(list_of_allele_probs)


conflicts = {'Bw4': 1, 'B27': 1, 'A2': 1.0, 'A*02:01': 0.99951694966802, 'A26': 1}




cag_prob_dict = {}

for tup in sorted_dict:
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
dq_cags = []
	

for i,j in cag_prob_dict.items():

	if i.startswith("A"):
		a_cags.append(j)

	if i.startswith("B"):
		b_cags.append(j)

	if i.startswith("C"):
		c_cags.append(j)

	if i.startswith("DR"):
		dr_cags.append(j)

	if i.startswith("DQ"):
		dq_cags.append(j)			

list_of_cag_probs = [a_cags] + [b_cags] + [bw_cags] + [c_cags] + [dr_cags] + [dq_cags]
#print(list_of_cag_probs) 

ag_probs = {'A2': 1.0, 'A26': 1.0, 'C05': 1.0, 'C01': 1.0, 'B72': 0.8535156141374103, 
			'B27': 1, 'Bw6': 1, 'Bw4': 1, 'B62': 0.1414434977701093, 'B75': 0.005040888092480535, 
			'DR12': 1.0, 'DR1': 1.0, 'DQ7': 1.0, 'DQ5': 1.0}


a_ags = [] 
b_ags = []
c_ags = []
dr_ags = []
dq_ags = []
bw_ags = []


for tup in sorted_dict:
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
	if ag_prob.startswith("A"):
		a_ags.append(ag_prob)

	if ag_prob.startswith("B"):
		b_ags.append(ag_prob)	


	if ag_prob.startswith("C"):
		c_ags.append(ag_prob)

	if ag_prob.startswith("DR"):
		dr_ags.append(ag_prob)

	if ag_prob.startswith("DQ"):
		dq_ags.append(ag_prob)				


list_of_ag_probs = [a_ags] + [b_ags] + [sorted(list(set(bw_ags)))] +[c_ags] + [dr_ags] + [dq_ags] 
#print(list_of_ag_probs)



ag_list = ['A2', 'Bw4', 'A*02:01']


cag_prob_dict = {}

for tup  in sorted_dict:
	i = tup[0]
	j = tup[1]
	ag = allele_to_ag_dict[i][0]
	if i in ag_list:
		cag_prob_dict[i] = i 

	if ag in ag_list:
		if i in cag_prob_dict.keys():
			cag_prob_dict[i] = cag_prob_dict[i] + " " + "(" + ag + ")"

		else:
			cag_prob_dict[i] = 	ag 

	else:
		cag_prob_dict[i] = ""
	
bw_cags = []
	
for i in ag_list:
	if i == "Bw4" or i == "Bw6":
		bw_prob_list_element = i 
		bw_cags.append(bw_prob_list_element)
	

a_cags = []
b_cags = []
c_cags = []
dr_cags = []
dq_cags = []
	

for i,j in cag_prob_dict.items():

	if i.startswith("A"):
		a_cags.append(j)

	if i.startswith("B"):
		b_cags.append(j)

	if i.startswith("C"):
		c_cags.append(j)

	if i.startswith("DR"):
		dr_cags.append(j)

	if i.startswith("DQ"):
		dq_cags.append(j)			

list_of_cag_probs = [a_cags] + [b_cags] + [bw_cags] + [c_cags] + [dr_cags] + [dq_cags]
print(list_of_cag_probs) 
