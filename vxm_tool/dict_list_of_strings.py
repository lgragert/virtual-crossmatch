#! usr/bin/python

import re

allele_probs = {'A*02:01': 0.0844, 'A*26:01': 0.0351, 'A*02:55': 0, 'A*26:07': 0.0001, 
			   'C*05:01': 0.0079, 'C*01:02': 0.11, 'B*15:01': 0.0312, 'B*15:02': 0.0422, 
			   'B*15:03': 0.0003, 'B*15:04': 0.0, 'B*27:05': 0.0082, 'DRB1*12:01': 0.0204, 
			   'DRB1*01:01': 0.0264, 'DQB1*03:01': 0.1849, 'DQB1*05:01': 0.0809}



a_alleles = []
b_alleles = []
c_alleles = []
dr_alleles = []
dq_alleles = []


for i,j in allele_probs.items():
        if re.findall("A", i):
            ji = str(j)
            fi = i + ": " + ji
            a_alleles.append(fi)

        if re.findall("B", i):
            ji = str(j)
            fi = i + ": " + ji 
            b_alleles.append(fi)


        if re.findall("C", i):
            ji = str(j)
            fi = i + ": " + ji 
            c_alleles.append(fi)


        if re.findall("DQ", i):
            ji = str(j)
            fi = i + ": " + ji 
            dq_alleles.append(fi)       

        if re.findall("DR", i):
            ji = str(j)
            fi = i + ": " + ji 
            dr_alleles.append(fi)



list_of_allele_probs = a_alleles + b_alleles + c_alleles + dr_alleles + dq_alleles 

print(list_of_allele_probs)
